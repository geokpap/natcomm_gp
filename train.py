import os
import shutil
import numpy as np
import torch
import scipy.signal
import logging
import argparse
import pickle

from utils.config_global import LOG_LEVEL, USE_CUDA, DEVICE
from utils.logger import Logger
from utils.config import BaseConfig
from utils.train_utils import init_model, init_env

def discount_cumsum(x, discount):
    """
    magic from rllab for computing discounted cumulative sums of vectors.

    input:
        vector x,
        [x0,
         x1,
         x2]

    output:
        [x0 + discount * x1 + discount^2 * x2,
         x1 + discount * x2,
         x2]
    """
    return scipy.signal.lfilter([1], [1, float(-discount)], x[::-1], axis=0)[::-1]


class VPGBuffer:
    """
    A buffer for storing trajectories experienced by a VPG agent interacting
    with the environment, and using Generalized Advantage Estimation (GAE-Lambda)
    for calculating the advantages of state-action pairs.
    """

    def __init__(self, size, gamma=0.99, lam=0.95):
        self.rew_buf = np.zeros(size, dtype=np.float32)
        self.val_buf = np.zeros(size, dtype=np.float32)

        self.adv_buf = np.zeros(size, dtype=np.float32)
        self.ret_buf = np.zeros(size, dtype=np.float32)

        self.gamma, self.lam = gamma, lam
        self.ptr, self.path_start_idx, self.max_size = 0, 0, size

    def store(self, rew, val):
        """
        Append one timestep of agent-environment interaction to the buffer.
        """
        assert self.ptr < self.max_size  # buffer has to have room so you can store
        self.rew_buf[self.ptr] = rew
        self.val_buf[self.ptr] = val
        self.ptr += 1

    def finish_path(self, last_val=0):
        """
        Call this at the end of a trajectory, or when one gets cut off
        by an epoch ending. This looks back in the buffer to where the
        trajectory started, and uses rewards and value estimates from
        the whole trajectory to compute advantage estimates with GAE-Lambda,
        as well as compute the rewards-to-go for each state, to use as
        the targets for the value function.

        The "last_val" argument should be 0 if the trajectory ended
        because the agent reached a terminal state (died), and otherwise
        should be V(s_T), the value function estimated for the last state.
        This allows us to bootstrap the reward-to-go calculation to account
        for timesteps beyond the arbitrary episode horizon (or epoch cutoff).
        """

        path_slice = slice(self.path_start_idx, self.ptr)
        rews = np.append(self.rew_buf[path_slice], last_val)
        vals = np.append(self.val_buf[path_slice], last_val)

        # the next two lines implement GAE-Lambda advantage calculation
        deltas = rews[:-1] + self.gamma * vals[1:] - vals[:-1]
        self.adv_buf[path_slice] = discount_cumsum(deltas, self.gamma * self.lam)

        # the next line computes rewards-to-go, to be targets for the value function
        self.ret_buf[path_slice] = discount_cumsum(rews, self.gamma)[:-1]

        self.path_start_idx = self.ptr

    def get(self):
        """
        Call this at the end of an epoch to get advantages and reward_togo from
        the buffer, with advantages appropriately normalized (shifted to have
        mean zero and std one). Also, resets some pointers in the buffer.
        """
        assert self.ptr == self.max_size    # buffer has to be full before you can get
        self.ptr, self.path_start_idx = 0, 0
        # the next two lines implement the advantage normalization trick
        adv_mean = np.mean(self.adv_buf)
        adv_std = np.std(self.adv_buf)
        self.adv_buf = (self.adv_buf - adv_mean) / adv_std
        return (torch.as_tensor(self.ret_buf, dtype=torch.float32).to(DEVICE),
                torch.as_tensor(self.adv_buf, dtype=torch.float32).to(DEVICE))

if __name__ == '__main__':
    if USE_CUDA:
        logging.info("training with GPU")
        
    # general configuration
    cfg = BaseConfig()
    
    # initialize logger
    logger = Logger(output_dir=cfg.save_path,
                    exp_name=cfg.experiment_name)

    env, obs_dim, act_dim = init_env(cfg.env) 
    ac = init_model(cfg, obs_dim, act_dim, mode='train')
    
    if cfg.optimizer_type == 'Adam':
        
        if cfg.model_type == 'Amemori':
            pi_optimizer = torch.optim.Adam([
                {'params': ac.pi.a1, 'lr': cfg.a1_lr},
                {'params': ac.pi.a2, 'lr': cfg.a2_lr}
            ])
        elif cfg.model_type == 'VariationAmemori':
            pi_optimizer = torch.optim.Adam([
                {'params': ac.pi.a1, 'lr': cfg.a1_lr},
                {'params': ac.pi.a2, 'lr': cfg.a2_lr},
                {'params': ac.pi.a4, 'lr': cfg.a4_lr}
            ])
        elif cfg.model_type == 'Naim':
            pi_optimizer = torch.optim.Adam([
                {'params': ac.pi.a1, 'lr': cfg.a1_lr},
                {'params': ac.pi.a2, 'lr': cfg.a2_lr},
                {'params': ac.pi.b1, 'lr': cfg.b1_lr},
                {'params': ac.pi.b2, 'lr': cfg.b2_lr}
            ])
        elif cfg.model_type == 'VariationNaim':
            pi_optimizer = torch.optim.Adam([
                {'params': ac.pi.a1, 'lr': cfg.a1_lr},
                {'params': ac.pi.a2, 'lr': cfg.a2_lr},
                {'params': ac.pi.a3, 'lr': cfg.a3_lr},
            ])
        else:
            pi_optimizer = torch.optim.Adam(ac.pi.parameters(), lr = cfg.pi_lr)
        # print(pi_optimizer)
        vf_optimizer = torch.optim.Adam(ac.v.parameters(), lr = cfg.vf_lr)
    else:
        raise NotImplementedError('optimizer not implemented')
    
    # Set up experience buffer
    buf = VPGBuffer(cfg.steps_per_epoch, cfg.gamma, cfg.lam)

    log_rewards = []
    tra_rewards = []
    tra_lens = []
    i_log = 0
    o, tra_rew, tra_len = env.reset(), 0, 0

    if os.path.exists(cfg.save_path):
        shutil.rmtree(cfg.save_path)
        os.makedirs(cfg.save_path)

    check_lr = 1
    count_rewards = 0
    apav_results = {}

    for ep in range(cfg.epochs):
        # save model
        if ep % cfg.save_every == 0:

            torch.save(ac.state_dict(),
                       os.path.join(cfg.save_path,
                                    'net_ep{}.pth'.format(ep)))
        logps = []
        values = []
        entropies = []

        # if ep in cfg.reversals:
        #     print('REVERSE')
        #     env.reverseTask()

        if count_rewards > cfg.count_reversal or ep in cfg.reversals:
            print(ep)
            print('REVERSE')
            env.reverseTask()
            # pi_optimizer = torch.optim.Adam([
            #     {'params': ac.pi.a1, 'lr': cfg.a1_lr*check_lr},
            #     {'params': ac.pi.a2, 'lr': cfg.a2_lr*check_lr}
            # ])
            # check_lr += .5
            # print(pi_optimizer)

        count_rewards = 0

        for t_ in range(cfg.steps_per_epoch):
            o = torch.as_tensor(o, dtype=torch.float32).unsqueeze(0).to(DEVICE)
            a, logp, entropy = ac.pi(o)
            v = ac.v(o)

            o, r, d, _ = env.step(a)
            tra_rew += r
            tra_len += 1

            logps.append(logp)
            values.append(v)
            entropies.append(entropy)
            buf.store(r, v.item())

            timeout = tra_len == cfg.max_tra_len
            terminal = d or timeout
            epoch_ended = t_ == (cfg.steps_per_epoch - 1)
            if terminal or epoch_ended:
                if epoch_ended and not terminal:
                    logging.debug('trajectory cut off by epoch at %d steps.' % tra_len)
                # if trajectory didn't reach terminal state, bootstrap value target
                if timeout or epoch_ended:
                    with torch.no_grad():
                        v = ac.v(torch.as_tensor(o, dtype=torch.float32).unsqueeze(0).to(DEVICE))
                    v = v.item()
                else:
                    v = 0.0
                buf.finish_path(v)
                if terminal:
                    # only save tra_rew / tra_len if trajectory finished
                    tra_rewards.append(tra_rew)
                    tra_lens.append(tra_len)
                o, tra_rew, tra_len = env.reset(), 0, 0

        # log performance
        if ep % cfg.log_every == 0:
            i_log += 1
            log_rew = np.mean(tra_rewards)
            count_rewards = log_rew
            log_rewards.append(log_rew)

            logger.log_tabular('Epoch', ep)
            logger.log_tabular('EnvInteracts', (ep+1)*cfg.steps_per_epoch)
            logger.log_tabular('AvgTraReward', log_rew)
            logger.log_tabular('TraLength', np.mean(tra_lens))
            logger.dump_tabular()
            tra_rewards = []
            tra_lens = []

            if False:
                # save the model with best testing loss
                if log_rew >= max(log_rewards):
                    torch.save(ac.state_dict(),
                            os.path.join(cfg.save_path, 'net_best.pth'))

        # get advantages, reward_to_gos
        reward_to_gos, advantages = buf.get()

        logps = torch.squeeze(torch.stack(logps, dim = 0))
        values = torch.squeeze(torch.stack(values, dim = 0))
        entropies = torch.squeeze(torch.stack(entropies, dim = 0))
        
        # Actor learning
        pi_optimizer.zero_grad()
        loss_pi = - (logps * advantages).mean() - cfg.beta * (entropies).mean()
        loss_pi.backward()
        pi_optimizer.step()

        # Critic learning
        vf_optimizer.zero_grad()
        loss_v = ((values - reward_to_gos) ** 2).mean()
        loss_v.backward()
        vf_optimizer.step()

        apav_results[str(ep)] = count_rewards

    # a_file = open(str(cfg.count_reversal) + ".pkl", "wb")
    # pickle.dump(apav_results, a_file)
    # a_file.close()