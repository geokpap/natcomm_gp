import os
import gym
import logging

import torch

from utils.models import *
from utils.environment import *
from utils.config_global import MAP_LOC, DEVICE

def make_input(obs):
    return torch.as_tensor(obs, dtype=torch.float32).unsqueeze(0).to(DEVICE)

def init_env(env_name):
    # initialize environment and get obs and act size
    if env_name == 'CostBenefit':
        environment = CostBenefitEnv()
        observation_dim = environment.observation_space.shape[0]
        action_dim = environment.action_space.n
    elif env_name == 'RewardNoReward':
        environment = RewardNoRewardEnv()
        observation_dim = environment.observation_space.shape[0]
        action_dim = environment.action_space.n
    else:
        raise NotImplementedError('Env not implemented')
    return environment, observation_dim, action_dim

def init_model(config, observation_dim, action_dim, mode):
    # initialize Actor Critic model
     
    print(config.model_type)
        
    if config.model_type == 'Amemori':
        ac_model = Amemori(observation_dim, action_dim)
    elif config.model_type == 'VariationAmemori':
        ac_model = VariationAmemori(observation_dim, action_dim)
    elif config.model_type == 'Naim':
        ac_model = Naim(observation_dim, action_dim)
    elif config.model_type == 'VariationNaim':
        ac_model = VariationNaim(observation_dim, action_dim)
    elif config.model_type == 'A2C':
        ac_model = A2C(observation_dim, action_dim)
    elif config.model_type == 'VariationA2C':
        ac_model = VariationA2C(observation_dim, action_dim)
    elif config.model_type == 'RNNActorCritic':
        ac_model = RNNActorCritic(observation_dim, action_dim, 100)
    else:
        raise NotImplementedError('AC not implemented')

    if mode == 'train':
        ac_model.train()
    elif mode == 'eval':
        model_path = os.path.join(config.save_path, 'net_best.pth')
        ac_model.load_state_dict(torch.load(model_path, map_location=MAP_LOC), strict=True)
        logging.info("Successfully loaded model: " + model_path)
        ac_model.eval()
    else:
        raise NotImplementedError('mode not specified')

    return ac_model.to(DEVICE)