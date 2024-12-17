import os.path as osp

class BaseConfig(object):
    def __init__(self):
        
        self.env = 'CostBenefit'

        self.experiment_name = 'test'
        self.save_path = osp.join('output', self.experiment_name)

        self.epochs = 1500
        self.save_every = 50
        self.log_every = 50
        self.optimizer_type = 'Adam'

        #self.model_type = 'Amemori'
        #self.model_type = 'VariationAmemori'
        self.model_type = 'Naim'
        #self.model_type = 'VariationNaim'
        #self.model_type = 'A2C'
        #self.model_type = 'VariationA2C'

        self.count_reversal = 100
        self.reversals = [4000]
        
        self.steps_per_epoch = 100
        self.max_tra_len = 100

        self.a1_lr = 1e-2
        self.a2_lr = 1e-2 / 2
        self.a3_lr = 1e-2
        self.a4_lr = 1e-2
        self.b1_lr = 1e-2
        self.b2_lr = 1e-2
        
        self.pi_lr = 1e-2
        self.vf_lr = 1e-2
        self.beta = 0.2

        # RNN
        self.gamma = 0.99
        self.lam = 0.97
        self.train_pi_iters = 40
        self.train_v_iters = 40
        self.clip_ratio = 0.2
        self.target_kl = 0.01

        self.deterministic_evaluation = False

    def update(self, new_config):
        self.__dict__.update(new_config.__dict__)

    def __str__(self):
        return str(self.__dict__)