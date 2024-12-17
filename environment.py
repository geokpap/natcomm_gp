import numpy as np
import gym

#### AP - AV

class CostBenefitEnv(gym.Env):
    def __init__(self):
        super().__init__()
        self.previous_ob = np.random.rand(2)
        self.observation_space = gym.spaces.Box(low=np.float32(np.zeros(2)), high=np.float32(np.ones(2)))
        self.action_space = gym.spaces.Discrete(2)
        # self.action_space = gym.spaces.Box(low=np.float32(np.zeros(1)), high=np.float32(np.ones(1)))
        self.reverse = False

    def step(self, action):

        ob = self.previous_ob
        done = False
        info = {}
        reward = 0

        if self.reverse:
            cost = ob[1]
            benefit = ob[0]
        else:
            cost = ob[0]
            benefit = ob[1]

        if action == 1:
            reward = benefit - cost

        ob = np.random.rand(2)  # New observation
        self.previous_ob = ob

        return ob, reward, done, info

    def reset(self):
        # PAVLOVIAN CASES
        self.previous_ob = np.random.rand(2)
        return self.previous_ob
    
    def reverseTask(self):
        self.reverse = not self.reverse
        
#### AP - AP - Reward, no reward

class RewardNoRewardEnv(gym.Env):
    def __init__(self):
        super().__init__()
        self.previous_ob = np.random.rand(2)
        self.observation_space = gym.spaces.Box(low=np.float32(np.zeros(2)), high=np.float32(np.ones(2)))
        self.action_space = gym.spaces.Discrete(2)
        # self.action_space = gym.spaces.Box(low=np.float32(np.zeros(1)), high=np.float32(np.ones(1)))
        self.reverse = False

    def step(self, action):

        ob = self.previous_ob
        done = False
        info = {}
        reward = 0

        if self.reverse:
            cost = ob[1]
            benefit = ob[0]
            if action == 1:
                reward = benefit
        else:
            cost = ob[0]
            benefit = ob[1]
            if action == 0:
                reward = benefit

        ob = np.random.rand(2)  # New observation
        self.previous_ob = ob

        return ob, reward, done, info

    def reset(self):
        self.previous_ob = np.random.rand(2)
        return self.previous_ob
    
    def reverseTask(self):
        self.reverse = not self.reverse