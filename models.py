import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as func
from torch.distributions.categorical import Categorical

class Actor(nn.Module):
    def __init__(self, inp_size, act_size, deterministic):
        super(Actor, self).__init__()
        assert inp_size == 2, 'only deal with input size of 2'
        assert act_size == 2, 'only deal with action size of 2'

        a1 = np.random.rand()
        a1 = torch.tensor(a1)

        a2 = np.random.rand()
        a2 = torch.tensor(a2)

        self.a1 = nn.Parameter(a1, requires_grad=True)
        self.a2 = nn.Parameter(a2, requires_grad=True)

        self.deterministic = deterministic

    def forward(self, obs):
                
        distance = (self.a1 * (obs[0][1]-0.5) + self.a2 * (obs[0][0]-0.5))
        probability_approach = 1 / (1 + torch.exp(-distance))

        if self.deterministic:
            if probability_approach > 0.5:
                a = 1
                logp = torch.log(probability_approach)
            else:
                a = 0
                logp = torch.log(1 - probability_approach)
        else:
            if torch.rand(1) < probability_approach:
                a = 1
                logp = torch.log(probability_approach)
            else:
                a = 0
                logp = torch.log(1 - probability_approach)

        entropy = -probability_approach * torch.log(probability_approach) - (1-probability_approach) * torch.log(1 - probability_approach)
        return int(a), logp, entropy

class Critic(nn.Module):
    def __init__(self, inp_size):
        super(Critic, self).__init__()
        self.l1 = nn.Linear(inp_size, 64)
        self.l2 = nn.Linear(64, 64)
        self.l3 = nn.Linear(64, 1)

    def forward(self, obs):
        x = func.relu(self.l1(obs))
        x = func.relu(self.l2(x))
        v = self.l3(x)
        return torch.squeeze(v, -1)

class Amemori(nn.Module):
    def __init__(self, inp_size, act_size, deterministic = False):
        super(Amemori, self).__init__()
        self.pi = Actor(inp_size, act_size, deterministic)
        self.v = Critic(inp_size)

####################

class Actor_with_parameter(nn.Module):
    def __init__(self, inp_size, act_size, deterministic):
        super(Actor_with_parameter, self).__init__()
        assert inp_size == 2, 'only deal with input size of 2'
        assert act_size == 2, 'only deal with action size of 2'

        a1 = np.random.rand()
        a1 = torch.tensor(a1)

        a2 = np.random.rand()
        a2 = torch.tensor(a2)

        a3 = np.random.rand()
        a3 = torch.tensor(a3)

        self.a1 = nn.Parameter(a1, requires_grad=True)
        self.a2 = nn.Parameter(a2, requires_grad=True)
        self.a3 = nn.Parameter(a3, requires_grad=True)

        self.deterministic = deterministic

    def forward(self, obs):
        
        distance = (self.a1 * (obs[0][1]-0.5) + self.a2 * (obs[0][0]-0.5))
        probability_approach = 1 / (1 + torch.exp(-self.a4*distance))

        if self.deterministic:
            if probability_approach > 0.5:
                a = 1
                logp = torch.log(probability_approach)
            else:
                a = 0
                logp = torch.log(1 - probability_approach)
        else:
            if torch.rand(1) < probability_approach:
                a = 1
                logp = torch.log(probability_approach)
            else:
                a = 0
                logp = torch.log(1 - probability_approach)

        entropy = -probability_approach * torch.log(probability_approach) - (1-probability_approach) * torch.log(1 - probability_approach)
        return int(a), logp, entropy

class VariationAmemori(nn.Module):
    def __init__(self, inp_size, act_size, deterministic = False):
        super(VariationAmemori, self).__init__()
        self.pi = Actor_with_parameter(inp_size, act_size, deterministic)
        self.v = Critic(inp_size)

####################

# Model parametrized with a1, a2, b1, b2

class ActorNaim(nn.Module):
    def __init__(self, inp_size, act_size, deterministic):
        super(ActorNaim, self).__init__()
        assert inp_size == 2, 'only deal with input size of 2'
        assert act_size == 2, 'only deal with action size of 2'

        a1 = np.random.rand()
        a1 = torch.tensor(a1)

        a2 = np.random.rand()
        a2 = torch.tensor(a2)

        self.a1 = nn.Parameter(a1, requires_grad=True)
        self.a2 = nn.Parameter(a2, requires_grad=True)

        b1 = np.random.rand()
        b1 = torch.tensor(b1)

        b2 = np.random.rand()
        b2 = torch.tensor(b2)

        self.b1 = nn.Parameter(b1, requires_grad=True)
        self.b2 = nn.Parameter(b2, requires_grad=True)

        self.deterministic = deterministic

    def forward(self, obs, microstimulation = None):
                
        distance = (self.a1 * (obs[0][1] + self.b1) + self.a2 * (obs[0][0] + self.b2))

        if microstimulation != None:
            distance = torch.tensor(microstimulation[0]) * self.a1 * (obs[0][1] + self.b1) + torch.tensor(microstimulation[1]) * self.a2 * (obs[0][0] + self.b2)

        probability_approach = 1 / (1 + torch.exp(-distance))

        if self.deterministic:
            if probability_approach > 0.5:
                a = 1
                logp = torch.log(probability_approach)
            else:
                a = 0
                logp = torch.log(1 - probability_approach)
        else:
            if torch.rand(1) < probability_approach:
                a = 1
                logp = torch.log(probability_approach)
            else:
                a = 0
                logp = torch.log(1 - probability_approach)

        entropy = -probability_approach * torch.log(probability_approach) - (1-probability_approach) * torch.log(1 - probability_approach)
        return int(a), logp, entropy

class Naim(nn.Module):
    def __init__(self, inp_size, act_size, deterministic = False):
        super(Naim, self).__init__()
        self.pi = ActorNaim(inp_size, act_size, deterministic)
        self.v = Critic(inp_size)

####################

# Model parametrized with a1, a2, a3

class ActorVariationNaim(nn.Module):
    def __init__(self, inp_size, act_size, deterministic):
        super(ActorVariationNaim, self).__init__()
        assert inp_size == 2, 'only deal with input size of 2'
        assert act_size == 2, 'only deal with action size of 2'

        a1 = np.random.rand()
        a1 = torch.tensor(a1)

        a2 = np.random.rand()
        a2 = torch.tensor(a2)

        self.a1 = nn.Parameter(a1, requires_grad=True)
        self.a2 = nn.Parameter(a2, requires_grad=True)

        a3 = np.random.rand()
        a3 = torch.tensor(a3)

        self.a3 = nn.Parameter(a3, requires_grad=True)

        self.deterministic = deterministic

    def forward(self, obs, microstimulation = None):
                
        distance = self.a1 * obs[0][1] + self.a2 * obs[0][0] + self.a3

        if microstimulation != None:
            distance = torch.tensor(microstimulation) * self.a1 * obs[0][1] + self.a2 * obs[0][0] + self.a3

        probability_approach = 1 / (1 + torch.exp(-distance))

        if self.deterministic:
            if probability_approach > 0.5:
                a = 1
                logp = torch.log(probability_approach)
            else:
                a = 0
                logp = torch.log(1 - probability_approach)
        else:
            if torch.rand(1) < probability_approach:
                a = 1
                logp = torch.log(probability_approach)
            else:
                a = 0
                logp = torch.log(1 - probability_approach)

        entropy = -probability_approach * torch.log(probability_approach) - (1-probability_approach) * torch.log(1 - probability_approach)
        return int(a), logp, entropy

class VariationNaim(nn.Module):
    def __init__(self, inp_size, act_size, deterministic = False):
        super(VariationNaim, self).__init__()
        self.pi = ActorVariationNaim(inp_size, act_size, deterministic)
        self.v = Critic(inp_size)
                
####################
        
class ClassicActor(nn.Module):
    def __init__(self, inp_size, act_size, deterministic):
        super(ClassicActor, self).__init__()
        self.l1 = nn.Linear(inp_size, 64)
        self.l2 = nn.Linear(64, 64)
        self.l3 = nn.Linear(64, act_size)

        self.deterministic = deterministic

    def forward(self, obs):
        x = func.relu(self.l1(obs))
        x = func.relu(self.l2(x))
        logits = self.l3(x)
        pi = Categorical(logits=logits)

        if self.deterministic == True:
            a = torch.argmax(pi.probs)
        else:
            a = pi.sample()
        
        logp = pi.log_prob(a) 
        entropy = pi.entropy()
        
        return int(a), logp, entropy

class A2C(nn.Module):
    def __init__(self, inp_size, act_size, deterministic = False):
        super(A2C, self).__init__()
        self.pi = ClassicActor(inp_size, act_size, deterministic)
        self.v = Critic(inp_size)

####################

class NotClassicActor(nn.Module):
    def __init__(self, inp_size, act_size, deterministic):
        super(NotClassicActor, self).__init__()
        self.l1 = nn.Linear(inp_size, 3)
        self.l2 = nn.Linear(3, 3)
        self.l3 = nn.Linear(3, act_size)

        a4 = np.random.rand()
        a4 = 0.5
        a4 = torch.tensor(a4)
        self.a4 = nn.Parameter(a4, requires_grad=True)

        self.deterministic = deterministic

    def forward(self, obs):
        x = func.relu(self.l1(obs))
        x = func.relu(self.l2(x))
        logits = self.l3(x)
        y1 = logits[0][0].clone()
        y2 = logits[0][1].clone()

        logits[0][0] = self.a4*y1
        logits[0][1] = self.a4*y2
    
        pi = Categorical(logits=logits)

        if self.deterministic == True:
            a = torch.argmax(pi.probs)
        else:
            a = pi.sample()
        
        logp = pi.log_prob(a) 
        entropy = pi.entropy()
        
        return int(a), logp, entropy

class VariationA2C(nn.Module):
    def __init__(self, inp_size, act_size, deterministic = False):
        super(VariationA2C, self).__init__()
        self.pi = NotClassicActor(inp_size, act_size, deterministic)
        self.v = Critic(inp_size)

####################

class RNNLSTM(nn.Module):
    def __init__(self, inp_size, out_size, hid_size):
        super(RNNLSTM, self).__init__()
        self.rnn = nn.LSTMCell(inp_size, hid_size)
        self.rnn_out = nn.Linear(hid_size, out_size)
        self.init = True

        self.h0 = nn.Parameter(torch.zeros((1, hid_size)),
                               requires_grad=True)
        self.c0 = nn.Parameter(torch.zeros((1, hid_size)),
                               requires_grad=True)
        self.h = None
        self.c = None

    def forward(self, inp):
        if self.init:
            hid_in = (self.h0, self.c0)
            self.init = False
        else:
            hid_in = (self.h, self.c)
        self.h, self.c = self.rnn(inp, hid_in)
        out = self.rnn_out(self.h)
        return out

    def reset(self):
        self.init = True
        self.h = None
        self.c = None

class RNNActor(nn.Module):
    def __init__(self, inp_size, act_size, hid_size):
        super(RNNActor, self).__init__()
        self.rnn = RNNLSTM(inp_size, act_size, hid_size)

    def _distribution(self, obs):
        logits = self.rnn(obs)
        return Categorical(logits=logits)

    def _log_prob_from_distribution(self, pi, act):
        return pi.log_prob(act)

    def forward(self, obs, act=None):
        pi = self._distribution(obs)
        logp_a = None
        if act is not None:
            logp_a = self._log_prob_from_distribution(pi, act)
        return pi, logp_a

    def reset(self):
        self.rnn.reset()

class RNNCritic(nn.Module):
    def __init__(self, inp_size, hid_size):
        super(RNNCritic, self).__init__()
        self.rnn = RNNLSTM(inp_size, 1, hid_size)

    def forward(self, obs):
        v = self.rnn(obs)
        return torch.squeeze(v, -1)

    def reset(self):
        self.rnn.reset()

class RNNActorCritic(nn.Module):
    def __init__(self, inp_size, act_size, hid_size):
        super(RNNActorCritic, self).__init__()
        self.pi = RNNActor(inp_size, act_size, hid_size)
        self.v = RNNCritic(inp_size, hid_size)

    def reset(self):
        self.pi.reset()
        self.v.reset()

    def step(self, obs):
        with torch.no_grad():
            pi = self.pi._distribution(obs)
            a = pi.sample()
            logp_a = self.pi._log_prob_from_distribution(pi, a)
            v = self.v(obs)
        return int(a), v.numpy(), logp_a.numpy()

    def act(self, obs):
        return self.step(obs)[0]