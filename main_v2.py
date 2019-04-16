import math
import random

T_MISSION = 10000
N_SIMS = 100
T_ZERO = 0
INITIAL_STATE = [0, 0, 0, 0, 0, 0]

# counters

TG_rates = {'TG_lambda01': 0.79e-3, 'TG_lambda02': 0.77e-3, 'TG_lambda12': 1.86e-3, 'TG_mu10': 0.032, 'TG_mu20': 0.038}
TC_rates = {'TC_lambda01': 0.67e-3, 'TC_lambda02': 0.74e-3, 'TC_lambda12': 2.12e-3, 'TC_mu10': 0.003, 'TC_mu20': 0.048}
EC_rates = {'EC_lambda': 0.17e-3, 'EC_mu': 0.032}
TEG_rates = {'TEG_lambda': 5.7e-5, 'TEG_mu': 0.333}

CHANGE_DICT = {0: {0: [[TG_rates['TG_lambda01'], 1], [TG_rates['TG_lambda02'], 2]], 1: [[TG_rates['TG_lambda12'], 2], [TG_rates['TG_mu10'], 0]], 2: [[TG_rates['TG_mu20'], 0]]}, 
               1: {0: [[TG_rates['TG_lambda01'], 1], [TG_rates['TG_lambda02'], 2]], 1: [[TG_rates['TG_lambda12'], 2], [TG_rates['TG_mu10'], 0]], 2: [[TG_rates['TG_mu20'], 0]]}, 
               2: {0: [[TC_rates['TC_lambda01'], 1], [TC_rates['TC_lambda02'], 2]], 1: [[TC_rates['TC_lambda12'], 2], [TC_rates['TC_mu10'], 0]], 2: [[TC_rates['TC_mu20'], 0]]}, 
               3: {0: [[TC_rates['TC_lambda01'], 1], [TC_rates['TC_lambda02'], 2]], 1: [[TC_rates['TC_lambda12'], 2], [TC_rates['TC_mu10'], 0]], 2: [[TC_rates['TC_mu20'], 0]]}, 
               4: {0: [[EC_rates['EC_lambda'], 2]], 2: [[EC_rates['EC_mu'], 0]]},
               5: {0: [[TEG_rates['TEG_lambda'], 2]], 2: [[TEG_rates['TEG_mu'], 0]]}}


def calculate_rate_out():  # create a dictionary of the rates of transitioning out of the state of a single component
    tg_0 = TG_rates['TG_lambda01'] + TG_rates['TG_lambda02']
    tg_1 = TG_rates['TG_lambda12'] + TG_rates['TG_mu10']
    tg_2 = TG_rates['TG_mu20']
    tc_0 = TC_rates['TC_lambda01'] + TC_rates['TC_lambda02']
    tc_1 = TC_rates['TC_lambda12'] + TC_rates['TC_mu10']
    tc_2 = TC_rates['TC_mu20']
    ec_0 = EC_rates['EC_lambda']
    ec_2 = EC_rates['EC_mu']
    teg_0 = TEG_rates['TEG_lambda']
    teg_2 = TEG_rates['TEG_mu']
    out_rates = {'TG_0': tg_0, 'TG_1': tg_1, 'TG_2': tg_2, 'TC_0': tc_0, 'TC_1': tc_1,
                 'TC_2': tc_2, 'EC_0': ec_0, 'EC_2': ec_2, 'TEG_0': teg_0, 'TEG_2': teg_2}
    return out_rates


def time_to_transition(t_now, state, rates):  # calculate at what time the new transitions is
    tg_state_rates = [rates['TG_0'], rates['TG_1'], rates['TG_2']]
    tc_state_rates = [rates['TC_0'], rates['TC_1'], rates['TC_2']]
    ec_state_rates = [rates['EC_0'], 0, rates['EC_2']]
    teg_state_rates = [rates['TEG_0'], 0, rates['TEG_2']]
    out_probability = tg_state_rates[state[0]] + tg_state_rates[state[1]] + tc_state_rates[state[2]] + \
                      tc_state_rates[state[3]] + ec_state_rates[state[4]] + teg_state_rates[state[5]]

    # print("Probability of transitioning out of the system: ", out_probability)
    t_out = t_now - 1/out_probability*math.log(1-random.random())
    # print("Time at which transition will happen: ", t_out)

    return t_out, out_probability


def sample_new_state(state, out_probability):
    component = 0
    for i in range(len(state)):
        
    
    return new_state