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


def create_states():  # create a list of lists. Each list is a unique state of the system
    states = []
    for a in range(3):  # states of TG_1 (0, 1, 2)
        for b in range(3):  # states of TG_2 (0, 1, 2)
            for c in range(3):  # states of TC_1 (0, 1, 2)
                for d in range(3):  # states of TC_2 (0, 1, 2)
                    for e in range(2):  # states of EC (0, 1)
                        for f in range(2):  # states of TEG (0, 1)
                            state = [a, b, c, d, e, f]
                            states.append(state)
    return states


def calculate_rate_out():  # create a dictionary of the rates of transitioning out of the state of a single component
    TG_0 = TG_rates['TG_lambda01'] + TG_rates['TG_lambda02']
    TG_1 = TG_rates['TG_lambda12'] + TG_rates['TG_mu10']
    TG_2 = TG_rates['TG_mu20']
    TC_0 = TC_rates['TC_lambda01'] + TC_rates['TC_lambda02']
    TC_1 = TC_rates['TC_lambda12'] + TC_rates['TC_mu10']
    TC_2 = TC_rates['TC_mu20']
    EC_0 = EC_rates['EC_lambda']
    EC_2 = EC_rates['EC_mu']
    TEG_0 = TEG_rates['TEG_lambda']
    TEG_2 = TEG_rates['TEG_mu']
    out_rates = {'TG_0': TG_0, 'TG_1': TG_1, 'TG_2': TG_2, 'TC_0': TC_0, 'TC_1': TC_1,
                 'TC_2': TC_2, 'EC_0': EC_0, 'EC_2': EC_2, 'TEG_0': TEG_0, 'TEG_2': TEG_2}
    return out_rates


def time_to_transition(t_now, state, rates):  # calculate at what time the new transitions is
    TG_state_rates = [rates['TG_0'], rates['TG_1'], rates['TG_2']]
    TC_state_rates = [rates['TC_0'], rates['TC_1'], rates['TC_2']]
    EC_state_rates = [rates['EC_0'], 0, rates['EC_2']]
    TEG_state_rates = [rates['TEG_0'], 0, rates['TEG_2']]
    out_probability = TG_state_rates[state[0]] + TG_state_rates[state[1]] + TC_state_rates[state[2]] + \
                      TC_state_rates[state[3]] + EC_state_rates[state[4]] + TEG_state_rates[state[5]]

    # print("Probability of transitioning out of the system: ", out_probability)
    t_out = t_now - 1/out_probability*math.log(1-random.random())
    # print("Time at which transition will happen: ", t_out)

    return t_out, out_probability


def sample_new_state(state, rates, out_prob):
    TG_out_rates = [[[TG_rates['TG_lambda01'], 1], [TG_rates['TG_lambda02'], 2]],
                    [[TG_rates['TG_lambda12'], 2], [TG_rates['TG_mu10'], 0]],
                    [[TG_rates['TG_mu20'], 0]]]
    TC_out_rates = [[[TC_rates['TC_lambda01'], 1], [TC_rates['TC_lambda02'], 2]],
                    [[TC_rates['TC_lambda12'], 2], [TC_rates['TC_mu10'], 0]],
                    [[TC_rates['TC_mu20'], 0]]]
    EC_out_rates = [[[EC_rates['EC_lambda'], 2]], [], [[EC_rates['EC_mu'], 0]]]
    TEG_out_rates = [[[TEG_rates['TEG_lambda'], 2]], [], [[TEG_rates['TEG_mu'], 0]]]

    sample_r = random.random()  # random sample between 0.0 and 1.0
    sum_probabilities = 0
    new_state = state

    # print("State before change: ", new_state)
    for j in range(2):
        changing_state = TG_out_rates[j][0][1]
        for rate in TG_out_rates[state[j]]:
            if sum_probabilities > sample_r:
                break
            probability = rate[0] / out_prob
            sum_probabilities += probability
            # print("Checking state {}-> prob at {} with sample at {} and sum of probabilities {}".format(j, probability, sample_r, sum_probabilities))
            if sum_probabilities > sample_r:
                new_state[j] = changing_state
                break
            else:
                changing_state = rate[1]

    if sum_probabilities < sample_r:
        for j in range(2, 4):
            if sum_probabilities > sample_r:
                break
            changing_state = TC_out_rates[j-2][0][1]
            for rate in TC_out_rates[state[j]]:
                probability = rate[0] / out_prob
                sum_probabilities += probability
                # print("Checking state {}-> prob at {} with sample at {} and sum of probabilities {}".format(j, probability, sample_r, sum_probabilities))
                if sum_probabilities > sample_r:
                    if j == 2:
                        new_state[j-1] = changing_state
                    new_state[j] = changing_state
                    break
                else:
                    changing_state = rate[1]

    if sum_probabilities < sample_r:
        changing_state = EC_out_rates[0][0][1]
        for rate in EC_out_rates[state[4]]:
            if sum_probabilities > sample_r:
                break
            probability = rate[0] / out_prob
            sum_probabilities += probability
            # print("Checking state {}-> prob at {} with sample at {} and sum of probabilities {}".format(4, probability, sample_r, sum_probabilities))
            if sum_probabilities > sample_r:
                new_state[3] = changing_state
                break
            else:
                changing_state = rate[1]

    if sum_probabilities < sample_r:
        changing_state = TEG_out_rates[0][0][1]
        for rate in TEG_out_rates[state[5]]:
            if sum_probabilities > sample_r:
                break
            probability = rate[0] / out_prob
            sum_probabilities += probability
            # print("Checking state {}-> prob at {} with sample at {} and sum of probabilities {}".format(5, probability, sample_r, sum_probabilities))
            if sum_probabilities > sample_r:
                new_state[4] = changing_state
                break
            else:
                changing_state = rate[1]
                new_state[5] = changing_state

    return new_state


def run_MC():
    t = 0
    current_state = INITIAL_STATE
    while t < T_MISSION:
        t_current, prob_trans = time_to_transition(t, current_state, out_transition_rates)
        new_state = sample_new_state(current_state, out_transition_rates, prob_trans)
        t = t_current
        current_state = new_state
        # print("Current time: ", t)
        print("Current state: ", current_state)


all_states = create_states()  # get list of all possible states
# print(all_states)
print("Amount of states: ", len(all_states))
out_transition_rates = calculate_rate_out()
print("Out transition rates of all components for each of its states: ", out_transition_rates)
t_transition, prob_transition = time_to_transition(T_ZERO, INITIAL_STATE, out_transition_rates)
newest_state = sample_new_state(INITIAL_STATE, out_transition_rates, prob_transition)
print("Newest state: ", newest_state)
print("Running MC simulation...")
run_MC()
