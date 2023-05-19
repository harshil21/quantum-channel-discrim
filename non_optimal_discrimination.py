"""Distinguish depolarizing channel from identity using brute force"""
import picos as pic
import numpy as np
from itertools import product


p1 = 0.5
p2 = 1 - p1


def find_lowest_error_prob_qubit():
    d = 2
    iMat = pic.Constant('I', np.eye(d,dtype='complex'))
    ket_0 = np.array([[1.], [0.]])
    ket_1 = np.array([[0.], [1.]])
    error_probs = []
    alphas = []
    for i in np.arange(0, 1, 0.01):
        alpha = pic.Constant('alpha', i, shape=(1, 1))
        alphas.append(alpha)
        beta = (1 - alpha**2) ** (1/2)
        beta = pic.Constant('beta', beta, shape=(1, 1))
        ket_0 = pic.Constant('ket_0', ket_0, shape=(2, 1))
        ket_1 = pic.Constant('ket_1', ket_1, shape=(2, 1))
        psi = (alpha*ket_0) + (beta*ket_1)

        inner = p1*psi*psi.H - p2*(iMat/d)
        trace_norm = pic.NuclearNorm(inner)
        error_p = 1/2*(1-trace_norm)
        error_probs.append(error_p.value)

    # print(error_probs)
    print(f"The lowest error probability is {min(error_probs)}")
    print(f"The maximum error prob is: {max(error_probs)}")
    print(f"The corresponding min alpha value is {alphas[error_probs.index(min(error_probs))]}")
    print(f"The corresponding max alpha value is: {alphas[error_probs.index(max(error_probs))]}")


def find_lowest_error_prob_qutrit():
    d = 3
    iMat = pic.Constant('I', np.eye(d,dtype='complex'))
    ket_0 = pic.Constant('ket_0', np.array([[1.], [0.], [0.]]), shape=(3, 1))
    ket_1 = pic.Constant('ket_1', np.array([[0.], [1.], [0.]]), shape=(3, 1))
    ket_2 = pic.Constant('ket_2', np.array([[0.], [0.], [1.]]), shape=(3, 1))
    error_probs = []
    alphas = betas = [i for i in np.arange(0, 1.05, 0.05)]
    for i in product(alphas, betas):  # Takes the cartesian product of the two lists
        alpha = i[0]
        beta = i[1]
        if alpha**2 + beta**2 > 1:  # Solution will not converge if this is the case
            continue
        # print(alpha, beta)
        gamma = (1 - alpha**2 - beta**2) ** (0.5)
        alpha = pic.Constant('alpha', alpha, shape=(1,1))
        beta = pic.Constant('beta', beta, shape=(1,1))
        gamma = pic.Constant('gamma', gamma, shape=(1,1))
        psi = (alpha*ket_0) + (beta*ket_1) + (gamma*ket_2)

        inner = p1*psi*psi.H - p2*(iMat/d)
        trace_norm = pic.NuclearNorm(inner)
        error_p = 1/2*(1-trace_norm)
        error_probs.append((error_p.value, i[0], i[1], gamma))
    
    error_probs.sort(key=lambda x: x[0])  # Sorts the probabilities in ascending order
    print(f"The lowest error probability is {error_probs[0][0]}")
    print(f"The maximum error prob is: {error_probs[-1][0]}")
    print(f"The corresponding min alpha, beta and gamma value is [{error_probs[0][1]} {error_probs[0][2]} {error_probs[-1][3]}]")
    print(f"The corresponding max alpha, beta and gamma value is [{error_probs[-1][1]} {error_probs[-1][2]} {error_probs[-1][3]}]")


# find_lowest_error_prob_qutrit()


def find_lowest_error_prob_qubit_depolarizing_factor(lmVal):
    """This function calculates the error probability for distinguishing a depolarizing channel 
    from an identity channel when a qubit is passed through them. It is also able to factor in the
    depolarizing factor to find the probability."""
    d = 2
    iMat = pic.Constant('I', np.eye(d, dtype='complex'))
    error_probs = []
    alphas = np.arange(0, 1, 0.01)
    ket_0 = pic.Constant('ket_0',  np.array([[1.], [0.]]), shape=(2, 1))
    ket_1 = pic.Constant('ket_1', np.array([[0.], [1.]]), shape=(2, 1))
    for alpha in alphas:
        beta = np.sqrt(1 - alpha**2)
        alpha = pic.Constant('alpha', alpha, shape=(1, 1))
        beta = pic.Constant('beta', beta, shape=(1, 1))
        psi = (alpha*ket_0) + (beta*ket_1)
        print(psi)
        inner = p1*psi*psi.H - p2*(iMat/d)
        print(psi*psi.H)
        trace_norm = pic.NuclearNorm(inner)
        error_p = (1 - lmVal)*0.25 + lmVal*trace_norm
        print(error_p.value)
        print()
        error_probs.append(error_p.value)

    min_error_prob = min(error_probs)
    max_error_prob = max(error_probs)
    min_alpha = alphas[error_probs.index(min_error_prob)]
    max_alpha = alphas[error_probs.index(max_error_prob)]

    # print(error_probs)
    print(f"The lowest error probability is {min_error_prob}")
    print(f"The maximum error probability is {max_error_prob}")
    print(f"The corresponding minimum alpha value is {min_alpha}")
    print(f"The corresponding maximum alpha value is {max_alpha}")

# find_lowest_error_prob_qubit()  # Calculates probability for lambda=0
print()
find_lowest_error_prob_qubit_depolarizing_factor(0.1)
