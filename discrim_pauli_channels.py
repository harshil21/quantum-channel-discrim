"""Discriminate between the three pauli channels"""

import numpy as np
import picos as pic
import random
import matplotlib.pyplot as plt
import matplotlib

# For displaying the plot in a separate window:
try:
    matplotlib.use("Qt5Agg")  # For Linux
except:
    matplotlib.use("TkAgg")  # For Windows / MacOS(?)


# Define 3 pauli matrices
pauli_x = np.array([[0, 1], [1, 0]])
pauli_y = np.array([[0, -1j], [1j, 0]])
pauli_z = np.array([[1, 0], [0, -1]])

# Define Identity matrix
identity = np.array([[1, 0], [0, 1]])

# Define |0> and |1>
ket_0 = np.array([[1], [0]])
ket_1 = np.array([[0], [1]])



# The commented out code below was used to generate random values for the qn lists. The only 
# difference this makes is that the error probability ranges will look different. However the plot
# will look the same (sinusoidal) and the optimal theta values will be the same.

# while True:
#     # Generate 4 random values for the first list
#     qn_1 = [random.random() for _ in range(4)]  # random.random() returns a random float [0, 1)
#     # Normalize the values so that the sum is 1
#     qn_1 /= np.sum(qn_1)
    
#     # Generate 4 random values for the second list
#     qn_2 = [random.random() for _ in range(4)]
#     qn_2 /= np.sum(qn_2)
    
#     # Check if the i-th elements are not the same
#     if all(qn_1[i] != qn_2[i] for i in range(4)):
#         break


def r(n):
    """Return r = p1*qn - p2*qn, where:

    p1: the probability of state passing through the first channel
    p2: the probability of state passing through the second channel
    qn: the probability value associated with the pauli matrix.

    Args:
        n (int): the index of the qn list. Corresponds to the index of the Pauli matrix. E.g.
            n=0 corresponds to identity, n=1 corresponds to pauli X, etc.
    """

    # The input state can pass through the channels equiprobably, i.e. p1=p2=1/2
    p1, p2 = 0.5, 0.5

    # The qn list should sum up to 1, e.g.
    qn_1 = [0.15, 0.25, 0.25, 0.35]
    qn_2 = [0.55, 0.2, 0.1, 0.15]
    assert np.sum(qn_1) == 1 == np.sum(qn_2)

    return p1*qn_1[n] - p2*qn_2[n]


def sigma(n):
    """Return the n-th Pauli matrix"""
    if n == 0:
        return identity
    elif n == 1:
        return pauli_x
    elif n == 2:
        return pauli_y
    elif n == 3:
        return pauli_z


def find_error_probability(theta):
    """Find the error probability of the three pauli channels, taking two at a time, i.e. 
    p1=p2=1/2. Uses a non-entangled input state.

    Args:
        theta (float): the angle of the input state in radians. 
            |psi> = cos(theta/2)|0> + sin(theta/2)|1>
    """

    # Define input state |psi>
    psi = pic.Constant(np.cos(theta / 2) * ket_0 + np.sin(theta / 2) * ket_1)

    inner = 0
    # Sum over the 4 pauli matrices
    for i in range(0, 4):
        inner += r(i) * sigma(i) * (psi * psi.H) * sigma(i)  # Apply the formula
    
    trace_norm = pic.NuclearNorm(inner)

    error_prob = 0.5 * (1 - trace_norm)
    return error_prob.value


def find_optimal_input_state():
    """Find the optimal input state that minimizes the error probability."""
    optimal_theta = 0
    lowest_error_p = 1
    for theta in np.arange(0, 2*np.pi, 0.05):
        # If the error probability is lower than the current lowest, update the lowest error
        if (error_p:=find_error_probability(theta)) < lowest_error_p:
            optimal_theta = theta
            lowest_error_p = error_p
    return lowest_error_p, optimal_theta


def plot_error_probability():
    """Plot the error probability of the three pauli channels, taking two at a time, i.e. 
    p1=p2=1/2. Uses a non-entangled input state.
    """

    theta = np.arange(0, 2*np.pi, 0.05)
    error_p = [find_error_probability(i) for i in theta]

    plt.plot(theta, error_p)
    # Change x axis to degrees
    plt.xticks(np.arange(0, 2*np.pi+0.1, np.pi/2), [f"{i:.0f}" for i in np.arange(0, 361, 90)])
    plt.title("Probabilites associated with discrimination of two Pauli channels\n(non-entangled)")
    plt.xlabel(f"\u03B8 (deg)")  # theta symbol
    plt.ylabel("Error Probability")
    plt.savefig("plots/error_probability_pauli_channel.png")
    plt.show()


# print(f"The lowest error probability and its theta value is: {find_optimal_input_state()}")
plot_error_probability()
