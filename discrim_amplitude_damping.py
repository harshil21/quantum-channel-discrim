"""
Perform one-shot discrimination of amplitude damping channels, from the paper: 
Discriminating qubit amplitude damping channels by M.Rexiti, Stefano Mancini
"""

import itertools
import numpy as np
import picos as pic
import matplotlib
import matplotlib.pyplot as plt


def gen_ket(n):
    """Generate the ket |n>"""
    return pic.Constant(np.array([[0]] * n + [[1]] + [[0]] * (4 - n - 1)))


def gen_bra(n):
    """Generate the bra <n|"""
    return gen_ket(n).H


def ladder_dag_a(n):
    """Ladder operator a, it increases the input state by 1.
    
    Args:
        n (int): the input state. E.g. |0> -> |1>, |1> -> |2>, etc.
    """
    return pic.Constant(np.sqrt(n+1) * gen_ket(n+1))


ladder_dag_b = ladder_dag_a


def ladder_a(n):
    """Ladder operator a, it decreases the input state by 1.
    
    Args:
        n (int): the input state. E.g. |1> -> |0>, |2> -> |1>, etc.
    """
    if n == 0:
        return pic.Constant(np.array([[0]] * 4))
    return pic.Constant(np.sqrt(n-1) * gen_ket(n-1))

ladder_b = ladder_a

original_bases = [
    (0, 0),
    (0, 1),
    (1, 0),
    (1, 1)
]


class AmplitudeChannelDiscrim:
    def __init__(self, x, eta_0, eta_1):
        """Initialize the class.

        `eta_0` and `eta_1` are in the range [0, pi/2]. `eta_0` must be greater than `eta_1`.
        
        Args:
            x (float): amplitude damping parameter, should be betweeen 0 and 1.
            eta_0 (float): the decay probability of the first channel. Must be greater than eta_1
            eta_1 (float): the decay probability of the second channel. Must be less than eta_0
        
        """
        self.x = x
        self.eta_0 = eta_0
        self.eta_1 = eta_1
        # assert self.eta_0 > self.eta_1 or self.eta_1 > self.eta_0  # without loss of generality
        self.gamma = np.cos(self.eta_1) + np.cos(self.eta_0)

    # The hamiltonian of the system
    def hamilitonian(self, eta) -> np.array:
        """Calculate the hamiltonian of the system, in the basis |00>, |01>, |10>, |11>"""
        H_mat = np.zeros((4, 4))

        # H = eta * (dag_a*b + a*dag_b)
        # Lets calculate <00|H|00>, <00|H|01>, and so on
        # This is done by expanding the H inside. So we get e.g.
        # <01|H|10> = (<0|dag_a|1> * <1|b|0>) + (<0|a|1> * <1|dag_b|0>)
        # the itertools.product is used to get all combinations of the bases (to generate the matrix)
        for bases in itertools.product(original_bases, original_bases):
            # the 0 and 1 correspond to the positions. E.g. <00| and |01> -> (0,0) and (0, 1)
            (bra_0, bra_1), (ket_0, ket_1) = bases
            inner = ((gen_bra(bra_0) * ladder_dag_a(ket_0)) * (gen_bra(bra_1) * ladder_b(ket_1))) \
                    + (gen_bra(bra_0) * ladder_a(ket_0)) * (gen_bra(bra_1) * ladder_dag_b(ket_1))
            # Substituting the values in the matrix:
            H_mat[bra_0 * 2 + bra_1, ket_0 * 2 + ket_1] = eta * inner
        
        assert H_mat[2, 1] == eta == H_mat[1, 2]  # Confirm we got the correct H matrix
        return H_mat

    def unitary(self, hamiltonian, eta) -> pic.Constant:
        """Generate the unitary matrix from the Hamiltonian using the formula: U = e^(-iH)"""

        # The e^(-iH) maclaurin series expansion is:
        # e^(-iH) = I - iH - H^2/2! + iH^3/3! + H^4/4! - iH^5/5! + ...
        I = np.eye(4)
        H = hamiltonian
        U = I - (1j * H) - ((H @ H) / 2) + ((1j * H @ H @ H) / 6) + ((H @ H @ H @ H) / 24) - ((1j * H @ H @ H @ H @ H) / 120)

        # Now let's confirm that this U matrix is the same as the one in the paper with sin and cos:
        final_U = np.array(
            [
                [1, 0, 0, 0],
                [0, np.cos(eta), 1j*-np.sin(eta), 0],
                [0, 1j*-np.sin(eta), np.cos(eta), 0],
                [0, 0, 0, 1]
            ]
        )
        # Confirm that the two matrices are indeed the same
        assert np.allclose(U, final_U, atol=0.1)  # atol is the absolute tolerance
        return pic.Constant(U)

    def input_state(self) -> pic.Constant:
        """The input state of the amplitude damping channel is given by:
        
        |psi> = sqrt(1-x)|0> + exp(-i*phi)*sqrt(x) |1>,
        where phi = 0, x = [0, 1], eta_0 > eta_1, thus it simplifies to:

        |psi> = sqrt(1-x)|0> + sqrt(x) |1>
        """
        psi = pic.Constant(
            [
                [(1-self.x) ** 0.5],
                [(self.x) ** 0.5],
            ],
            shape=(2, 1)
        )
        return psi

    def density_matrix(self, input_state, unitary_mat, eta) -> pic.Constant:
        """Calculate the density matrix of the system, by tracing out the environment."""
        # To calculate: Tr2[U(eta_i)*|psi>|0><psi|<0|U(eta_i)^dag]
        ket_0 = pic.Constant(np.array([[1], [0]]), shape=(2, 1))
        psi = input_state
        part_1 = unitary_mat * (psi @ ket_0)
        part_2 = (psi.H @ ket_0.H) * unitary_mat.H
        total = part_1 * part_2
        traced_out = total.partial_trace(1)

        # Now let's check if what we calculated is correct by checking the hand solved matrix
        rho = np.array(
            [
                [1-self.x + np.sin(eta)**2, np.sqrt((1-self.x) * self.x) * np.cos(eta)],
                [np.sqrt(self.x * (1-self.x)) * np.cos(eta), self.x * np.cos(eta)**2],
            ]
        )
        assert np.allclose(traced_out, rho, atol=1.0)  # Thus, the two matrices are the same
        return pic.Constant(traced_out, shape=(2, 2))

    def calc_success_prob_from_density_mat(self, input_density_mat_1, input_density_mat_2) -> float:
        """We can calculate the probability of discriminating two channels via the formula:
        
        P_success = 1/2 * (1 + 1/2(||rho_1 - rho_2||_1)
        """
        difference = input_density_mat_1 - input_density_mat_2
        trace_norm = pic.NuclearNorm(difference)
        P_success = 0.5 * (1 + 0.5 * trace_norm)
        return P_success
    
    def calc_success_prob_from_etas_and_x(self) -> float:
        """We can also calculate the success probability of discriminating channels by the formula:
        
        P_success = 1/2 * {1 + (cos(eta_1) - cos(eta_0)) * sqrt(x * [1-x*(1-gamma^2)]) },
        where gamma = cos(eta_1) + cos(eta_0)

        Ref: see formula (13) in the paper.
        """
        larger = max(self.eta_0, self.eta_1)
        smaller = min(self.eta_0, self.eta_1)
        P_success = 0.5 * (1 + (np.cos(smaller) - np.cos(larger)) * np.sqrt(self.x * (1 - self.x * (1 - self.gamma**2))))
        return P_success

    def generate_density_matrices_for_different_eta(self) -> tuple[pic.Constant, pic.Constant]:
        """Return the density matrices for two different eta values, keeping x constant"""
        h1 = self.hamilitonian(self.eta_0)
        h2 = self.hamilitonian(self.eta_1)
        u1 = self.unitary(h1, self.eta_0)
        u2 = self.unitary(h2, self.eta_1)
        psi = self.input_state()
        rho1 = self.density_matrix(psi, u1, self.eta_0)
        rho2 = self.density_matrix(psi, u2, self.eta_1)
        return rho1, rho2


def find_success_prob_specific_values(x, eta_0, eta_1):
    amp_channel = AmplitudeChannelDiscrim(x=x, eta_0=eta_0, eta_1=eta_1)
    rho1, rho2 = amp_channel.generate_density_matrices_for_different_eta()
    print(f"Success probability from density matrices: {amp_channel.calc_success_prob_from_density_mat(rho1, rho2)}")
    print(f"Success probability from x, eta_0 & eta_1: {amp_channel.calc_success_prob_from_etas_and_x()}")
    print(f"Error prob: {1 - amp_channel.calc_success_prob_from_etas_and_x()}")


def find_optimal_value_of_x():
    # Case 1:
    # Let's take the case of gamma >= 1/sqrt(2), we should always get x=1 as optimal:
    eta_0, eta_1 = 0.2, 0.1 # Gives gamma of 1.97, which is greater than 1/sqrt(2) = 0.707
    gamma = np.cos(eta_0) + np.cos(eta_1)
    assert gamma >= 1/np.sqrt(2)
    x_vals = np.arange(0, 1.1, 0.1)
    probs = []
    for x in x_vals:
        amp_channel = AmplitudeChannelDiscrim(x=x, eta_0=eta_0, eta_1=eta_1)
        rho1, rho2 = amp_channel.generate_density_matrices_for_different_eta()
        # print(x)
        # print(f"Success probability from density matrices: {amp_channel.calc_success_prob_from_density_mat(rho1, rho2)}")
        # print(f"Success probability from x, eta_0 & eta_1: {amp_channel.calc_success_prob_from_etas_and_x()}")
        # print()
        probs.append(amp_channel.calc_success_prob_from_density_mat(rho1, rho2))
    print(f"Given gamma = {gamma}, we should always get x=1 as optimal")
    print(f"Optimal value of x for gamma >= 1/sqrt(2): {x_vals[np.argmax(probs)]}")
    print(f"The corresponding maximum success probability is: {np.max(probs)}")

    print()
    # Case 2:
    # Let's take the case of gamma < 1/sqrt(2), optimal x = 1/(2(1-gamma^2)):
    eta_0, eta_1 = 1.5, 1.0 # Gives gamma of 0.62, which is less than 1/sqrt(2) = 0.707
    gamma = np.cos(eta_0) + np.cos(eta_1)
    assert gamma < 1/np.sqrt(2)
    x_vals = np.arange(0, 1.1, 0.1)
    probs = []
    for x in x_vals:
        amp_channel = AmplitudeChannelDiscrim(x=x, eta_0=eta_0, eta_1=eta_1)
        rho1, rho2 = amp_channel.generate_density_matrices_for_different_eta()
        probs.append(amp_channel.calc_success_prob_from_density_mat(rho1, rho2))

    print(f"Given gamma = {gamma}:")
    print(f"Optimal value of x for gamma < 1/sqrt(2): {x_vals[np.argmax(probs)]}")
    print(f"The corresponding maximum success probability is: {np.max(probs)}")
    print(f"The theoretical optimal value of x for gamma < 1/sqrt(2): {np.round(1/(2*(1-gamma**2)), decimals=1)}")


def contour_plot_vary_etas():
    range_of_etas = np.arange(0, 1.6, 0.1)

    for x in (0.25, 0.5, 0.75, 1):
        error_probs = []
        for eta in itertools.product(range_of_etas, range_of_etas):
            eta_0, eta_1 = eta
            amp_channel = AmplitudeChannelDiscrim(x=x, eta_0=eta_0, eta_1=eta_1)
            success_prob = amp_channel.calc_success_prob_from_etas_and_x()
            error_prob = 1 - success_prob
            error_probs.append(error_prob)

        error_probs = np.array([error_probs])
        error_probs = error_probs.reshape(len(range_of_etas), len(range_of_etas))
        eta_list = list(range_of_etas)
        cs = plt.contourf(eta_list, eta_list, error_probs)
        l = [f'{a:.2f}' for a in cs.levels]
        proxy = [plt.Rectangle((1, 1), 2, 2, fc=pc.get_facecolor()[0]) for pc in
                cs.collections]
        plt.legend(proxy, l)
        plt.xlabel("n0")
        plt.ylabel("n1")
        plt.title(f"Error probabilities for {x=}")
        plt.savefig(f"contour_{x=}.png")


def contour_plot_vary_x_and_eta_1():
    range_of_eta_1 = np.linspace(1.4, 3.1415/2, 30)
    x_vals = np.linspace(0, 1, 30)

    eta_0 = 1.1
    error_probs = []
    for x, eta_1 in itertools.product(x_vals, range_of_eta_1):
        amp_channel = AmplitudeChannelDiscrim(x=x, eta_0=eta_0, eta_1=eta_1)
        success_prob = amp_channel.calc_success_prob_from_etas_and_x()
        error_prob = 1 - success_prob
        error_probs.append(error_prob)

    error_probs = np.array([error_probs])
    error_probs = error_probs.reshape(len(range_of_eta_1), len(x_vals))
    x_list = list(x_vals)
    eta_1_list = list(range_of_eta_1)
    cs = plt.contourf(x_list, eta_1_list, error_probs.T)
    l = [f'{a:.2f}' for a in cs.levels]
    proxy = [plt.Rectangle((1, 1), 2, 2, fc=pc.get_facecolor()[0]) for pc in
            cs.collections]
    plt.legend(proxy, l)
    plt.xlabel("x")
    plt.ylabel("n1")
    plt.title(f"Error probabilities for {eta_0=}")
    plt.savefig(f"Error probs with {eta_0=}.png")

# find_success_prob_specific_values(x=0.3, eta_0=1.1, eta_1=1.40)
# contour_plot_vary_etas()
# contour_plot_vary_x_and_eta_1()
find_optimal_value_of_x()
