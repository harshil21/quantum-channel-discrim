"""Distingusighing deploarizing from identity channel using picos"""

import picos as pic
import numpy as np
import matplotlib.pyplot as plt


def calc_lambda_from_p(p):
    return 1 - 4*p/3


def calc_p_from_lambda(lm):
    return 3*(1 - lm)/4


def cj_mat(p) -> np.array:
    lmVal = calc_lambda_from_p(p)
    cjMat = lambda lm : np.array([[(1. + lm)/2, 0., 0., lm],
                        [0., (1. - lm)/2, 0., 0.],
                        [0., 0., (1. - lm)/2, 0.],
                        [lm, 0., 0., (1 + lm)/2]])
    return cjMat(lmVal)


# Example 2. Using the dual formulation
def find_success_prob_distinguish_depolarizing_with_identity(p):
    upsBA = cj_mat(1)

    thetaBA = cj_mat(p)

    gamma2BA = upsBA - thetaBA  #Primal SDP
    #Constants
    #----------
    gamma2Pic = pic.Constant("gamma2BA", gamma2BA)
    
    shpBA = np.shape( gamma2Pic )
    shpB = np.shape(pic.partial_trace(gamma2Pic, subsystems=(1),dimensions=2))
    shpA = np.shape(pic.partial_trace(gamma2Pic, subsystems=(0),dimensions=2))
    iMatB = pic.Constant('Ib', np.eye(shpB[0]))
    
    #Variables
    #----------
    rhoPic = pic.HermitianVariable("rhoA", shpA)
    sigPic = pic.HermitianVariable("sigA", shpA)
    XPic = pic.ComplexVariable("X", shpBA)
    
    prob2P = pic.Problem()

    #Constraint
    #----------
    prob2P.add_constraint(((iMatB @ rhoPic & XPic ) // (XPic.H & iMatB @ sigPic)) >> 0)
    prob2P.add_constraint(pic.trace(rhoPic) == 1)
    prob2P.add_constraint(pic.trace(sigPic) == 1)
    
    #Objective
    #----------
    obj = pic.trace(gamma2Pic | XPic).real
    
    prob2P.set_objective('max',obj)

    #Solve the problem using mosek as a cvxopt
    prob2P.solve(verbosity=False,solver='cvxopt')
    #Solver claims to have found optimal saolution
    dNorm2P =  prob2P.value
    #Example 2 Dual Formulation
    #Constants
    #----------
    iMatA = pic.Constant('Ia', np.eye(shpA[0]))

    #Variables
    #----------
    NPicBA = pic.HermitianVariable("Nba", shpBA)
    MPicBA = pic.HermitianVariable("Mba", shpBA)
    mu = pic.RealVariable("mu")
    nu = pic.RealVariable("nu")
    prob2D = pic.Problem()

    # Constraint
    #----------
    prob2D.add_constraint(((NPicBA & -gamma2Pic) // (-gamma2Pic.H & MPicBA)) >> 0)
    
    NPicA = pic.partial_trace(NPicBA,subsystems=(0),dimensions=2)
    MPicA = pic.partial_trace(MPicBA,subsystems=(0),dimensions=2)
    
    prob2D.add_constraint(MPicA<<mu*iMatA)
    prob2D.add_constraint(NPicA<<nu*iMatA)

    # Objective
    #----------
    obj = (mu + nu)/2
    
    prob2D.set_objective('min',obj)
    # Solve the problem using mosek as a cvxopt
    prob2D.solve(verbosity=False,solver='cvxopt')
    # Solver claims to have found optimal solution
    dNorm2D =  prob2D.value
    # print('Diamond Norm distance between identity and equal probability Pauli error Channel')
    # print('Using Primal SDP = ', dNorm2P)
    # print('Using DualSDP = ', dNorm2D)
    # print('Difference between primal and dual values', abs(dNorm2D - dNorm2P))
    pE = (1 + dNorm2D/2)/2
    # print('Probability of distinguishing with an entangled input s* = ', pE)
    lmVal = 1.-4.*p/3.
    # print('Probability of distinguishing without an entangled input q* = ', (3 - lmVal)/4)
    return pE, (3 - lmVal)/4


def plot_example_2():
    entangled_probs = []
    non_entangled_probs = []
    p_vals = []

    for p in np.arange(0.0, 1.0, 0.05):
        p_vals.append(p)

    for p in p_vals:
        entangled_prob, non_entangled_prob = find_success_prob_distinguish_depolarizing_with_identity(p)
        entangled_probs.append(entangled_prob)
        non_entangled_probs.append(non_entangled_prob)

    print(entangled_probs)
    print(non_entangled_probs)
    print(len(p_vals))
    x = np.array(p_vals)
    y = np.array(non_entangled_probs)
    print(y)
    print(x.shape, y.shape)
    plt.plot(x, y)
    plt.show()


def run_example_2_specific_p(p):
    """p is calculated from lambda and has a range from 0 to 1"""
    print(f"For p = {p}:")
    entangled_prob, non_entangled_prob = find_success_prob_distinguish_depolarizing_with_identity(p)
    print(f"The entangled success probability is: {entangled_prob}")
    print(f"The non-entangled success probability is: {non_entangled_prob}")


def run_example_2_specific_depolarizing_factor(lmVal):
    """The depolarizing factor lambda should range from -0.333 to 1"""
    print(f"For lmVal = {lmVal}:")
    p = calc_p_from_lambda(lmVal)
    entangled_prob, non_entangled_prob = find_success_prob_distinguish_depolarizing_with_identity(p)
    print(f"The entangled success, and error probability is: {entangled_prob}, {1-entangled_prob}")
    print(f"The non-entangled success, and error probability is: {non_entangled_prob}, {1-non_entangled_prob}")


# run_example_2_specific_p(0.75)
run_example_2_specific_depolarizing_factor(0)
# plot_example_2()
