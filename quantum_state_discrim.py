
# obj = 

# Example 1
# sig1 = np.array([[1.,0],[0.,0.]])
# sig2 = np.array([[0.,0],[0.,1.]])

# qVal = np.random.rand()

# Sgs1 = pic.Constant("sg1", sig1)
# Sgs2 = pic.Constant("sg2", sig2)

# q = pic.Constant("q", qVal)

# #Identity matrix
# shp = np.shape(sig1)

# #Variables
# #----------
# eVars = pic.HermitianVariable("E1", shp)

# prob1 = pic.Problem()
    
# #Constraint
# #----------
# prob1.add_constraint(eVars >> 0)
# prob1.add_constraint(iMat - eVars >> 0)

# #Objective
# #----------
# obj = q*(Sgs1 | eVars) + (1-q)*(Sgs2 | iMat - eVars)

# prob1.set_objective('max',obj)
# print(prob1)
# prob1.solve(verbosity=False,solver='cvxopt')
# pvm1 = np.array(eVars.value)
# prb = prob1.value
# print('Probability of discriminating the inputs is = ', prb)
# print(eVars)

# # Example 2

# d = 6

# mt1 = np.random.rand(d,d) + 1j*np.random.randn(d,d)
# mt1 = np.dot(mt1,mt1.conj().T)
# sig1 = mt1/np.trace(mt1)

# mt2 = np.random.rand(d,d) + 1j*np.random.randn(d,d)
# mt2 = np.dot(mt2,mt2.conj().T)
# sig2 = mt2/np.trace(mt2)

# qVal = 0.0
# #Constants
# #----------
# Sgs1 = pic.Constant("sg1", sig1)
# Sgs2 = pic.Constant("sg2", sig2)

# q = pic.Constant("q", qVal)

# #Identity matrix
# shp = np.shape(sig1)
# iMat = pic.Constant('I', np.eye(shp[0],dtype='complex'))

# #Variables
# #----------
# eVars = pic.HermitianVariable("E1", shp)

# prob2 = pic.Problem()
    
# #Constraint
# #----------
# prob2.add_constraint(eVars >> 0)
# prob2.add_constraint(iMat - eVars >> 0)

# #Objective
# #----------
# obj = q*(Sgs1 | eVars) + (1-q)*(Sgs2 | iMat - eVars)

# prob2.set_objective('max',obj.real)
# print(prob2)
# prob2.solve(verbosity=False,solver='cvxopt')
# #Helstrom result computes using numpy
# delOpt = qVal*sig1 - (1-qVal)*sig2
# # print(np.linalg.eigh(delOpt))
# (u,s,vH) = np.linalg.svd(delOpt, hermitian=False)
# pStarAlg = (np.sum(s) + 1)/2
# prb2 = prob2.value
# print('Probability of discriminating the input')
# print('Using SDP = ', prb2)
# print('Helstrom result', pStarAlg)
# print('Difference between these values', abs(prb2 - pStarAlg))
# print("Sig1 and Sig2 is:")
# print(sig1.shape)
# print(sig1)
# print()
# print(sig2)


#Example 3

# sig1 = np.array([[1.,0.],[0.,0.]])
# sig2 = np.array([[1.0, 1.0],[1.,1.]])/2
# sig3 = np.array([[1.,0.],[0.,1.]])/2

# qVal = 0.0
# #Constants
# #----------
# Sgs1 = pic.Constant("sg1", sig1)
# Sgs2 = pic.Constant("sg2", sig2)
# Sgs3 = pic.Constant("sg3", sig3)

# q = pic.Constant("q", qVal)

# #Identity matrix
# shp = np.shape(sig1)
# iMat = pic.Constant('I', np.eye(shp[0],dtype='complex'))

# #Variables
# #----------
# eVars1 = pic.HermitianVariable("E1", shp)
# eVars2 = pic.HermitianVariable("E2", shp)

# prob3 = pic.Problem()
    
# #Constraint
# #----------
# prob3.add_constraint(eVars1 >> 0)
# prob3.add_constraint(eVars2 >> 0)
# prob3.add_constraint(iMat - eVars1 -eVars2 >> 0)

# #Objective
# #----------
# obj3 = (q/2)*(Sgs1 | eVars1) + (q/2)*(Sgs2 | eVars2) + (1-q)*(Sgs3 | iMat - eVars1 - eVars2)

# prob3.set_objective('max',obj3)
# #User readable view of the problem being composed in PICOS
# print(prob3)
# prob3.solve(verbosity=False,solver='cvxopt')
# pov1 = np.array(eVars1.value)
# pov2 = np.array(eVars2.value)
# prb3 = prob3.value
# print('Probability of discriminating the input given by SDP = ', prb3)
# print('Algebraic solution = ', max(qVal,1-qVal))
# print(eVars1)
# print()
# print(eVars2)