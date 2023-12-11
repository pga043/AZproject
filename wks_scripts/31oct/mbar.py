import numpy as np
import pymbar 
from pymbar import *


T = 298.15 # kelvin
beta = 1/(0.001987204259*T)

# energies are in reduced units (kT)
matrix  = np.loadtxt('matrix.dat') 
vector = np.loadtxt('vector.dat')

mbar = MBAR(u_kn=matrix, N_k=vector)

result = mbar.compute_free_energy_differences()

# Free energy and uncertainity in kcal/mol
fe = result['Delta_f'][0,1]/beta
unc = result['dDelta_f'][0,1]/beta

print('kT = ' + str(1/beta))

print(fe, unc)

