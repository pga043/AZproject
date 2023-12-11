import numpy as np
import pymbar 
from pymbar import *
import netCDF4 as nc


# read simulation .nc file (all the information is stored here)
file =  nc.Dataset('simulation.nc')


#--------- ipython stuff -------------
#file
#file.variables
#file.variables['energies']
#file.variables['energies'][0:,0:,0].shape
#file.variables['energies'][0:,0:,0].T.shape

#------------------------------------

nstates = file.variables['energies'].shape[2]
nframes = file.variables['energies'].shape[0]

T = 298.15 # kelvin
beta = 1/(0.001987204259*T)


# Initialize an empty matrix
matrix = np.empty((nstates,0))

'''
units: kT
long_name: energies[iteration][replica][state] is the reduced (unitless) energy of replica 'replica' from iteration 'iteration' evaluated at the thermodynamic state 'state'.
'''
for i in range(nstates):
    energies = file.variables['energies'][0:,0:,i].T 
    matrix = np.concatenate((matrix, energies), axis=1)


vector = np.empty(nstates)
vector.fill(nframes)

mbar = MBAR(u_kn=matrix, N_k=vector)

result = mbar.compute_free_energy_differences()

# Free energy and uncertainity in kcal/mol
fe = result['Delta_f'][0,1]/beta
unc = result['dDelta_f'][0,1]/beta

print('kT = ' + str(1/beta))

print(fe, unc)

print('Free energy difference is {:.3f} +- {:.3f} kT'.format(result['Delta_f'][0,1], result['dDelta_f'][0,1]))

print('Free energy difference is {:.3f} +- {:.3f} kcal/mol'.format(result['Delta_f'][0,1]/beta, result['dDelta_f'][0,1]/beta))


