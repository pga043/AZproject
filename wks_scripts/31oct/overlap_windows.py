import numpy as np
import sys
import matplotlib.pyplot as plt

m = np.genfromtxt('matrix.dat')

# number of windows and number of frames
n  = 11
fr = 5001


'''
matrix format:

lam1  U(from win1 frames @lam1) U(from win2 frames @lam1) U(from win3 frames @lam1) ........ U(from winN frames @lam1)
lam2  U(from win1 frames @lam2) U(from win2 frames @lam2) U(from win3 frames @lam2) ........ U(from winN frames @lam2)
.
.
.
lamN U(from win1 frames @lamN) U(from win2 frames @lamN) U(from win3 frames @lamN) ........ U(from winN frames @lamN)


Diagonal elements belong to the actual simulations that were run.
Non-diagonal elements are obtained via post processing of digonal trajectories for each lambda.

'''

# difference always final-initial

#------- test -----------

#win_ij = m[0, 5001:10001] - m[0, 0:5000]

#win_ji = m[1, 0:5000] - m[1, 5001:10001]

#------------------------

for i in range(n-1):
    
    win_ij = m[i, (i+1)*fr:((i+1)*fr+fr-1)] - m[i, i*fr:((i+1)*fr-1)]

    win_ji = m[i+1, i*fr:((i+1)*fr-1)] - m[i+1, (i+1)*fr:((i+1)*fr+fr-1)]

    plt.hist(win_ij, bins=50, density=True, label='Forward')
    plt.hist(-(win_ji), bins=50, density = True, label='Backward')
    plt.title('win_' + str(i+1) + '_to_' + str(i+2))
    plt.legend(frameon=False) 
    plt.show()


print('Done plotting overlap between neighbouring windows')
quit()
