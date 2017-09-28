#%% Loading python modules
import numpy as np
from numpy import genfromtxt
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
pi=np.pi


#%% Measuring the start time
time_start = time.time()


#%% Loading the dat file
simul_params1 = genfromtxt("s1.csv", delimiter=',')
simul_params2 = genfromtxt("s2.csv", delimiter=',')
E_over_V0 = genfromtxt("E_over_V0.csv", delimiter=',')
V0_star = genfromtxt("V0_star.csv", delimiter=',')
Tunneling_prob = genfromtxt("Tunneling_prob.csv", delimiter=',')
Ain = genfromtxt("Ain.csv", delimiter=',')
Aout = genfromtxt("Aout.csv", delimiter=',')

end = simul_params1[0]
step_size = simul_params1[1]
E_over_V0_max = simul_params1[2]
V0_star_max = simul_params1[3]

num_E = int(simul_params2[0])
num_V = int(simul_params2[1])



#%% With the post processing
#==============================================================================
# prob1 = 0.095
# prob2 = 0.105
# 
# E_over_V0_coordinate = np.where( (prob1 < Tunneling_prob) & (Tunneling_prob < prob2) == True )[0]
# V0_star_coordinate = np.where( (prob1 < Tunneling_prob) & (Tunneling_prob < prob2) == True )[1]
# 
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# 
# ax.scatter(E_over_V0[E_over_V0_coordinate, V0_star_coordinate], V0_star[E_over_V0_coordinate, V0_star_coordinate], Tunneling_prob[E_over_V0_coordinate, V0_star_coordinate], c='b', marker='o')
# #ax.scatter(Ain[E_over_V0_coordinate, V0_star_coordinate], Aout[E_over_V0_coordinate, V0_star_coordinate], Tunneling_prob[E_over_V0_coordinate, V0_star_coordinate], c='b', marker='o')
# 
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('Tunneling probability')
# #ax.set_zlim(0,1)
# 
# plt.show()
# 
# #plt.scatter(np.ravel(Ain[E_over_V0_coordinate, V0_star_coordinate]), np.ravel(Tunneling_prob[E_over_V0_coordinate, V0_star_coordinate]))
#==============================================================================


#%% Without the post processing (3D)
#==============================================================================
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# #ax.scatter(E_over_V0, V0_star, Tunneling_prob, c='r', marker='o')
# ax.plot_wireframe(E_over_V0, V0_star, Tunneling_prob)
# #ax.plot_surface(E_over_V0, V0_star, Tunneling_prob, cmap=cm.coolwarm)
# #ax.contour(E_over_V0, V0_star, Tunneling_prob, np.linspace(0, 1, 11))
#  
# ax.set_xlabel('E_over_V0')
# ax.set_ylabel('V0_star')
# ax.set_zlabel('Tunneling probability')
# #ax.set_zlim(0,1)
#  
# plt.show()
#  
# #plt.scatter(np.ravel(Ain), np.ravel(Tunneling_prob))
# #plt.scatter(np.ravel(Aout), np.ravel(Tunneling_prob))
#  
# #plt.scatter(np.ravel(np.exp(-2*Ain)), np.ravel(Tunneling_prob))
# #plt.scatter(np.ravel(np.exp(-2*Aout)), np.ravel(Tunneling_prob))
#  
# #plt.scatter(np.ravel(np.exp(Ain)), np.ravel(Tunneling_prob))
#==============================================================================


#%% Without the post processing (2D)
#==============================================================================
# plt.figure()
# CS = plt.contour(E_over_V0, V0_star, Tunneling_prob, np.linspace(0, 1, 21), colors='k')
# #CS = plt.contour(Ain, Aout, Tunneling_prob, np.linspace(0, 1, 21))
# #CS = plt.contour(np.exp(-2*Ain), np.exp(-2*Aout), Tunneling_prob, np.linspace(0, 1, 21))
# plt.clabel(CS, inline=1, fontsize=10)
# plt.title('Simplest default with labels')
#==============================================================================


#%% Ain, Aout 비교
#==============================================================================
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(Ain, Aout, Tunneling_prob, c='r', marker='o')
# ax.set_yscale('log')
# ax.set_xscale('log')
# 
# ax.set_xlabel('Ain')
# ax.set_ylabel('Aout')
# ax.set_zlabel('Tunneling probability')
# ax.set_zlim(0,1)
# 
# plt.show()
#==============================================================================
     

#%% Measuring the end time
time_end = time.time()


#%% Representing the simulation time                    
print('########################################################################')
print('Time Difference : %f' %(time_end-time_start))
