# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:41:58 2018

@author: sadra
"""

# External imports
import numpy as np

# My modules
from pypolycontain.lib.zonotope import zonotope
from pypolycontain.lib.polytope import polytope


# Internal imports
from PWA_lib.trajectory.system import system
from PWA_lib.trajectory.poly_trajectory import point_trajectory

sys=system()
sys.name="inverted pendulum with wall"


sys.A[1,0]=np.array([[1,0.01],[0.1,1]])
sys.A[1,1]=np.array([[1,0.01],[-9.9,1]])
sys.B[1,0]=np.array([[0,0.01]]).T
sys.B[1,1]=np.array([[0,0.01]]).T
sys.c[1,0]=np.array([[0,0]]).T
sys.c[1,1]=np.array([[0,1]]).T

sys.A[1,1]=sys.A[1,0]
sys.B[1,1]=sys.B[1,0]
sys.c[1,1]=sys.c[1,0]

H=np.array([[1,0,0],[0,1,0],[-1,0,0],[0,-1,0],[0,0,1],[0,0,-1]])
h=np.array([[0.10,1,0.12,1,4,4]]).T   
sys.C[1,0]=polytope(H,h)


H=np.array([[1,0,0],[0,1,0],[-1,0,0],[0,-1,0],[0,0,1],[0,0,-1]])
h=np.array([[0.12,1,-0.1,1,4,4]]).T  
sys.C[1,1]=polytope(H,h)

sys.goal=zonotope(np.array([0,0]).reshape(2,1),np.array([[0.0,0],[0,0]]))


sys.n=2
sys.m=1
sys.list_of_sum_indices=[1]
sys.list_of_modes={}
sys.list_of_modes[1]=[0,1]

#sys.build()


import matplotlib.pyplot as plt

x0=np.array([0.11,-0.4]).reshape(2,1)
T=50
(x,u,delta_PWA,mu)=point_trajectory(sys,x0,[sys.goal],T)

plt.plot([x[t][0,0] for t in range(T+1)],[x[t][1,0] for t in range(T+1)])
plt.plot([x[t][0,0] for t in range(T+1)],[x[t][1,0] for t in range(T+1)],'+')
plt.plot([0],[0],'o')
print np.array([delta_PWA[t,1,1] for t in range(T)])