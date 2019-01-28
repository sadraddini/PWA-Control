#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 08:36:24 2019

@author: sadra
"""

# External imports
import numpy as np

# My modules
from pypolycontain.lib.zonotope import zonotope
from pypolycontain.lib.polytope import polytope


# Internal imports
from PWA_lib.trajectory.system import system,linear_cell
from PWA_lib.trajectory.poly_trajectory import point_trajectory,polytopic_trajectory_given_modes

sys=system()
sys.name="inverted pendulum with two walls"


sys.A[0,0]=np.array([[1,0.01],[0.1,1]])
sys.B[0,0]=np.array([[0,0.01]]).T
sys.c[0,0]=np.array([[0,0]]).T

sys.A[1,0]=np.zeros((2,2))
sys.B[1,0]=np.zeros((2,1))
sys.c[1,0]=np.zeros((2,1))
sys.A[1,1]=np.array([[0,0],[-10,0]])
sys.B[1,1]=np.array([[0,0]]).T
sys.c[1,1]=np.array([[0,1]]).T


sys.A[2,0]=np.zeros((2,2))
sys.B[2,0]=np.zeros((2,1))
sys.c[2,0]=np.zeros((2,1))
sys.A[2,1]=np.array([[0,0],[-10,0]])
sys.B[2,1]=np.array([[0,0]]).T
sys.c[2,1]=np.array([[0,-1]]).T

H=np.array([[1,0,0],[0,1,0],[-1,0,0],[0,-1,0],[0,0,1],[0,0,-1]])
h=np.array([[0.12,1,0.12,1,4,4]]).T   
sys.C[0,0]=polytope(H,h)

H=np.array([[1,0,0],[0,1,0],[-1,0,0],[0,-1,0],[0,0,1],[0,0,-1]])
h=np.array([[0.1,1,0.12,1,4,4]]).T 
sys.C[1,0]=polytope(H,h)

H=np.array([[1,0,0],[0,1,0],[-1,0,0],[0,-1,0],[0,0,1],[0,0,-1]])
h=np.array([[0.12,1,0.1,1,4,4]]).T 
sys.C[2,0]=polytope(H,h)

H=np.array([[1,0,0],[0,1,0],[-1,0,0],[0,-1,0],[0,0,1],[0,0,-1]])
h=np.array([[0.12,1,-0.1,1,4,4]]).T  
sys.C[1,1]=polytope(H,h)

H=np.array([[1,0,0],[0,1,0],[-1,0,0],[0,-1,0],[0,0,1],[0,0,-1]])
h=np.array([[-0.1,1,0.12,1,4,4]]).T  
sys.C[2,1]=polytope(H,h)

sys.goal=zonotope(np.array([0,0.0]).reshape(2,1),np.array([[0,0],[0,0]]))

sys.n=2
sys.m=1
sys.list_of_sum_indices=[0,1,2]
sys.list_of_modes={}
sys.list_of_modes[0]=[0]
sys.list_of_modes[1]=[0,1]
sys.list_of_modes[2]=[0,1]

sys.build()
sys.scale=np.array([0.12,1])    

sys.build_cells()


from PWA_lib.polytree.tree import tree



mytree=tree(sys)
mytree.fill="continous"
mytree.system.name="Inverted Pendulum with Two Walls"

#list_of_samples=[np.random.random((2,1))*np.array([-0.12,1]).reshape(2,1)*(-1)**i for i in range(2)]
#mytree.construct_tree(list_of_samples,T_start=50,T_max=51,eps_point=0.0,eps_poly=0.1)
#mytree.visualize(axis_limit=[-0.12,0.12,-1,1])   
#
#assert 1==0
x_sample=np.array([0.06,0.75]).reshape(2,1)
T=50
mytree.extend_RRT(x_sample,T,eps_point=0.0,eps_poly=0.1)
mytree.visualize(axis_limit=[-0.12,0.12,-1,1])   

x_sample=np.array([0.0,0.9]).reshape(2,1)
T=10
mytree.extend_RRT(x_sample,T,eps_point=0.0,eps_poly=0.1)

mytree.visualize(axis_limit=[-0.12,0.12,-1,1])   


x_sample=np.array([0.01,0.7]).reshape(2,1)
T=30
mytree.extend_RRT(x_sample,T,eps_point=0.0,eps_poly=0.1)
mytree.visualize(axis_limit=[-0.12,0.12,-1,1])   

list_of_samples=[np.random.random((2,1))*np.array([-0.12,1]).reshape(2,1)*(-1)**i for i in range(27)]
mytree.construct_tree(list_of_samples,T_start=50,T_max=51,eps_point=0.0,eps_poly=0.1)

for x in mytree.states:
    x.p.color=(1,0,1)
mytree.visualize(axis_limit=[-0.12,0.12,-1,1]) 
assert 1==0

x_sample=np.array([-0.05,0.0]).reshape(2,1)
T=80
mytree.extend_RRT(x_sample,T,eps=0.1)

x_sample=np.array([0.0,-0.9]).reshape(2,1)
T=40
mytree.extend_RRT(x_sample,T,eps=0.01)

x_sample=np.array([0.0,0]).reshape(2,1)
T=60
mytree.extend_RRT(x_sample,T,eps=0.01)


mytree.visualize(axis_limit=[-0.12,0.12,-1,1])   