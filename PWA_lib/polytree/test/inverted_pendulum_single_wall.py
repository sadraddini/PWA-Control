#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 08:45:48 2019

@author: sadra
"""

# External imports
import numpy as np
import pickle


# My modules
from pypolycontain.lib.zonotope import zonotope
from pypolycontain.lib.polytope import polytope
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes_ax as visZ,plt


# Internal imports
from PWA_lib.trajectory.system import system,linear_cell
from PWA_lib.trajectory.poly_trajectory import point_trajectory,polytopic_trajectory_given_modes



sys=system()
sys.name="inverted pendulum with wall"


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
sys.list_of_sum_indices=[0,1]
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
mytree.system.name="Inverted Pendulum with Wall"
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 20,
        }

N=50
list_of_samples=[np.random.random((2,1))*np.array([0.24,2]).reshape(2,1)+np.array([-0.12,-1]).reshape(2,1) for i in range(N)]
mytree.construct_tree(list_of_samples,T_start=20,T_max=51,eps_point=0.0,eps_poly=0.1)

def visualize_tree(tree):
    fig,ax=plt.subplots()
    axis_limit=[-0.15,0.15,-1.2,1.2]
    visZ(ax,[x.p for x in tree.states],axis_limit=axis_limit,title="%s - %d zonotopes %d branches" %(tree.name,len(tree.states),len(tree.branches)))
    ax.set_xlabel(r"$\theta$",fontdict=font)
    h=ax.set_ylabel(r"$\dot{\theta}$",fontdict=font)
    h.set_rotation(0)
    ax.plot([0.1,0.1],[-1,1],color=(0,0,0),linewidth=2,linestyle='dashed')
    ax.plot([0.12,0.12,-0.12,-0.12,0.12],[-1,1,1,-1,-1],color=(0,0,0), linewidth=2)

visualize_tree(mytree)    

""" Save File"""
pickle.dump(mytree,open("inverted_pendulum.pkl","w"))

"""Coverage Evaluation"""
if False:
    N=1000
    list_of_samples=[np.random.random((2,1))*np.array([0.24,2]).reshape(2,1)+np.array([-0.12,-1]).reshape(2,1) for i in range(N)]
    inside_MILP=0
    inside_tree=0
    inside_tree_and_MILP=0
    i=0
    T=80
    for x_sample in list_of_samples:
        i+=1
        print "point :", i
        (x_nom,u_nom,delta_PWA,mu_nom,flag_MILP)=point_trajectory(sys,x_sample,[sys.goal],T,eps=0)
        flag_tree=mytree.inside_tree(x_sample)
        inside_MILP+=flag_MILP
        inside_tree+=flag_tree
        inside_tree_and_MILP+=flag_MILP*flag_tree
    print N,inside_MILP,inside_tree,inside_tree_and_MILP
