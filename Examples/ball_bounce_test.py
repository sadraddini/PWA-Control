#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 00:43:30 2018

@author: sadraddini
"""

# Internal imports:
from ball_bounce import *

from main.visualization import visualize_set_tube,visualize_subset_tree,visualize_X_eps_time,visualize_X_eps_cost
from main.tree import intitialize_tree,Random_Tree_of_Polytopes,tree_value_function
from main.tree_locator import tree_locator
from main.simulate import simulate_vanilla

intitialize_tree(s,T=15,alpha_start=0)
visualize_set_tube(s.X,-dmax,xmax,-vmax,vmax,tube_size=0.001)

Random_Tree_of_Polytopes(s,T_max=15)
    
visualize_X_eps_time(s,-dmax,xmax,-vmax,vmax,"Height","velocity","steps-to-go")
visualize_X_eps_time(s,-0.2,1.3,-5.5,5.5,"Height","velocity","steps-to-go")

s.Q=np.array([[1,-1],[-1,0]])
s.R=np.eye(1)*20
tree_value_function(s)
visualize_X_eps_cost(s,-0.2,1.3,-5.5,5.5,"Height","velocity","cost-to-go")

#
#from auxilary_methods import sample
#
#visualize_subset_tree(s,100,-dmax,xmax,-vmax,vmax)
#
def animate(d):
    return visualize_subset_tree(s,d,-0.2,1.3,-5.5,5.5,"Height","velocity")
#
#i=0
#x_sample=sample(s.l[i],s.u[i])

import pickle
s=pickle.load(open("ball_bounce.pkl","rb"))


### EVALUATION:
from main.trajectory import state_trajectory
from main.tree_locator import inside_tree,array_tree
from main.auxilary_methods import sample

N=500
T_max=80
sample_points=np.zeros((N,3))
x_sample={}
array_tree(s)
for iterations in range(N):
    print(iterations,"iterations")
    x_sample[iterations]=sample(s.l[0],s.u[0])
    sample_points[iterations,0]=int(inside_tree(s,x_sample[iterations]))
    for T in range(1,T_max):
        (x,u,z,flag)=state_trajectory(s,x_sample[iterations],s.goal,T)
        if flag==True:
            sample_points[iterations,1]=1
            sample_points[iterations,2]=T
            break

T_sample=np.zeros((N,1))
for iterations in range(N):
    print(iterations)
    T_sample[iterations]=tree_locator_time(s,x_sample[iterations]).time_to_go

appearence_time=np.ones((N,1))*100        
for iterations in range(1,N):
    print(iterations,"iterations")  
    STATES=zero_distance_states(s,x_sample[iterations])

x0=np.array([0.02,0]).reshape(2,1)
simulate_vanilla(s,x0)
visualize_set_tube_simulation(s,-dmax,xmax,-vmax,vmax,"Height","velocity","time_optimal")