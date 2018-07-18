#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 00:04:20 2018

@author: sadraddini
"""
import numpy as np
from random import random,choice,randint

from main.trajectory import polytope_trajectory,make_state_trajectory_state_end
from main.auxilary_methods import sample
from main.tree_locator import array_tree,inside_tree,sorted_distance_states,all_vertices_out_of_tree
from main.simulate import simulate_0


def intitialize_tree(s,T=20,x0=np.array([0,0]).reshape(2,1),alpha_start=0):
    goal=s.goal
    (x,u,G,theta,z,flag)=polytope_trajectory(s,x0,goal,T,alpha_start)
    if flag==True:
        make_state_trajectory_state_end(s,x,u,z,G,theta,T,goal)
        s.X.append(s.goal)
        array_tree(s)
        s.tree_iterations+=1
        s.tree_size[s.tree_iterations]=len(s.X)

def extend_RRT(s,T,alpha_start=10**5,eps=0.1):
    i=choice(s.modes)
    x_sample=sample(s.l[i],s.u[i])
    array_tree(s)
    if inside_tree(s,x_sample)==True:
        print("inside X_tree")
        return False
    else:
        (x_zero,T)=simulate_0(s,x_sample,T)
        STATES=sorted_distance_states(s,x_zero)
        for state_end in STATES:
            (x,u,G,theta,z,flag)=polytope_trajectory(s,x_sample,state_end,T,alpha_start,eps)
            if flag==True:
                if all_vertices_out_of_tree(s,x[0],G[0])==True:
                    make_state_trajectory_state_end(s,x,u,z,G,theta,T,state_end)
                    return True
                
def Random_Tree_of_Polytopes(s,T_max=10):
    for t in range(0,500):
        print("*"*100,"iteration",t)
        T=randint(1,T_max)
        flag=extend_RRT(s,T,alpha_start=-1,eps=random())
        if flag==True:
            s.tree_iterations+=1
            s.tree_size[s.tree_iterations]=len(s.X)
                       
def rewiring(s,x):
    pass

def tree_construct(s,goal,max_iterations):
    pass