#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 00:04:20 2018

@author: sadraddini
"""
import numpy as np
from random import random,choice,randint

from ana_system import tree
from trajectory import polytope_trajectory,make_state_trajectory_state_end
from auxilary_methods import sample

from tree_locator import array_tree,inside_tree,sorted_distance_states,all_vertices_out_of_tree
from simulate import simulate_0


def intitialize_tree(s,T=20,alpha_start=0):
    x0=np.array([0.5,0]).reshape(2,1)
    goal=s.goal
    (x,u,G,theta,z,flag)=polytope_trajectory(s,x0,goal,T,alpha_start)
    if flag==True:
        make_state_trajectory_state_end(s,x,u,z,G,theta,T,goal)
        s.tree=tree(s.goal)
        s.X.append(s.goal)
        tree.nodes=s.X
        array_tree(s)

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
        extend_RRT(s,T,alpha_start=-1,eps=random())
                       
def rewiring(s,x):
    pass

def tree_construct(s,goal,max_iterations):
    pass