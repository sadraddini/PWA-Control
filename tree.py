#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 00:04:20 2018

@author: sadraddini
"""
import numpy as np
import random as random

from ana_system import tree
from trajectory import polytope_trajectory,make_state_trajectory_state_end
from auxilary_methods import sample

from tree_locator import array


def intitialize_tree(s,T=20,alpha_start=0):
    x0=np.array([0.5,0]).reshape(2,1)
    goal=s.goal
    (x,u,G,theta,z,flag)=polytope_trajectory(s,x0,goal,T,alpha_start)
    if flag==True:
        make_state_trajectory_state_end(s,x,u,z,G,theta,T,goal)
        s.tree=tree(s.goal)
        s.X.append(s.goal)
        tree.nodes=s.X
        array(s,s.X)

def extend_RRT(s,tree):
    i=random.choice(s.modes)
    x_sample=sample(s.l[i],s.u[i])
    array(s,s.X)
    if inside_tree(s,x_sample)==True:
        print("inside X_tree")
        return False
    else:
        
    
        

def tree_construct(s,goal,max_iterations):
    pass