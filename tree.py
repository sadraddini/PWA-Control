#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 00:04:20 2018

@author: sadraddini
"""
import numpy as np

from ana_system import tree
from trajectory import polytope_trajectory,make_state_trajectory_state_end


def intitialize_tree(s,T=20,alpha_start=0):
    x0=np.array([0.5,0]).reshape(2,1)
    goal=s.goal
    (x,u,G,theta,z,flag)=polytope_trajectory(s,x0,goal,T,alpha_start)
    if flag==True:
        make_state_trajectory_state_end(s,x,u,z,G,theta,T,goal)
        s.tree=tree(s.goal)
        s.X.append(s.goal)
        tree.nodes=s.X
        for 

def extend_RRT(s,tree):
    pass

def tree_construct(s,goal,max_iterations):
    pass