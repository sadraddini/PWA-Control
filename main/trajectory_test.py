#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 14:37:00 2018

@author: sadra
"""

from ball_bounce import *


from main.trajectory import polytope_trajectory,make_state_trajectory_state_end
from main.visualization import visualize_set_tube

from tree_locator import array_tree,d_tree
from tree import extend_RRT

x0=np.array([0.2,0]).reshape(2,1)
T=25
alpha_start=0
(x,u,G,theta,z,flag)=polytope_trajectory(s,x0,s.goal,T,alpha_start)
if flag==True:
    make_state_trajectory_state_end(s,x,u,z,G,theta,T,s.goal)
visualize_set_tube(s.X,-dmax,xmax,-vmax,vmax,tube_size=0.001)

for t in range(0,100):
    extend_RRT(s,5,alpha_start=10**5)
    