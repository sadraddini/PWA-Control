#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 14:37:00 2018

@author: sadra
"""

from ball_bounce import *


from trajectory import polytope_trajectory,make_state_trajectory_state_end
from visualization import visualize_set_tube

x0=np.array([0.5,0]).reshape(2,1)
T=20
alpha_start=-1
(x,u,G,theta,z,flag)=polytope_trajectory(s,x0,goal,T,alpha_start)
goal=s.goal
if flag==True:
    make_state_trajectory_state_end(s,x,u,z,G,theta,T,goal)
visualize_set_tube(s.X,-dmax,xmax,-vmax,vmax,tube_size=0.001)