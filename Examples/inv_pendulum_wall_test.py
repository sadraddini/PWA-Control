#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 12:11:15 2018

@author: sadra
"""

# Internal imports:
from inv_pendulum_wall import *
from main.visualization import visualize_set_tube,visualize_X_eps
from main.tree import intitialize_tree,Random_Tree_of_Polytopes
from main.tree_locator import tree_locator
from main.simulate import simulate_vanilla

intitialize_tree(s,T=55,alpha_start=0)
visualize_set_tube(s.X,-0.12,0.12,-1,1,tube_size=0.001)

Random_Tree_of_Polytopes(s,T_max=30)
    
visualize_X_eps(s,-0.12,0.12,-1,1,"angle","angular velocity","steps-to-go")

def animate(d):
    return visualize_subset_tree(s,d,-0.12,0.12,-1,1,"angle","angular velocity")
#visualize_subset_tree(s,10,-0.12,0.12,-1,1)

x0=np.array([-0.11,0.45]).reshape(2,1)
(x,u,G,theta,z,flag)=polytope_trajectory(s,x0,s.goal,10,-1)
