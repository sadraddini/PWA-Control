#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 12:11:15 2018

@author: sadra
"""

# Internal imports:
from inv_pendulum_wall import *
from visualization import visualize_set_tube,visualize_X_eps
from tree import intitialize_tree,Random_Tree_of_Polytopes
from tree_locator import tree_locator
from simulate import simulate_vanilla

intitialize_tree(s,T=100,alpha_start=0)
visualize_set_tube(s.X,-0.12,0.12,-1,1,tube_size=0.001)

Random_Tree_of_Polytopes(s,T_max=15)
    
visualize_X_eps(s,-0.12,0.12,-1,1)

#visualize_subset_tree(s,10,-0.12,0.12,-1,1)