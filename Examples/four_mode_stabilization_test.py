#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 13:50:45 2018

@author: sadra
"""

from four_mode_stabilization import *

from visualization import visualize_set_tube,visualize_X_eps
from tree import intitialize_tree,Random_Tree_of_Polytopes
from tree_locator import tree_locator
from simulate import simulate_vanilla

intitialize_tree(s,T=15,alpha_start=0)
visualize_set_tube(s.X,-5,5,-5,5,tube_size=0.001)

Random_Tree_of_Polytopes(s,T_max=10)
    
visualize_X_eps(s,-5,5,-5,5)