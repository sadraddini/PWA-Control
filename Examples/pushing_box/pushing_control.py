#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 18:31:43 2018

@author: sadra
"""

### Internal imports
import sys
sys.path.append('../..')
from main.trajectory import state_trajectory
from main.tree_locator import inside_tree,array_tree
from main.auxilary_methods import sample
from main.visualization import visualize_subset_tree,visualize_set_tube_simulation
from main.simulate import simulate_vanilla,simulate_convex

array_tree(s)

x0=np.array([4,0.5,0,0,1,0]).reshape(6,1)
simulate_convex(s,x0)
visualize_set_tube_simulation(s,-0.12,0.12,-1,1,"Angle","Angular Velocity","time_optimal")