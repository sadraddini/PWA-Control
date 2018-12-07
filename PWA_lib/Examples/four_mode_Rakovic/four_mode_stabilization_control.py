#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 16:59:00 2018

@author: sadra
"""

from four_mode_stabilization import *

import pickle
s=pickle.load(open("Rakovic4.pkl","rb"))


### Internal imports
import sys
sys.path.append('../..')
from main.trajectory import state_trajectory
from main.tree_locator import inside_tree,array_tree
from main.auxilary_methods import sample
from main.visualization import visualize_subset_tree,visualize_set_tube_simulation
from main.simulate import simulate_vanilla

array_tree(s)

x0=np.array([3,2]).reshape(2,1)
simulate_vanilla(s,x0)
visualize_set_tube_simulation(s,-5,5,-5,5,"Angle","Angular Velocity","time_optimal")