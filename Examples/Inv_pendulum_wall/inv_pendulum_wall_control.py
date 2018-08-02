#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 16:47:31 2018

@author: sadra
"""

from inv_pendulum_wall import *

import pickle
s=pickle.load(open("inv_pendulum_wall.pkl","rb"))


### Internal imports
import sys
sys.path.append('../..')
from main.trajectory import state_trajectory
from main.tree_locator import inside_tree,array_tree
from main.auxilary_methods import sample
from main.visualization import visualize_subset_tree,visualize_set_tube_simulation
from main.simulate import simulate_vanilla

array_tree(s)

x0=np.array([0.05,0.75]).reshape(2,1)
simulate_vanilla(s,x0)
visualize_set_tube_simulation(s,-0.12,0.12,-1,1,"Angle","Angular Velocity","time_optimal")