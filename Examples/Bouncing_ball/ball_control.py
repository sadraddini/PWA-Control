#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 16:38:18 2018

@author: sadra
"""
from ball_bounce import *

import pickle
s=pickle.load(open("ball_bounce.pkl","rb"))


### Internal imports
import sys
sys.path.append('..')
from main.trajectory import state_trajectory
from main.tree_locator import inside_tree,array_tree
from main.auxilary_methods import sample
from main.visualization import visualize_subset_tree,visualize_set_tube_simulation
from main.simulate import simulate_vanilla

array_tree(s)

def animate(d):
    return visualize_subset_tree(s,d,-0.2,1.3,-5.5,5.5,"Height","velocity")

x0=np.array([0.02,0]).reshape(2,1)
simulate_vanilla(s,x0)
visualize_set_tube_simulation(s,-dmax,xmax,-vmax,vmax,"Height","velocity","time_optimal")