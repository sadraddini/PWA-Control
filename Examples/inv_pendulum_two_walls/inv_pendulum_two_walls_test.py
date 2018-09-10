#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 14:32:13 2018

@author: sadra
"""

from inv_pendulum_two_walls import *

### Internal imports
import sys
sys.path.append('../..')


# Internal imports:
from main.visualization import visualize_set_tube,visualize_X_eps_time,visualize_X_eps_cost,visualize_subset_tree,visualize_X_time_hull_eps
from main.tree import intitialize_tree,Random_Tree_of_Polytopes
from main.tree_locator import tree_locator_time
from main.simulate import simulate_vanilla
from main.gurobi_m_library import trajectory_model


s.library={}
Tmax=15
for T in range(1,Tmax+1):
    print(T)
    trajectory_model(s,T)

intitialize_tree(s,T=Tmax ,alpha_start=0)
visualize_set_tube(s.X,-0.12,0.12,-1,1,tube_size=0.001)

Random_Tree_of_Polytopes(s,T_max=Tmax,eps_max=0)
    
visualize_X_eps_time(s,-0.12,0.12,-1,1,"angle","angular velocity","steps-to-go")
visualize_X_time_hull_eps(s,-0.12,0.12,-1,1,"angle","angular velocity","steps-to-go")

def animate(d):
    return visualize_subset_tree(s,d,-0.12,0.12,-1,1,"angle","angular velocity")
#visualize_subset_tree(s,10,-0.12,0.12,-1,1)