#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 00:43:30 2018

@author: sadraddini
"""


# Internal imports:
from ball_bounce import *

### Internal imports
import sys
sys.path.append('../..')


from main.visualization import visualize_set_tube,visualize_subset_tree,visualize_X_eps_time,visualize_X_eps_cost,visualize_X_time_hull_eps
from main.tree import intitialize_tree,Random_Tree_of_Polytopes,tree_value_function
from main.gurobi_m_library import trajectory_model


Tmax=20
s.library={}
s.failure_tree_extentions=[]
for T in range(1,Tmax+1):
    print(T)
    trajectory_model(s,T)
    
    
intitialize_tree(s,T=Tmax,alpha_start=0)
visualize_set_tube(s.X,-dmax,xmax,-vmax,vmax,tube_size=0.001)

Random_Tree_of_Polytopes(s,T_max=Tmax,eps_max=0.1)
    
visualize_X_eps_time(s,-dmax,xmax,-vmax,vmax,"Height","velocity","")
visualize_X_time_hull_eps(s,-dmax,xmax,-vmax,vmax,r"Height",r"velocity","")

def animate(d):
    return visualize_subset_tree(s,d,-0.2,1.3,-5.5,5.5,"Height","velocity")