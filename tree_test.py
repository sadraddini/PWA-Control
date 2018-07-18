#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 00:43:30 2018

@author: sadraddini
"""
from random import randint

# Internal imports:
from ball_bounce import *
from visualization import visualize_set_tube
from tree import intitialize_tree,Random_Tree_of_Polytopes
from tree_locator import tree_locator
from simulate import simulate_vanilla

intitialize_tree(s,T=15,alpha_start=0)
visualize_set_tube(s.X,-dmax,xmax,-vmax,vmax,tube_size=0.001)

Random_Tree_of_Polytopes(s,T_max=6)
    
visualize_X_eps(s,-dmax,xmax,-vmax,vmax)
#x=np.array([0.5,3]).reshape(2,1)
#x0=tree_locator(s,x)
#u=simulate_vanilla(s,x)