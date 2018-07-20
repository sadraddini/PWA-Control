#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 00:43:30 2018

@author: sadraddini
"""

# Internal imports:
from ball_bounce import *

from main.visualization import visualize_set_tube,visualize_X_eps,visualize_subset_tree
from main.tree import intitialize_tree,Random_Tree_of_Polytopes
from main.tree_locator import tree_locator
from main.simulate import simulate_vanilla

intitialize_tree(s,T=15,alpha_start=0)
visualize_set_tube(s.X,-dmax,xmax,-vmax,vmax,tube_size=0.001)

Random_Tree_of_Polytopes(s,T_max=15)
    
visualize_X_eps(s,-dmax,xmax,-vmax,vmax,"Height","velocity","steps-to-go")
visualize_X_eps(s,-0.2,1.3,-5.5,5.5,"Height","velocity","steps-to-go")

#
#from auxilary_methods import sample
#
#visualize_subset_tree(s,100,-dmax,xmax,-vmax,vmax)
#
def animate(d):
    return visualize_subset_tree(s,d,-0.2,1.3,-5.5,5.5,"Height","velocity")
#
#i=0
#x_sample=sample(s.l[i],s.u[i])

import pickle
s=pickle.load(open("ball_bounce.pkl","rb"))