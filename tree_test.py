#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 00:43:30 2018

@author: sadraddini
"""
from ball_bounce import *

from visualization import visualize_set_tube

from tree import intitialize_tree

from tree_locator import tree_locator

from simulate import simulate_vanilla

intitialize_tree(s,T=20,alpha_start=0)
visualize_set_tube(s.X,-dmax,xmax,-vmax,vmax,tube_size=0.001)


x=np.array([0.5,3]).reshape(2,1)
x0=tree_locator(s,x)
u=simulate_vanilla(s,x)