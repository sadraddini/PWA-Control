#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 15:57:10 2018

@author: sadra
"""


### Internal imports
import sys
sys.path.append('../..')


# Internal imports:
from pushing import *
from main.tree import intitialize_tree,Random_Tree_of_Polytopes
from main.gurobi_m_library import trajectory_model


s.library={}
Tmax=3
for T in range(1,Tmax+1):
    print(T)
    trajectory_model(s,T)

intitialize_tree(s,T=Tmax,alpha_start=0)

Random_Tree_of_Polytopes(s,T_max=Tmax,eps_max=10)