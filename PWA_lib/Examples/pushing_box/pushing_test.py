#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 15:57:10 2018

@author: sadra
"""

sfda
### Internal imports
import sys
sys.path.append('../..')


# Internal imports:
from pushing import *
from main.tree import intitialize_tree,Random_Tree_of_Polytopes
from main.gurobi_m_library import trajectory_model
from main.visualization_high import visualize_proj_eps,visualize_proj_eps_states_simulation,visualize_proj_eps_states

                   

s.failure_tree_extentions=[]
s.library={}
Tmax=10
for T in range(1,Tmax+1):
    print(T)
    trajectory_model(s,T)

intitialize_tree(s,T=Tmax,alpha_start=0)

#x_sample_list=[np.array([7,-0.3,0.2,0,1,0]).reshape(6,1)]
Random_Tree_of_Polytopes(s,T_max=Tmax,eps_max=5)
#Random_Tree_of_Polytopes(s,T_max=Tmax,eps_max=5,x_sample_list=x_sample_list)

visualize_proj_eps(s,0,1,xmin=0,xmax=10,ymin=-1.5,ymax=1.5,xlabel=r'$x$',ylabel=r'$y$',title=r"planar pushing")
visualize_proj_eps(s,1,2,xmin=-1.5,xmax=1.5,ymin=-0.4,ymax=0.4,xlabel=r'$y$',ylabel=r'$\theta$',title="planar pushing")
