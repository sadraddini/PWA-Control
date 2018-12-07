#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 18:31:43 2018

@author: sadra
"""

### Internal imports
import sys
sys.path.append('../..')
from main.trajectory import state_trajectory
from main.tree_locator import inside_tree,array_tree
from main.auxilary_methods import sample
from main.visualization import visualize_subset_tree,visualize_set_tube_simulation
from main.simulate import simulate_vanilla,simulate_convex,simulate_MPC
from main.visualization_high import visualize_proj_eps,visualize_proj_eps_states_simulation,visualize_proj_eps_states

#visualize_proj_eps(s,0,1,xmin=0,xmax=10,ymin=-1,ymax=1,xlabel='x',ylabel='y',title="planar pushing")
#visualize_proj_eps(s,1,2,xmin=-1,xmax=1,ymin=-0.5,ymax=0.5,xlabel='y',ylabel='theta',title="planar pushing")

array_tree(s)

x0=np.array([0,1.5,0.1,0,1,0]).reshape(6,1)
simulate_convex(s,x0)

visualize_proj_eps_states_simulation(s,0,1,xmin=0,xmax=10,ymin=-1.7,ymax=1.7,xlabel=r'$x$',ylabel=r'$y$',title="interesting plot")
visualize_proj_eps_states_simulation(s,1,2,xmin=-1.7,xmax=1.7,ymin=-0.5,ymax=0.5,xlabel=r'$y$',ylabel=r'$\theta$',title="interesting plot")
visualize_proj_eps_states_simulation(s,1,3,xmin=-1.7,xmax=1.7,ymin=-0.15,ymax=0.15,xlabel=r'$y$',ylabel=r'$\delta$',title="interesting plot")
visualize_proj_eps_states(s.X,0,1,xmin=0,xmax=10,ymin=-1,ymax=1,xlabel='x_1',ylabel='x_2',title="interesting plot")