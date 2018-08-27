#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 14:11:58 2018

@author: sadra
"""

from scipy.spatial import ConvexHull
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.animation as animation

### Internal imports
import sys
sys.path.append('..')

from main.auxilary_methods import vertices_cube
from main.ana_system import states_time_order,states_cost_order

def vertices_hull(s,list_of_states):
    PV=np.empty((0,s.n))
    for p in list_of_states:
        PV=np.vstack((PV,p.vertices+p.x.T))
    hull=ConvexHull(PV)
    return PV[hull.vertices,:]

def visualize_X_time_hull(s,list_of_states,xmin=-1,xmax=1,ymin=-1,ymax=1,xlabel='x_1',ylabel='x_2',title="interesting plot"):
    ax1 = plt.subplot(111)
    plt.figure(figsize=(20,20),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    # Trajectories
    p_list=[]
    for state in list_of_states:
        p_list.append(patches.Polygon(vertices_hull(s,[state]+[state.successor[0]]), True))
    p=PatchCollection(p_list,color=[0.4,0.8,0.4])
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_title(title)
    plt.show()
    return plt