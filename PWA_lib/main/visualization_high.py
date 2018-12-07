#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 14:11:58 2018

@author: sadra
"""
import numpy as np


from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

### Internal imports
import sys
sys.path.append('..')

from main.auxilary_methods import vertices_cube
from main.ana_system import states_time_order,states_cost_order

def visualize_proj_eps(s,di,dj,xmin=-1,xmax=1,ymin=-1,ymax=1,xlabel='x_1',ylabel='x_2',title="interesting plot"):
    """
    di: x dimension
    dj: y dimension
    """
    ax1 = plt.subplot(111)
    plt.figure(figsize=(20,20),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel(xlabel,fontsize=20)
    ax1.set_ylabel(ylabel,fontsize=20)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    vertices_0=vertices_cube(s.n)
    # Trajectories
    p_list=[]
    STATES=states_time_order(s)
    for state_X in STATES[::-1]:
        v=np.dot(state_X.G_eps,vertices_0.T)+state_X.x
        v=v[[di,dj],:].T
        hull=ConvexHull(v)
        p_list.append(patches.Polygon((v[hull.vertices,:]), True))
    max_T=max([state.time_to_go for state in STATES])
    p=PatchCollection(p_list,color=[(state.time_to_go/max_T,1-state.time_to_go/max_T,0) for state in STATES[::-1]])
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_title(title,fontsize=20)
    return plt

def visualize_proj_eps_states(list_of_states,di,dj,xmin=-1,xmax=1,ymin=-1,ymax=1,xlabel='x_1',ylabel='x_2',title="interesting plot"):
    """
    di: x dimension
    dj: y dimension
    """
    ax1 = plt.subplot(111)
    plt.figure(figsize=(20,20),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel(xlabel,fontsize=20)
    ax1.set_ylabel(ylabel,fontsize=20)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    n=list_of_states[0].x.shape[0]
    vertices_0=vertices_cube(n)
    # Trajectories
    p_list=[]
    for state_X in list_of_states:
        v=np.dot(state_X.G_eps,vertices_0.T)+state_X.x
        v=v[[di,dj],:].T
        hull=ConvexHull(v)
        p_list.append(patches.Polygon((v[hull.vertices,:]), True))
    p=PatchCollection(p_list,color=(0.4,0.4,0.3))
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_title(title,fontsize=20)
    return plt

def visualize_proj_eps_states_simulation(s,di,dj,xmin=-1,xmax=1,ymin=-1,ymax=1,xlabel='x_1',ylabel='x_2',title="interesting plot"):
    """
    di: x dimension
    dj: y dimension
    """
    ax1 = plt.subplot(111)
    plt.figure(figsize=(20,20),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel(xlabel,fontsize=20)
    ax1.set_ylabel(ylabel,fontsize=20)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    n=s.n
    vertices_0=vertices_cube(n)
    # Trajectories
    p_list=[]
    for state_X in s.X:
        v=np.dot(state_X.G_eps,vertices_0.T)+state_X.x
        v=v[[di,dj],:].T
        hull=ConvexHull(v)
        p_list.append(patches.Polygon((v[hull.vertices,:]), True))
    p=PatchCollection(p_list,color=(0.3,0.3,0.4))
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_title(title,fontsize=20)
    trajectory=np.array(s.traj)
    ax1.plot(trajectory[:,di,:],trajectory[:,dj,:],'r',linewidth=0.4)
    ax1.plot(trajectory[:,di,:],trajectory[:,dj,:],'ro')
    return plt

def finder_successor(s,N):
    for n in range(len(s.X)):
        if s.X[N].successor[0]==s.X[n]:
            return n
 