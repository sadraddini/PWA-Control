#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 22 18:03:08 2018

@author: sadra
"""

import numpy as np

import matplotlib
#matplotlib.use('Agg')

from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

from main.auxilary_methods import vertices_cube

def visualize_funnel(s,iterations,xmin=-1,xmax=1,ymin=-1,ymax=1,xlabel='x_1',ylabel='x_2',title="interesting plot"):
    fig,ax1 = plt.subplots()
    #fig=plt.figure(figsize=(20,20),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    p_list=[]
    list_of_states=s.X[0:s.tree_size[iterations]]
    for state in list_of_states:
        if state.mode==state.successor[0].mode:
            p_list.append(patches.Polygon(vertices_hull_eps(s,[state]+[state.successor[0]]), True))
        else:
            p_list.append(patches.Polygon(vertices_hull_eps(s,[state]), True))
    p=PatchCollection(p_list,color=(0.4,0.4,0.4))
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_title(title)
    ax1.plot([0,0],[-5,5],'black',linewidth=0.5)
    fig.savefig('figures/bouncing_ball_funnel_%d.png'%iterations, dpi=100)
    return fig

def vertices_hull_eps(s,list_of_states):
    """
    Convex hull computation of list of states
    """
    PV=np.empty((0,s.n))
    for p in list_of_states:
        PV=np.vstack((PV,p.vertices_eps+p.x.T))
    hull=ConvexHull(PV)
    return PV[hull.vertices,:]

def show_ball(s,t):
    x=0
    y=s.traj[t][0]
    p=[]
    p1=[]
    p2=[]
    fig,ax1 = plt.subplots()
    r=0.1
    ax1.set_aspect('equal')
    # Ground
    p1.append(patches.Polygon(np.array([[-0.5,-0.5,0.5,0.5],[0,-0.5,-0.5,0]]).T, True))
    p1=PatchCollection(p1,color=(0.5,0.7,1))
    # Goal
    p2.append(patches.Polygon(np.array([[-0.2,-0.2,0.2,0.2],[1.2,1,1,1.2]]).T, True))
    p2=PatchCollection(p2,color=(0.5,1,0.5))    
    # Ball
    circle = patches.Circle([x,y], r)
    #fig=plt.figure(figsize=(10,10),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlim([-0.5,0.5])
    ax1.set_ylim([-0.5,1.2])
    p.append(circle)
    p=PatchCollection(p,color=(0.99,0.5,0.5))
    ax1.add_collection(p2)
    ax1.add_collection(p1)
    ax1.add_collection(p)
    ax1.set_facecolor((0, 0.0, 0.0))
    ax1.set_title("",fontsize=20)
    fig.savefig('figures/ball_%d.png'%t, dpi=100)
    return fig

if False:
    for iteration in range(1,s.tree_iterations):
        fig=visualize_funnel(s,iteration,xmin=-0.1,xmax=1.2,ymin=-5,ymax=5,xlabel=r'Height',ylabel=r'Velocity',title="Bouncing Ball Iteration - %d"%iteration)
    
if True:
    for t in range(len(s.traj)):
        fig=show_ball(s,t)
    