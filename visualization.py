#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 10:59:22 2018

@author: sadra
"""

from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from auxilary_methods import vertices_cube
import numpy as np
import random as random

def visualize_set_tube(list_of_states,xmin=-1,xmax=1,ymin=-1,ymax=1,tube_size=0.01):
    from matplotlib.collections import PatchCollection
    ax1 = plt.subplot(111)
    plt.figure(figsize=(30,10))
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    vertices_0=np.array([[1,1,-1,-1],[1,-1,-1,1]])
    # Trajectories
    p_list=[]
    for state in list_of_states:
        vertices=vertices_cube(state.G.shape[1]).T
        v=np.dot(state.G,vertices)
        v=Minkowski_hull(v.T,(vertices_0.T*tube_size)).T
        p_list.append(patches.Polygon(v.T+state.x.T, True))
    p=PatchCollection(p_list, alpha=0.4,color=(0.5,0.2,0.7))
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    plt.show()
    
def Minkowski_hull(p1,p2):
    """
    Inputs: 
        p1: n1points * ndim matrix
        p2: n2points * ndim matrix
    Decription:
        Compute the convex hull of minkowski sum of p1 and p2
    Output:
        p: points 
    """
    (n1,d)=p1.shape
    (n2,d)=p2.shape
    out=np.ones((n1*n2,d))
    for i in range(n2):
        out[n1*i:n1*(i+1),:]=p1+p2[i,:]
    return out[ConvexHull(out).vertices,:]