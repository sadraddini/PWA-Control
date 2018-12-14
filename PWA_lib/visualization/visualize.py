# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 17:47:23 2018

@author: sadra
"""

import numpy as np

from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon


from pyinpolytope.utilities.utils import vertices_cube as vcube


def add_tube(ax,x,G,eps=0,list_of_dimensions=[0,1],axis=0):
    """
    Only when G is a zonotope generator
    """
    T=len(x.keys())
    p_list=[]
    n=x[0].shape[0]
    points_convex_all=np.empty((0,2))
    for t in range(0,T-1):
        y0=x[t].T+np.dot(G[t]+eps*np.eye(n),vcube(G[t].shape[1]).T).T
        y1=x[t+1].T+np.dot(G[t+1]+eps*np.eye(n),vcube(G[t+1].shape[1]).T).T
        y=np.vstack((y0,y1))
#        print y.shape
        points=y[:,list_of_dimensions]
        points_convex=points[ConvexHull(points).vertices,:]
#        print points_convex
        p=Polygon(points_convex)
        p_list.append(p)
        points_convex_all=np.vstack((points_convex_all,points_convex))
    p_patch = PatchCollection(p_list, color=(0.5,0.5,0.5),alpha=0.75)
    ax.add_collection(p_patch)
    ax.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    if axis==1:
        ax.set_xlim([min([x[t][0,0] for t in range(T)]),max([x[t][0,0] for t in range(T)])])
        ax.set_ylim([min([x[t][1,0] for t in range(T)]),max([x[t][1,0] for t in range(T)])])
    elif axis==0:
        ax.set_xlim([min(points_convex_all[:,0]),max(points_convex_all[:,0])])
        ax.set_ylim([min(points_convex_all[:,1]),max(points_convex_all[:,1])])
    return ax