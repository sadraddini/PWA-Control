#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 18:54:28 2018

@author: sadra
"""

import numpy as np

import matplotlib
matplotlib.use('Agg')

import pylab

from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

### Internal imports
import sys
sys.path.append('..')

def visualize_push_box(x,t,xmin=0,xmax=10,ymin=-3,ymax=3,xlabel='x',ylabel='y',title="planar pushing"):
    fig,ax1 = plt.subplots()
    #fig=plt.figure(figsize=(8,5),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    v=np.array([[-0.2, -0.2],[-0.2,  0.2],[0.2, 0.2],[0.2,  -0.2]])*3
    theta=x[2]
    R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]).reshape((2,2))
    v_box=np.dot(R,v.T)+x[0:2,:]
    # Trajectories
    p=patches.Polygon(v_box[:,:].T, True)
    p=PatchCollection([p],color=(0.9,0.4,0.3))
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_title(title)
    alpha=2
    c_point=x[0:2,:]+np.dot(R,np.array([-0.7,6*np.asscalar(x[3])]).reshape((2,1)))
    f_point=x[4:,:]*alpha
    ax1.arrow(np.asscalar(c_point[0]-f_point[0]),np.asscalar(c_point[1]-f_point[1]),np.asscalar(f_point[0]),np.asscalar(f_point[1]),head_width=0.2)
    fig.savefig('figures/planar_pushing_%d.png'%t, dpi=100)
    return fig

for t in range(len(s.traj)):
    fig=visualize_push_box(s.traj[t],t,xmin=0,xmax=10,ymin=-3,ymax=3,xlabel='x',ylabel='y',title="planar pushing")
    #fig.savefig('planar_pushing_bizim_%d.png'%t, dpi=100)
    