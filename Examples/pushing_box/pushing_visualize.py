#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 18:54:28 2018

@author: sadra
"""

import numpy as np

import matplotlib
matplotlib.use('Agg')

from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

### Internal imports
import sys
sys.path.append('../..')

from main.auxilary_methods import vertices_cube,find_mode_control

def visualize_push_box(s,t,xmin=0,xmax=10,ymin=-3,ymax=3,xlabel='x',ylabel='y',title="planar pushing"):
    fig,ax1 = plt.subplots()
    #fig=plt.figure(figsize=(8,5),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    v=np.array([[-0.2, -0.2],[-0.2,  0.2],[0.2, 0.2],[0.2,  -0.2]])*3
    p_list=[]
    for tau in range(t):
        x=s.traj[tau]
        theta=x[2]
        R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]).reshape((2,2))
        v_box=np.dot(R,v.T)+x[0:2,:]
    # Trajectories
        p=patches.Polygon(v_box[:,:].T, True)
        p_list.append(p)
    # End
    p_list=PatchCollection(p_list,color=(0.9,0.4,0.3),alpha=0.01)
    ax1.add_collection(p_list)
    ax1.add_collection(PatchCollection([p],color=(0.9,0,0),alpha=0.5))
    theta=0
    R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]).reshape((2,2))
    v_box=np.dot(R,v.T)+np.array([10,0]).reshape(2,1)
    p=patches.Polygon(v_box[:,:].T, True)
    ax1.add_collection(PatchCollection([p],color=(0,0.9,0),alpha=0.5))
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_title(title)
    props=dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    i=find_mode_control(s,x)
    if i==2:
        textstr="sticking"
    elif i==3:
        textstr="sliding up"
    elif i==4:
        textstr="sliding down"
    ax1.text(0.05,0.15, textstr, transform=ax1.transAxes, fontsize=14,verticalalignment='top', bbox=props)
    alpha=2
    c_point=x[0:2,:]+np.dot(R,np.array([-0.7,4*np.asscalar(x[3])]).reshape((2,1)))
    f_point=np.dot(R,np.array([x[4,0]*alpha,x[5,0]*alpha*0.3]).reshape(2,1))
    ax1.arrow(np.asscalar(c_point[0]-f_point[0]),np.asscalar(c_point[1]-f_point[1]),np.asscalar(f_point[0]),np.asscalar(f_point[1]),head_width=0.2)
    fig.savefig('figures/more_planar_pushing_%d.png'%t, dpi=100)
    plt.close()
    return fig


def visualize_push_funnel(s,iterations,di,dj,xmin=-1,xmax=1,ymin=-1,ymax=1,xlabel='x_1',ylabel='x_2',title="interesting plot"):
    """
    di: x dimension
    dj: y dimension
    """
    fig,ax1 = plt.subplots()
    ax1.set_xlabel(xlabel,fontsize=20)
    ax1.set_ylabel(ylabel,fontsize=20)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    list_of_states=s.X[0:s.tree_size[iterations]]
    n=list_of_states[0].x.shape[0]
    vertices_0=vertices_cube(n)
    # Trajectories
    p_list=[]
    for state_X in list_of_states:
        v=np.dot(state_X.G_eps,vertices_0.T)+state_X.x
        v=v[[di,dj],:].T
        hull=ConvexHull(v)
        p_list.append(patches.Polygon((v[hull.vertices,:]), True))
    p=PatchCollection(p_list,color=(0.5,0.5,0.4))
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_title(title,fontsize=20)
    fig.savefig('figures/funnel_pushing_%d.png'%iterations, dpi=100)
    plt.close()
    return fig

def visualize_trajectory_in_funnel(s,t,di,dj,xmin=-1,xmax=1,ymin=-1,ymax=1,xlabel='x_1',ylabel='x_2',title="interesting plot"):
    """
    di: x dimension
    dj: y dimension
    """
    fig,ax1 = plt.subplots()
    ax1.set_xlabel(xlabel,fontsize=20)
    ax1.set_ylabel(ylabel,fontsize=20)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    list_of_states=s.X
    n=list_of_states[0].x.shape[0]
    vertices_0=vertices_cube(n)
    # Trajectories
    p_list=[]
    for state_X in list_of_states:
        v=np.dot(state_X.G_eps,vertices_0.T)+state_X.x
        v=v[[di,dj],:].T
        hull=ConvexHull(v)
        p_list.append(patches.Polygon((v[hull.vertices,:]), True))
    p=PatchCollection(p_list,color=(0.4,0.4,0.5))
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_title(title,fontsize=20)
    trajectory=np.array(s.traj)
    ax1.plot(trajectory[0:t,di,:],trajectory[0:t,dj,:],'r',linewidth=0.5)
    ax1.plot(trajectory[0:t,di,:],trajectory[0:t,dj,:],'ro')
    fig.savefig('figures/funnel_pushing_traj_ytheta_%d.png'%t, dpi=100)
    plt.close()
    return fig

if False:
    for iterations in range(1,473,5):
        fig=visualize_push_funnel(s,iterations,0,1,xmin=0,xmax=10,ymin=-1.5,ymax=1.5,xlabel=r'$x$',ylabel=r'$y$',title="planar pushing tree - %d iterations"%iterations)
    

if False:
    for t in range(1,len(s.traj)):
        fig=visualize_push_box(s,t,xmin=0,xmax=10,ymin=-3,ymax=3,xlabel=r'$x$',ylabel=r'$y$',title="planar pushing")
    #fig.savefig('planar_pushing_bizim_%d.png'%t, dpi=100)
    
if True:
    for t in range(1,len(s.traj)):
        #fig=visualize_trajectory_in_funnel(s,t,0,1,xmin=0,xmax=10,ymin=-2,ymax=2,xlabel=r'$x$',ylabel=r'$y$',title="planar pushing")
        fig=visualize_trajectory_in_funnel(s,t,1,2,xmin=-2,xmax=2,ymin=-0.5,ymax=0.5,xlabel=r'$y$',ylabel=r'$\theta$',title="planar pushing")

    