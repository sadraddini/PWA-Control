#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 10:59:22 2018

@author: sadra
"""
import numpy as np

from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import pylab

from main.ana_system import states_time_order,states_cost_order

def visualize_set_tube(list_of_states,xmin=-1,xmax=1,ymin=-1,ymax=1,tube_size=0.01):
    from matplotlib.collections import PatchCollection
    ax1 = plt.subplot(111)
    plt.figure(figsize=(20,20),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel('Height')
    ax1.set_ylabel('Velocity')
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    ax1.autoscale_view()
    vertices_0=np.array([[1,1,-1,-1],[1,-1,-1,1]])
    # Trajectories
    p_list=[]
    for state in list_of_states:
#        vertices=vertices_cube(state.G.shape[1]).T
        v=np.dot(state.G,vertices_0)
        v=Minkowski_hull(v.T,(vertices_0.T*tube_size)).T
        p_list.append(patches.Polygon(v.T+state.x.T, True))
    p=PatchCollection(p_list, alpha=0.4,color=(0.2,0.2,0.7))
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_xlabel('Height')
    ax1.set_ylabel('Velocity')
    plt.show()
    
def visualize_set_tube_simulation(s,xmin=-1,xmax=1,ymin=-1,ymax=1,xlabel='x_1',ylabel='x_2',title="interesting plot"):
    ax1 = plt.subplot(111)
    plt.figure(figsize=(20,20),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel(xlabel,fontsize=20)
    ax1.set_ylabel(ylabel,fontsize=20)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    vertices_0=np.array([[1,-1,-1,1],[1,1,-1,-1]])
    # Trajectories
    p_list=[]
    for state in s.X:
#        vertices=vertices_cube(state.G.shape[1]).T
        v=np.dot(state.G_eps,vertices_0)
        p_list.append(patches.Polygon(v.T+state.x.T, True))
    p=PatchCollection(p_list,color=(0.4,0.4,0.4))
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    trajectory=np.array(s.traj)
    ax1.plot(trajectory[:,0,:],trajectory[:,1,:],'r',linewidth=0.5)
    ax1.plot(trajectory[:,0,:],trajectory[:,1,:],'ro')
    ax1.plot([0.1,0.1],[-1,1],'black')
    plt.show()
    


def visualize_X_eps_time(s,xmin=-1,xmax=1,ymin=-1,ymax=1,xlabel='x_1',ylabel='x_2',title="interesting plot"):
    ax1 = plt.subplot(111)
    plt.figure(figsize=(10,10),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel(xlabel,fontsize=20)
    ax1.set_ylabel(ylabel,fontsize=20)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    vertices_0=np.array([[1,-1,-1,1],[1,1,-1,-1]])
    # Trajectories
    p_list=[]
    STATES=states_time_order(s)
    for state in STATES[::-1]:
#        vertices=vertices_cube(state.G.shape[1]).T
        v=np.dot(state.G_eps,vertices_0)
        p_list.append(patches.Polygon(v.T+state.x.T, True))
    max_T=max([state.time_to_go for state in STATES])
    p=PatchCollection(p_list,color=[(state.time_to_go/max_T,1-state.time_to_go/max_T,0) for state in STATES[::-1]])
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
#    ax1.set_title(title)￼
    plt.show()
    return plt

def visualize_X_eps_cost(s,xmin=-1,xmax=1,ymin=-1,ymax=1,xlabel='x_1',ylabel='x_2',title="interesting plot"):
    ax1 = plt.subplot(111)
    plt.figure(figsize=(20,20),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)
    vertices_0=np.array([[1,-1,-1,1],[1,1,-1,-1]])
    # Trajectories
    p_list=[]
    STATES=states_cost_order(s)
    for state in STATES[::-1]:
#        vertices=vertices_cube(state.G.shape[1]).T
        v=np.dot(state.G_eps,vertices_0)
        p_list.append(patches.Polygon(v.T+state.x.T, True))
    max_J=max([state.cost_to_go for state in STATES])
    p=PatchCollection(p_list,color=[(state.cost_to_go/max_J,1-state.cost_to_go/max_J,0) for state in STATES[::-1]])
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_title(title)
    plt.show()
    return plt

def vertices_hull_eps(s,list_of_states):
    """
    Convex hull computation of list of states
    """
    PV=np.empty((0,s.n))
    for p in list_of_states:
        PV=np.vstack((PV,p.vertices_eps+p.x.T))
    hull=ConvexHull(PV)
    return PV[hull.vertices,:]

def visualize_X_time_hull_eps(s,xmin=-1,xmax=1,ymin=-1,ymax=1,xlabel='x_1',ylabel='x_2',title="interesting plot"):
    ax1 = plt.subplot(111)
    plt.figure(figsize=(20,20),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel(xlabel,fontsize=20)
    ax1.set_ylabel(ylabel,fontsize=20)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    STATES=states_time_order(s)
    max_T=max([state.time_to_go for state in STATES])
    p_list=[]
    for state in STATES[::-1]:
        if state.mode==state.successor[0].mode:
            p_list.append(patches.Polygon(vertices_hull_eps(s,[state]+[state.successor[0]]), True))
        else:
            p_list.append(patches.Polygon(vertices_hull_eps(s,[state]), True))
    p=PatchCollection(p_list,color=[(state.time_to_go/max_T,1-state.time_to_go/max_T,0) for state in STATES[::-1]])
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_title(title)
    ax1.plot([0.1,0.1],[-1,1],'black')
    plt.show()
    return plt

def visualize_X_time_hull_eps_simulation(s,xmin=-1,xmax=1,ymin=-1,ymax=1,xlabel='x_1',ylabel='x_2',title="interesting plot"):
    ax1 = plt.subplot(111)
    plt.figure(figsize=(20,20),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel(xlabel,fontsize=20)
    ax1.set_ylabel(ylabel,fontsize=20)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    STATES=states_time_order(s)
    max_T=max([state.time_to_go for state in STATES])
    p_list=[]
    for state in STATES[::-1]:
        if state.mode==state.successor[0].mode:
            p_list.append(patches.Polygon(vertices_hull_eps(s,[state]+[state.successor[0]]), True))
        else:
            p_list.append(patches.Polygon(vertices_hull_eps(s,[state]), True))
    p=PatchCollection(p_list,color=[(state.time_to_go/max_T,1-state.time_to_go/max_T,0) for state in STATES[::-1]])
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_title(title)
    plt.show()
    return plt



def visualize_X_hull_eps(s,list_of_states,xmin=-1,xmax=1,ymin=-1,ymax=1,xlabel='x_1',ylabel='x_2',title="interesting plot"):
    ax1 = plt.subplot(111)
    #fig=plt.figure(figsize=(20,20),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    p_list=[]
    for state in list_of_states:
        p_list.append(patches.Polygon(vertices_hull_eps(s,[state]+[state.successor[0]]), True))
    p=PatchCollection(p_list,color=(0.4,0.4,0.4))
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_title(title)
    return plt

    
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

def visualize_subset_tree(s,iterations,xmin=-1,xmax=1,ymin=-1,ymax=1,xlabel='x_1',ylabel='x_2'):
    ax1 = plt.subplot(111)
    plt.figure(figsize=(20,20),dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel(xlabel,fontsize=24)
    ax1.set_ylabel(ylabel,fontsize=24)
    ax1.legend(str(iterations)+"iterations")
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    vertices_0=np.array([[1,-1,-1,1],[1,1,-1,-1]])
    # Trajectories
    p_list=[]
    for state in s.X[0:s.tree_size[iterations]]:
#        vertices=vertices_cube(state.G.shape[1]).T
        v=np.dot(state.G_eps,vertices_0)
        p_list.append(patches.Polygon(v.T+state.x.T, True))
    p=PatchCollection(p_list,color=(0.3,0.3,0.3))
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_title(str(iterations)+" iterations,  "+str(s.tree_size[iterations])+" polytopes",fontsize=24)
    ax1.plot([0,0,0],[-8,0,8],color=(0,0,0),linewidth=1)
    plt.show()
    return plt