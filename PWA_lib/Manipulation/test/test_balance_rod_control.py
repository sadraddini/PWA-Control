#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 11:07:22 2019

@author: sadra
"""

import numpy as np
from PWA_lib.trajectory.poly_trajectory import point_trajectory,polytopic_trajectory_given_modes
import pickle

sys=pickle.load(open("sadra_bar.pkl","r"))


gravity=9.8
x0=np.array([0,0,0.0,1,0.00,-1,0.0,1.3,0,0.0]).reshape(10,1)
#u0=np.array(([0,0,0,0,gravity/2,0,gravity/2,0])).reshape(8,1)
u0=np.array(([0,0,0,0,0,0,gravity/2,0])).reshape(8,1)

T=60
(x_n,u,delta_PWA,mu,flag)=point_trajectory(sys,x0,[sys.goal],T,eps=0,optimize_controls_indices=[])

list_of_cells=[]
list_of_modes={}
for t in range(T):
    mode=tuple([i for n in sys.list_of_sum_indices for i in sys.list_of_modes[n] if abs(delta_PWA[t,n,i]-1)<10**-6]) 
    list_of_modes[t]=mode
    list_of_cells.append(sys.cell[mode])

#A=A0+sys.A[1,"ST"]+sys.A[2,"ST"]
#B=B0+sys.B[1,"ST"]+sys.B[2,"ST"]
#c=c0+sys.c[1,"ST"]+sys.c[2,"ST"]
#np.dot(A,x0)
#x_future=np.dot(A,x0)+np.dot(B,u0)+c
#
#xu0=np.vstack((x0,u0))
#sys.C[1,"ST"].if_inside(xu0)


import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

plt.plot([x_n[t][2] for t in range(T)])

def animate(X,ax1):
    x,y,theta,x_1,y_1,x_2,y_2=X[0:7]
    L=1.5
    b=0.3
    # x
    x_left=x-L*np.cos(theta)
    x_left_bar=x-L*np.cos(theta)-b*np.sin(theta)
    x_right=x+L*np.cos(theta)
    x_right_bar=x+L*np.cos(theta)-b*np.sin(theta)
    # y
    y_left=y-L*np.sin(theta)
    y_left_bar=y_left+b*np.cos(theta)
    y_right=y+L*np.sin(theta)
    y_right_bar=y_right+b*np.cos(theta)
    # Good now plot
    ax1.set_xlabel("x",fontsize=20)
    ax1.set_ylabel("y",fontsize=20)
    ax1.set_xlim([-2*L,2*L])
    ax1.set_ylim([-0.4,0.8])
    fig.gca().set_aspect('equal')
    bar=[patches.Polygon(np.array([[x_left,x_left_bar,x_right_bar,x_right],[y_left,y_left_bar,y_right_bar,y_right]]).reshape(2,4).T, True)]
    ax1.add_collection(PatchCollection(bar,color=(0.8,0.3,0.4),alpha=0.8,edgecolor=(0,0,0)))
#    xf_1=x+x_1*np.cos(theta)+y_1*np.sin(theta)
#    xf_2=x+x_2*np.cos(theta)+y_2*np.sin(theta)
#    yf_1=y+x_1*np.sin(theta)-y_1*np.cos(theta)-0.07
#    yf_2=y+x_2*np.sin(theta)-y_2*np.cos(theta)-0.07
#    ax1.plot([xf_1,xf_2],[yf_1,yf_2],'^',linewidth=3,markersize=10)
    ax1.plot([x_1,x_2],[y_1,y_2],'^',linewidth=3,markersize=10)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.5)


mode_name={"NC":"No Contact", "SP": "Contact >", "SN":"Contact <", "ST":"Contact -"}
    
for t in range(T):
    fig,ax = plt.subplots()
    ax.set_title("Balancing t=%d \n Right finger: %s \n Left Finger: %s"%(t,mode_name[list_of_modes[t][1]],mode_name[list_of_modes[t][2]]))
    animate(x_n[t],ax)
    if t==0:
        ax.arrow(-1.8, 0.1, 0.3, 0.0, head_width=0.3, head_length=0.3, linewidth=10, fc='k', ec='k')
    fig.savefig('figures/bar_%d.png'%t, dpi=100)
    
(x,u,G,theta)=polytopic_trajectory_given_modes(x0,list_of_cells,sys.goal,eps=1,order=1)

from PWA_lib.visualization.visualize import add_tube
fig,ax=plt.subplots()
ax.set_xlim([-0.5,0.5])
ax.set_ylim([-1,1])
list_of_dimensions=[0,1]
i,j=list_of_dimensions
add_tube(ax,x,G,eps=0.0001,list_of_dimensions=list_of_dimensions)
ax.plot([x[t][i,0] for t in range(T+1)],[x[t][j,0] for t in range(T+1)])
#ax.plot([x[t][0,0] for t in range(T+1)],[x[t][1,0] for t in range(T+1)],'+')
ax.plot([0],[0],'o')
plt.plot([x_n[t][i,0] for t in range(T+1)],[x_n[t][j,0] for t in range(T+1)])
plt.plot([x_n[t][i,0] for t in range(T+1)],[x_n[t][j,0] for t in range(T+1)],'+')
#plt.plot([-0.05],[0.5],'o')

from PWA_lib.polytree.tree import tree
mytree=tree(sys)
mytree.fill="continous"
mytree.add_branch(x,u,G,theta,sys.goal) 
mytree.visualize(axis_limit=[-0.5,0.5,-1,1])