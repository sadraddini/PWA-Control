# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 19:29:22 2018

@author: sadra
"""
import numpy as np
import pickle
from pypolycontain.utils.redundancy_reduction import canonical_polytope
from pypolycontain.lib.zonotope import zonotope
from PWA_lib.trajectory.system import linear_cell
from PWA_lib.trajectory.poly_trajectory import polytopic_trajectory_given_modes
from PWA_lib.visualization.visualize import add_tube
import matplotlib.pyplot as plt


q=pickle.load(open("carrot_linear_blocks.pkl","r"))
X_traj=pickle.load(open("X_traj.pkl","r"))
U_traj=pickle.load(open("U_traj.pkl","r"))
T=len(q)
T1=420
T2=480
x0=X_traj[:,T1].reshape(10,1)
goal=zonotope(X_traj[:,T2+1].reshape(10,1),np.eye(10)*0.001)


for t in range(T1,T2-1):
    cell=q[t]
    A,B,c,p=cell.A,cell.B,cell.c,cell.p
    x,u,x_plus=X_traj[:,t].reshape(10,1),U_traj[:,t].reshape(4,1),X_traj[:,t+1].reshape(10,1)
    e=x_plus-np.dot(A,x)-np.dot(B,u)-c
    q[t].c=q[t].c+e
    print e.T,e.shape
    print all(p.h-np.dot(p.H,np.vstack((x,u)))>=-10**-8)

for t in range(T1,T2-1):
    print t
    cell=q[t]
    A,B,c,p=cell.A,cell.B,cell.c,cell.p
    x,u,x_plus=X_traj[:,t].reshape(10,1),U_traj[:,t].reshape(4,1),X_traj[:,t+1].reshape(10,1)
    e=x_plus-np.dot(A,x)-np.dot(B,u)-c
    q[t].c=q[t].c+e
    q[t].p.h=q[t].p.h+10**-5
    print e.T,e.shape
    print all(p.h-np.dot(p.H,np.vstack((x,u)))>=-10**-8)

(x,G)=polytopic_trajectory_given_modes(x0,q[T1:T2],goal,eps=1,order=1)
#assert 1==0



#t=100
#x,u=X_traj[:,t].reshape(10,1),U_traj[:,t].reshape(4,1)
#x_plus=X_traj[:,t+1].reshape(10,1)
#q=find_the_dynamics(cell[t],x,u,A0,B0,c0,H0,h0,x_plus,eps=1,delta_t=0.0025)


#(x,G)=polytopic_trajectory_given_modes(x0,list_of_cells,sys.goal,eps=0.05,order=1)


fig,ax=plt.subplots()
add_tube(ax,x,G,eps=0.0001,list_of_dimensions=[0,1],axis=0)
ax.plot([x[t][0,0] for t in range(T2-T1+1)],[x[t][1,0] for t in range(T2-T1+1)])
pickle.dump(x,open("x_funnel_4.pkl","w"))
pickle.dump(G,open("G_funnel_4.pkl","w"))
#ax.plot([x[t][0,0] for t in range(T+1)],[x[t][1,0] for t in range(T+1)],'+')