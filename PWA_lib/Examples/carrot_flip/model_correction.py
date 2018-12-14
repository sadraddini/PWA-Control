# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 11:44:18 2018

@author: sadra
"""

from carrot import *
from pypolycontain.utils.redundancy_reduction import canonical_polytope
from pypolycontain.lib.polytope import polytope
from pypolycontain.lib.zonotope import zonotope

from PWA_lib.trajectory.system import linear_cell
from PWA_lib.trajectory.poly_trajectory import polytopic_trajectory_given_modes


import numpy as np
import pickle
cell=pickle.load(open("cell.pkl","r"))  
X_traj=pickle.load(open("X_traj.pkl","r"))
U_traj=pickle.load(open("U_traj.pkl","r"))

A0=np.eye(10)
B0=np.vstack((np.zeros([6,4]),h*np.eye(4)))
c0=np.zeros((10,1))
c0[4,0]=-g*h
H0=np.vstack((np.eye(14),-np.eye(14)))
h0=np.ones((28,1))*10


T=len(cell)

#for t in range(T-1):
#    d_numerical=X_traj[:,t+1]-X_traj[:,t]
#    find_the_dynamics(cell[t],x,u,A0,B0,c0,H0,h0,eps,h=0.0025)


def point_in_polytope(p,x):
    """
    arguments point x and polyhedron p
    return x in p
    """
    e=p.h-np.dot(p.H,x)
    e.reshape(p.h.shape[0])
    return all(e>=0)
    
def find_the_dynamics(pwa_cell,x,u,A0,B0,c0,H0,h0,x_plus,eps=1,delta_t=0.0025):
    A=A0
    B=B0
    c=c0
    H=H0
    h=h0
    xu=np.vstack((x,u))
    for contact in ["ground","left","right"]:
        for i in range(3):
            p=pwa_cell[i,contact].p
#            print p.H.shape,p.h.shape,xu.shape
            flag=point_in_polytope(p,xu)
#            print p.H,p.h,p.h-np.dot(p.H,xu)
            print p.h-np.dot(p.H,xu)
            print contact,i,flag
            A=A+pwa_cell[i,contact].A*flag*delta_t    
            B=B+pwa_cell[i,contact].B*flag*delta_t
            c=c+pwa_cell[i,contact].c*flag*delta_t
            if flag==True:
                H=np.vstack((H,p.H))
                h=np.vstack((h,p.h))
#    print c0.T
    (H,h)=canonical_polytope(H,h)
    e=x_plus-np.dot(A,x)-np.dot(B,u)-c
    c=c+e
    print e.T,(x_plus-np.dot(A,x)-np.dot(B,u)-c).T
    return linear_cell(A,B,c,polytope(H,h))
    
q=[0]*T
for t in range(T):
    print t,"*"*50
    x,u=X_traj[:,t].reshape(10,1),U_traj[:,t].reshape(4,1)
    x_plus=X_traj[:,t+1].reshape(10,1)
    q[t]=find_the_dynamics(cell[t],x,u,A0,B0,c0,H0,h0,x_plus,eps=1,delta_t=0.0025)
    
pickle.dump(q,open("carrot_linear_blocks.pkl","w"))
    