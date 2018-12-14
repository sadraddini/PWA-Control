# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 11:14:13 2018

@author: sadra
"""

import numpy as np

imported=True

if imported==False:
    import pickle
    
    cell=pickle.load(open("cell.pkl","r"))
    X_traj=pickle.load(open("X_traj.pkl","r"))
    U_traj=pickle.load(open("U_traj.pkl","r"))
    
t=200
x=X_traj[:,t].reshape(10,1)
u=U_traj[:,t].reshape(4,1)
x_plus=X_traj[:,t+1]

d={}
for i in [0,1,2]:
    for contact in ["ground","left","right"]:
        d[i,contact]=cell[t][0,contact].p.h-np.dot(cell[t][0,contact].p.H,np.vstack((x,u)))