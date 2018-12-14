# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 19:55:39 2018

@author: sadra
"""
from carrot_PWA_generation import PWA
import pickle
X_traj=pickle.load(open("X_traj.pkl","r"))
U_traj=pickle.load(open("U_traj.pkl","r"))

cell={}
T=490
for t in range(T):
    print t
    x=X_traj[:,t].reshape(10,1)
    u=U_traj[:,t].reshape(4,1)
    cell[t]=PWA(x,u)
    
pickle.dump(cell,open("cell.pkl","w"))