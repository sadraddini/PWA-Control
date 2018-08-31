#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 11:14:36 2018

@author: sadra
"""
#%%
x=np.zeros(6)
for t in range(100):
    print(t," ",x.T)
    Ky=1
    Ktheta=3
    u=np.array([1,Ky*x[1]+Ktheta*x[2]])
    x=pushing_evolve(x,u)+(np.random.random(6)-0.5)*0.6*(t<30)
#%%    
x=np.zeros(6)
for t in range(42):
    print(t," ",x.T)
    Ky=-1
    Ktheta=40
    u=np.array([1,0])
    x=pushing_evolve(x,u)+(np.array(2)-0.5)*0.0000
#%%  Controllability
    mode=3
    C=s.B[mode]
    for n in range(1,6):
        C=np.hstack((C,np.dot(np.linalg.matrix_power(s.A[mode], n),s.B[mode])))