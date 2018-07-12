#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 15:36:54 2018

@author: sadra
"""
# External imports
import numpy as np

# Internal imports
from ana_system import system,state
from auxilary_methods import vertices_cube

# Code:
s=system(2,1)
s.modes=[0,1]
dt=0.05

K=40
g=9.8
xmax=1.2
vmax=5
umax=3
dmax=1
y_p=1 # I want to get into this height


s.A[0]=np.array([[1,dt],[0,1]])
s.A[1]=np.array([[1,dt],[-K*dt,1]])
s.A[1]=np.array([[0,0],[0,-1]])

s.B[0]=np.array([[0,dt]]).T
s.B[1]=np.array([[0,dt]]).T
s.c[0]=np.array([[0,-g*dt]]).T
s.c[1]=np.array([[0,-g*dt]]).T
s.c[1]=np.array([[0.01*dt,0]]).T
s.B[1]=np.array([[0,0]]).T

s.H[0]=np.array([[1,0],[0,1],[-1,0],[0,-1]])
s.h[0]=np.array([[xmax,vmax,0,vmax]]).T   

s.H[1]=np.array([[1,0],[0,1],[-1,0],[0,-1]])
s.h[1]=np.array([[0,vmax,xmax,vmax]]).T

s.F[0]=np.array([[1,-1]]).T
s.f[0]=np.array([[umax,umax]]).T

s.F[1]=np.array([[1,-1]]).T
s.f[1]=np.array([[umax,umax]]).T


s.Pi=np.array([[1,0],[0,1],[-1,0],[0,-1]])

s.l[0]=np.array([0,-vmax]).reshape(2,1)
s.u[0]=np.array([y_p,vmax]).reshape(2,1)

s.l[1]=np.array([-dmax,-vmax]).reshape(2,1)
s.u[1]=np.array([0,vmax]).reshape(2,1)

s.vertices=vertices_cube(2)

s.W={}

for i in s.modes:
    s.W[i]=np.array([[0,0],[0,1]])

s.weight={}
s.weight[0]=16
s.weight[1]=1

s.goal=state(np.array([(xmax+y_p)/2,0]).reshape(2,1),np.array([[(xmax-y_p)/2,0],[0,1]]),0,0,0,10)