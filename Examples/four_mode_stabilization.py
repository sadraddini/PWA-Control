#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 13:48:40 2018

@author: sadra
"""
# External imports
import numpy as np
import sys
sys.path.append('..')

# Internal imports
from main.ana_system import system,state
from main.auxilary_methods import vertices_cube


s=system(2,1)
s.modes=[1,2,3,4]

s.A[1]=np.array([[1,1],[0,1]])
s.A[2]=np.array([[1,1],[0,1]])
s.A[3]=np.array([[1,-1],[0,1]])
s.A[4]=np.array([[1,-1],[0,1]])

s.B[1]=np.array([[1,0.5]]).T
s.B[2]=np.array([[-1,-0.5]]).T
s.B[3]=np.array([[-1,0.5]]).T
s.B[4]=np.array([[1,-0.5]]).T

s.H[1]=np.array([[1,0],[0,1],[-1,0],[0,-1]])
s.h[1]=np.array([[5,5,0,0]]).T   

s.H[2]=np.array([[1,0],[0,1],[-1,0],[0,-1]])
s.h[2]=np.array([[0,0,5,5]]).T  

s.H[3]=np.array([[1,0],[0,1],[-1,0],[0,-1]])
s.h[3]=np.array([[0,5,5,0]]).T 

s.H[4]=np.array([[1,0],[0,1],[-1,0],[0,-1]])
s.h[4]=np.array([[5,0,0,5]]).T 

for i in s.modes:
    s.F[i]=np.array([[1,-1]]).T
    s.f[i]=np.array([[1,1]]).T
    s.c[i]=np.array([[0,0]]).T

s.Pi=np.array([[1,0],[0,1],[-1,0],[0,-1]])

s.l[1]=np.array([0,0]).reshape(2,1)
s.u[1]=np.array([5,5]).reshape(2,1)

s.l[2]=np.array([-5,-5]).reshape(2,1)
s.u[2]=np.array([0,0]).reshape(2,1)

s.l[3]=np.array([-5,0]).reshape(2,1)
s.u[3]=np.array([0,5]).reshape(2,1)

s.l[4]=np.array([0,-5]).reshape(2,1)
s.u[4]=np.array([5,0]).reshape(2,1)

s.vertices=vertices_cube(2)

s.W={}

for i in s.modes:
    s.W[i]=np.array([[0,0],[0,1]])
    
s.weight={}
s.weight[0]=0.2
s.weight[1]=0.2

s.vertices=vertices_cube(2)

s.W={}

for i in s.modes:
    s.W[i]=np.array([[1,1],[1,1]])
    
s.goal=state(np.array([0,0]).reshape(2,1),np.array([[0.1,0],[0,0.1]]),0,0,0,10)