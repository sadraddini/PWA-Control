#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 12:16:10 2018

@author: sadra
"""

# External imports
import numpy as np
import sys
sys.path.append('../..')
# Internal imports
from main.ana_system import system,state
from main.auxilary_methods import vertices_cube
from main.polytope import polytope,anchor_point



s=system(2,1)
s.modes=[0,1]

s.A[0]=np.array([[1,0.01],[0.1,1]])
s.A[1]=np.array([[1,0.01],[-9.9,1]])

s.B[0]=np.array([[0,0.01]]).T
s.B[1]=np.array([[0,0.01]]).T
s.c[0]=np.array([[0,0]]).T
s.c[1]=np.array([[0,1]]).T

s.H[0]=np.array([[1,0],[0,1],[-1,0],[0,-1]])
s.h[0]=np.array([[0.1,1,0.12,1]]).T   

s.H[1]=np.array([[1,0],[0,1],[-1,0],[0,-1]])
s.h[1]=np.array([[0.12,1,-0.1,1]]).T  


s.F[0]=np.array([[1,-1]]).T
s.f[0]=np.array([[4,4]]).T

s.F[1]=np.array([[1,-1]]).T
s.f[1]=np.array([[4,4]]).T

s.Pi=np.array([[1,0],[0,1],[-1,0],[0,-1]])

s.l[0]=np.array([-0.12,-1]).reshape(2,1)
s.u[0]=np.array([0.1,1]).reshape(2,1)

s.l[1]=np.array([0.1,-1]).reshape(2,1)
s.u[1]=np.array([0.12,1]).reshape(2,1)

s.vertices=vertices_cube(2)

"""
These are polytopes for each mode 
"""
s.mode_polytope={}
for mode in s.modes:
    p=polytope(s.H[mode],s.h[mode])
    p.anchor=anchor_point(p)
    s.mode_polytope[mode]=p


s.weight={}
s.weight[0]=100/12
s.weight[1]=1


    
s.goal=state(np.array([0,0]).reshape(2,1),np.array([[0.,0],[0,0.0]]),0,0,0,10)