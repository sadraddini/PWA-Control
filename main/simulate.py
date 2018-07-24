#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 11:09:51 2018

@author: sadra
"""

import numpy as np

from main.auxilary_methods import find_mode,inside_X
from main.controller import control_vanilla

def simulate_vanilla(s,x):
    t=0
    s.traj=[]
    s.control_traj=[]
    while t<100:
        t+=1
        s.traj.append(x)
        print("state:",x.T)
        u=control_vanilla(s,x)
        if u=="END":
            print("END Control")
            return
        print("control:",u.T)
        x=evolve(s,x,u)
    
def evolve(s,x,u):
    i=find_mode(s,x)
#    print("x=",x.T,"u=",u,"i=",i,"\n")
    return np.dot(s.A[i],x)+np.dot(s.B[i],u)+s.c[i]#+(np.random.random((2,1))-0.5)*0.00001

def simulate_0(s,x,T):
    for t in range(T):
        u=np.zeros((s.m,1))
        x_next=evolve(s,x,u)
        if inside_X(s,x_next)==False:
            return (x,T)
        else:
            x=x_next
    return (x,T)