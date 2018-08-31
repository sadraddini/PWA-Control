#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 10:44:58 2018

@author: sadra
"""

# External imports
import numpy as np
import sys
sys.path.append('../..')

# Internal imports
from main.ana_system import system,state
from main.auxilary_methods import vertices_cube

def pushing_evolve(x,u):
    """
    x[0]= x: world frame
    x[1]= y: world frame
    x[2]= theta: world frame
    x[3]= b: arm of contact point from center-line, body frame
    x[4]= rx: normal force, body frame
    x[5]= ry: tangential force, body frame
    """
    dt=0.01
    alpha=1
    beta=1
    gamma=1
    l=1
    mu=0.5
    y=np.empty((6))
    if x[4]>=0 and abs(x[5])<=mu*x[4]: # Sticking
        print("sticking")
        y[0]=x[0]+x[4]*alpha*dt*np.cos(x[2])-x[5]*alpha*dt*np.sin(x[2])
        y[1]=x[1]+x[4]*alpha*dt*np.sin(x[2])+x[5]*alpha*dt*np.cos(x[2])
        y[2]=x[2]-beta*l*dt*x[5]-beta*x[3]*x[4]*dt
        y[3]=x[3]
        y[4]=u[0]
        y[5]=u[1]
    elif x[4]>=0 and x[5]>mu*x[4]: # Sliding up
        print("sliding up")
        y[0]=x[0]+x[4]*alpha*dt*np.cos(x[2])-mu*x[4]*alpha*dt*np.sin(x[2])
        y[1]=x[1]+x[4]*alpha*dt*np.sin(x[2])+mu*x[4]*alpha*dt*np.cos(x[2])
        y[2]=x[2]-beta*l*dt*mu*x[4]-beta*x[3]*x[4]*dt
        y[3]=x[3]+gamma*x[5]*dt
        y[4]=u[0]
        y[5]=u[1]
    elif x[4]>=0 and x[5]<-mu*x[4]: # Sliding down
        print("sliding down")
        y[0]=x[0]+x[4]*alpha*dt*np.cos(x[2])+mu*x[4]*alpha*dt*np.sin(x[2])
        y[1]=x[1]+x[4]*alpha*dt*np.sin(x[2])-mu*x[4]*alpha*dt*np.cos(x[2])
        y[2]=x[2]+beta*l*dt*mu*x[4]-beta*x[3]*x[4]*dt
        y[3]=x[3]+gamma*x[5]*dt
        y[4]=u[0]
        y[5]=u[1]  
    elif x[4]<0: # No Contact
        y[0]=x[0]
        y[1]=x[1]
        y[2]=x[2]
        y[3]=x[3]
        y[4]=u[0]
        y[5]=u[1]
    return y