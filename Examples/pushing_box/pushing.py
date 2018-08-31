#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 10:42:55 2018

@author: sadra
"""

# External imports
import numpy as np
import sys
sys.path.append('../..')

# Internal imports
from main.ana_system import system,state
from main.auxilary_methods import vertices_cube,PI


s=system(6,2)
s.modes=[1,2,3,4]

dt=0.1 # time-step
alpha=1 # translational coefficient
beta=1 # angular coefficient
gamma=1 # sliding coefficient
l=1 # arm for fT
mu=0.5 # friction coefficient

rx_nominal=1

"""
 Mode 1: No Contact x[4]<=0
"""

s.A[1]=np.zeros((6,6))
for row in range(4):
    s.A[1][row,row]=1
s.c[1]=np.zeros((6,1))

xmax=10
ymax=10
theta_max=np.pi
bmax=1
rx_max=10
ry_max=rx_max

s.H[1]=np.vstack((np.eye(6),-np.eye(6)))
s.h[1]=np.array([xmax,ymax,theta_max,bmax,0,ry_max,xmax,ymax,theta_max,bmax,rx_max,ry_max]).reshape(12,1)

""" 
    Mode 2: Sticking, limited change in b, limited change in theta, limited change in rx
    rx>0, mu*rx<=ry<=mu*rx
    |b|<b_epsilon
    |theta|<=theta_eps
"""
s.A[2]=np.zeros((6,6))
for row in range(4):
    s.A[2][row,row]=1
s.A[2][0,4]=alpha*dt
s.A[2][1,2]=alpha*dt*rx_nominal
s.A[2][1,5]=alpha*dt
s.A[2][2,3]=-beta*dt*rx_nominal
s.A[2][2,5]=-beta*dt*l

s.c[2]=np.array([0,0,0,0,0,0]).reshape(6,1)

theta_eps=0.3
b_eps=0.1
s.H[2]=np.vstack((np.eye(6),-np.eye(6)))
s.h[2]=np.array([xmax,ymax,theta_eps,b_eps,rx_nominal*1.2,ry_max,xmax,ymax,theta_eps,b_eps,rx_nominal*0.8,ry_max]).reshape(12,1)

s.H[2]=np.vstack((s.H[2],np.array([0,0,0,0,-mu,-1]).reshape(1,6),np.array([0,0,0,0,-mu,1]).reshape(1,6)))
s.h[2]=np.vstack((s.h[2],np.zeros((2,1))))

"""
    Mode 3: Sliding up ry>=mu*rx
        y[0]=x[0]+x[4]*alpha*dt*np.cos(x[2])-mu*x[4]*alpha*dt*np.sin(x[3])
        y[1]=x[1]+x[4]*alpha*dt*np.sin(x[2])+mu*x[4]*alpha*dt*np.cos(x[3])
        y[2]=x[2]-beta*l*dt*mu*x[4]-beta*x[3]*x[4]*dt
"""
s.A[3]=np.zeros((6,6))
for row in range(4):
    s.A[3][row,row]=1
    
s.A[3][0,4]=alpha*dt
s.A[3][1,2]=alpha*dt*rx_nominal
s.A[3][1,4]=mu*alpha*dt*rx_nominal
s.A[3][2,3]=-beta*dt*rx_nominal
s.A[3][2,4]=-beta*l*dt*mu
s.A[3][3,5]=gamma*dt

s.c[3]=np.array([0,0,0,0,0,0]).reshape(6,1)

s.H[3]=np.vstack((np.eye(6),-np.eye(6)))
s.h[3]=np.array([xmax,ymax,theta_eps,b_eps,rx_nominal*1.2,ry_max,xmax,ymax,theta_eps,b_eps,rx_nominal*0.8,ry_max]).reshape(12,1)

s.H[3]=np.vstack((s.H[2],np.array([0,0,0,0,mu,-1]).reshape(1,6)))
s.h[3]=np.vstack((s.h[2],np.zeros((1,1))))


"""
    Mode 4: Sliding down ry<=-mu*rx
        y[0]=x[0]+x[4]*alpha*dt*np.cos(x[2])-mu*x[4]*alpha*dt*np.sin(x[3])
        y[1]=x[1]+x[4]*alpha*dt*np.sin(x[2])+mu*x[4]*alpha*dt*np.cos(x[3])
        y[2]=x[2]-beta*l*dt*mu*x[4]-beta*x[3]*x[4]*dt
"""
s.A[3]=np.zeros((6,6))
for row in range(4):
    s.A[3][row,row]=1
    
s.A[4][0,4]=alpha*dt
s.A[4][1,2]=alpha*dt*rx_nominal
s.A[4][1,4]=-mu*alpha*dt*rx_nominal
s.A[4][2,3]=-beta*dt*rx_nominal
s.A[4][2,4]=beta*l*dt*mu
s.A[4][3,5]=gamma*dt

s.c[4]=np.array([0,0,0,0,0,0]).reshape(6,1)

s.H[4]=np.vstack((np.eye(6),-np.eye(6)))
s.h[4]=np.array([xmax,ymax,theta_eps,b_eps,rx_nominal*1.2,ry_max,xmax,ymax,theta_eps,b_eps,rx_nominal*0.8,ry_max]).reshape(12,1)

s.H[4]=np.vstack((s.H[2],np.array([0,0,0,0,mu,1]).reshape(1,6)))
s.h[4]=np.vstack((s.h[2],np.zeros((1,1))))

"""
Control limits and bounding-boxes
"""
    
for mode in s.modes:
    s.B[mode]=np.zeros((6,2))
    s.B[mode][4,0]=1
    s.B[mode][5,1]=1

for i in s.modes:
    s.F[i]=PI(2)
    s.f[i]=10*np.ones((4,1))

s.Pi=PI[6]

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