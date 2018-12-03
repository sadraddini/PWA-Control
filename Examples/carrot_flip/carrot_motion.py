#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 07:56:45 2018

@author: sadraddini
"""
from carrot import *


if True:
    X=np.array([0,3,0.0,0,0,0,2.5,1,7,1.57])
    U=np.zeros(4)
    T=500
    X_traj=np.zeros((10,T))
    U_traj=np.zeros((4,T))
    gripper_traj=np.zeros(T)
    alpha=0.17
    for t in range(T):
        gripper="on"
        x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi=X[0:10]
        d=D-R-c*np.sin(psi)-d*np.cos(psi)
        U=np.array([0.0,-3,0.0,0.0])
        if X[7]<=-0*alpha:
            U[3]=-0.5+0.3*(theta>1)
        if X[7]<=-0.10*alpha:
            U[1]=0 # Left finger guard
        if X[7]<0 and c>-0.1:
            U[0]=-0.5
        if psi<-0.01: # Psi guard
            U[3]=0
            U[2]=-8
        if d<=-0.06*alpha: # right finger penetration guard
            U[2]=0
        elif X[7]<=0:
            U[2]=-3
        if y>3.2: # Pickup! Go to right!
            U[0]=0
            if x>=0:
                if d>=-0.1*alpha:
                    U[2]=-5
                else:
                    U[2]=5
        if theta>1.6:
            U=np.array([-40,40,40,-10*(psi)])
        if theta>2.0:
            gripper=0
        else:
            gripper=1
        X_traj[:,t]=X
        U_traj[:,t]=U  
        gripper_traj[t]=gripper
        X=evolve_carrot(X,U,gripper)
        print(t,X.T,"\n")
        #visualize(X,U,gripper)

if True: # Open-loop
    T=500
    X=np.array([0,3,0,0,0,0,2.5,1,7,1.57])
    X_traj_new=np.zeros((10,T))
    for t in range(T):
        X_traj_new[:,t]=X
        U=U_traj[:,t]
        gripper=gripper_traj[t]
        X=evolve_carrot(X,U,gripper)
        visualize(X,U,gripper)    
    e=X_traj-X_traj_new
    plt.plot(range(T),e.T)

if False:
    X=np.array([0,3,0,0,0,0,2.5,1,7,1.57])
    U=np.zeros(4)
    T=600
    X_traj=np.zeros((10,T))
    alpha=0.4
    for t in range(T):
        gripper="on"
        x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi=X[0:10]
        d=D-R-c*np.sin(psi)-d*np.cos(psi)
        U=np.array([0.0,-2,0.0,0.0])
        if X[7]<=0*alpha:
            U[3]=-0.6
        if X[7]<=-0.1*alpha:
            U[1]=0 # Left finger guard
        if X[7]<0 and c>-0.1:
            U[0]=-0.5
        if psi<-0.01: # Psi guard
            U[3]=0
            U[2]=-8
        if d<=-0.11*alpha: # right finger penetration guard
            U[2]=0
        elif X[7]<=0:
            U[2]=-3
        if y>3.2: # Pickup! Go to right!
            U[0]=0
            if x>=0:
                if d>=-0.1*alpha:
                    U[2]=-5
                else:
                    U[2]=5
        if theta>1.6:
            U=np.array([-10,5,10,0])
        if theta>2.3:
            gripper="off"
        X=evolve_carrot(X,U,gripper)
        X_traj[:,t]=X
        print(t,X.T,"\n")
        visualize(X,U,gripper)
        
        
        
if False:
    X=np.array([0,3,0,0,0,0,2,1,7,1.57])
    U=np.zeros(4)
    T=600
    X_traj=np.zeros((10,T))
    alpha=0.3
    for t in range(T):
        gripper="on"
        x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi=X[0:10]
        d=D-R-c*np.sin(psi)-d*np.cos(psi)
        U=np.array([0.0,-2,0.0,0.0])
        if X[7]<=0*alpha:
            U[3]=-0.6
        if X[7]<=-0.1*alpha:
            U[1]=0 # Left finger guard
        if X[7]<0 and c>-0.1:
            U[0]=-0.5
        if psi<-0.01: # Psi guard
            U[3]=0
            U[2]=-8
        if d<=-0.15*alpha: # right finger penetration guard
            U[2]=0
        elif X[7]<=0:
            U[2]=-3
        if y>3.2: # Pickup! Go to right!
            U[0]=0
            if x>=0:
                if d>=-0.1*alpha:
                    U[2]=-5
                else:
                    U[2]=5
        if theta>1.57:
            U=np.array([-10,5,10,0])
        if theta>2.3:
            gripper="off"
        X=evolve_carrot(X,U,gripper)
        X_traj[:,t]=X
        print(t,X.T,"\n")
        visualize(X,U,gripper)
        
if False:
    T2=300
    X=np.array([6,4,1.6,-1,-0.5,0,2,1,7,1.57])
    for t in range(T2):
        x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi=X[0:10]
        d=D-R-c*np.sin(psi)-d*np.cos(psi)
        U=np.array([0,0,0,0])
        X=evolve_carrot(X,U,gripper=0)
        print(t,X.T,"\n")
        visualize(X,U,gripper=0)
        
if False: # Stand up!
    T=100
    X=np.array([0,4,np.pi/2,0,0,0,0,-0.3,2.4,0])
    for t in range(T):
        x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi=X[0:10]
        d=D-R-c*np.sin(psi)-d*np.cos(psi)
        U=np.array([1,0,0,0])
        X=evolve_carrot(X,U)
        print(t,X.T)
        visualize(X,U)
        
if False:
    X=np.array([0,4,np.pi+0.1,0,0,0,200,1,300,1.57])
    U=np.zeros(4)
    T=300
    X_traj=np.zeros((10,T))
    for t in range(T):
        x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi=X[0:10]
        d=D-R-c*np.sin(psi)-d*np.cos(psi)
        U=np.array([-1,-2,-2,-2])
        if d>0:
            U[0]=0
            U[3]=0
        if c<=-0.1:
            U[0]=0
        if psi<-0.01: # Psi guard
            U[3]=0
            U[2]=-3
        if d<=-0.1: # right finger penetration guard
            U[2]=0
        else:
            U[2]=-2
        if X[7]<=-0.1:
            U[1]=0 # Left finger guard
        X=evolve_carrot(X,U)
        X_traj[:,t]=X
        print(t,X.T)
        visualize(X,U)