#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 15:27:23 2018

@author: sadra
"""

# Primary imports
import numpy as np
from gurobipy import Model,GRB,LinExpr,QuadExpr
from random import choice as rchoice
from random import random

import sys
sys.path.append('..')

# Secondary imports
from main.auxilary_methods import find_mode,valuation,mode_sequence
from main.ana_system import state,cost_state
from main.trajectory import subset_MILP

def trajectory_model(s,T):
    model=Model("polytopic trajectory of PWA systems")
    x={}
    u={}
    theta={}
    z={}
    G_bound=100
    # Mode 1:
    for t in range(T):
        x[t]=np.empty((s.n,1),dtype='object') # n*1
        u[t]=np.empty((s.m,1),dtype='object') # m*1
        theta[t]=np.empty((s.m,s.n),dtype='object') # n*m
        for row in range(s.n):
            x[t][row,0]=model.addVar(lb=-G_bound,ub=G_bound)
        for row in range(s.m):
            u[t][row,0]=model.addVar(lb=-G_bound,ub=G_bound)
        for row in range(s.m):
            for column in range(s.n):
                theta[t][row,column]=model.addVar(lb=-G_bound,ub=G_bound)   
    for t in range(T+1):
        for i in s.modes:
            z[t,i]=model.addVar(vtype=GRB.BINARY)
    x[T]=np.empty((s.n,1),dtype='object') # Final state in Mode i
    for row in range(s.n):
        x[T][row,0]=model.addVar(lb=-G_bound,ub=G_bound)
    G={}
    for t in range(T+1):
        G[t]=np.empty((s.n,s.n),dtype='object')
        for row in range(s.n):
            for column in range(s.n):
                G[t][row,column]=model.addVar(lb=-G_bound,ub=G_bound)
    model.update()
    # Trajectory Constraints:
    # Mode i:    x={}
    u={}
    theta={}
    z={}
    G_bound=100
    # Mode 1:
    for t in range(T):
        x[t]=np.empty((s.n,1),dtype='object') # n*1
        u[t]=np.empty((s.m,1),dtype='object') # m*1
        theta[t]=np.empty((s.m,s.n),dtype='object') # n*m
        for row in range(s.n):
            x[t][row,0]=model.addVar(lb=-G_bound,ub=G_bound)
        for row in range(s.m):
            u[t][row,0]=model.addVar(lb=-G_bound,ub=G_bound)
        for row in range(s.m):
            for column in range(s.n):
                theta[t][row,column]=model.addVar(lb=-G_bound,ub=G_bound)   
    for t in range(T+1):
        for i in s.modes:
            z[t,i]=model.addVar(vtype=GRB.BINARY)
    x[T]=np.empty((s.n,1),dtype='object') # Final state in Mode i
    for row in range(s.n):
        x[T][row,0]=model.addVar(lb=-G_bound,ub=G_bound)
    G={}
    for t in range(T+1):
        G[t]=np.empty((s.n,s.n),dtype='object')
        for row in range(s.n):
            for column in range(s.n):
                G[t][row,column]=model.addVar(lb=-G_bound,ub=G_bound)
    model.update()
    # Trajectory Constraints:
    # Mode i:
    bigM=10
    for i in s.modes:
        for t in range(T):
            for row in range(s.n):
                Ax=LinExpr()
                for k in range(s.n):
                    Ax.add(s.A[i][row,k]*x[t][k,0])
                for k in range(s.m):
                    Ax.add(s.B[i][row,k]*u[t][k,0])
                model.addConstr(x[t+1][row,0]<=Ax+s.c[i][row]+bigM-bigM*z[t,i])
                model.addConstr(x[t+1][row,0]>=Ax+s.c[i][row]-bigM+bigM*z[t,i])
    # Generator Dynamics Constraints:
    for i in s.modes:
        for t in range(T):
            for row in range(s.n):
                for column in range(s.n):
                    AG=LinExpr()
                    for k in range(s.n):
                        AG.add(s.A[i][row,k]*G[t][k,column])
                    for k in range(s.m):
                        AG.add(s.B[i][row,k]*theta[t][k,column])
                    model.addConstr(G[t+1][row,column]<=AG+bigM-bigM*z[t,i])
                    model.addConstr(G[t+1][row,column]>=AG-bigM+bigM*z[t,i])
    # Constraints of modes:
    for t in range(T+1):
        sum_z=LinExpr()
        for i in s.modes:
            sum_z.add(z[t,i])
        model.addConstr(sum_z==1)
    # Constraints of mode subsets
    for t in range(T):
        for i in s.modes:
            subset_MILP(model,G[t],s.Pi,s.H[i],s.h[i],x[t],z[t,i])
            subset_MILP(model,theta[t],s.Pi,s.F[i],s.f[i],u[t],z[t,i])   
    # Terminal Constraint 
    bigM=100
    for i in s.modes:
        for t in range(T):
            for row in range(s.n):
                Ax=LinExpr()
                for k in range(s.n):
                    Ax.add(s.A[i][row,k]*x[t][k,0])
                for k in range(s.m):
                    Ax.add(s.B[i][row,k]*u[t][k,0])
                model.addConstr(x[t+1][row,0]<=Ax+s.c[i][row]+bigM-bigM*z[t,i])
                model.addConstr(x[t+1][row,0]>=Ax+s.c[i][row]-bigM+bigM*z[t,i])
    # Generator Dynamics Constraints:
    for i in s.modes:
        for t in range(T):
            for row in range(s.n):
                for column in range(s.n):
                    AG=LinExpr()
                    for k in range(s.n):
                        AG.add(s.A[i][row,k]*G[t][k,column])
                    for k in range(s.m):
                        AG.add(s.B[i][row,k]*theta[t][k,column])
                    model.addConstr(G[t+1][row,column]<=AG+bigM-bigM*z[t,i])
                    model.addConstr(G[t+1][row,column]>=AG-bigM+bigM*z[t,i])
    # Constraints of modes:
    for t in range(T+1):
        sum_z=LinExpr()
        for i in s.modes:
            sum_z.add(z[t,i])
        model.addConstr(sum_z==1)
    # Constraints of mode subsets
    for t in range(T):
        for i in s.modes:
            subset_MILP(model,G[t],s.Pi,s.H[i],s.h[i],x[t],z[t,i])
            subset_MILP(model,theta[t],s.Pi,s.F[i],s.f[i],u[t],z[t,i])   
    s.library[T]=(model,x,u,G,theta,z)