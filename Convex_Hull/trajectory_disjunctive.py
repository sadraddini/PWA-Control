#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 14:38:17 2018

@author: sadra
"""

### Internal imports
import sys
sys.path.append('..')

# Primary imports
import numpy as np
from gurobipy import Model,GRB,LinExpr,QuadExpr
from random import choice as rchoice
from random import random

# Secondary imports
from main.auxilary_methods import find_mode,valuation,mode_sequence
from main.ana_system import state

from main.trajectory import subset_MILP

def trajectory_to_p_set_big_M(s,x0,T,list_of_polytopes):
    model=Model("trajectory to a list of polytopes using big-M method")
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
    # set objective
    J_area=LinExpr()
    d_min=model.addVar(lb=0.0001)
    beta=10**2 # Weight of infinity norm
    model.update()
    for row in range(s.n):
        for column in range(s.n):
            if coin<0.1:
                if row<column:
                    model.addConstr(G[0][row,column]==0)
            elif coin>0.9:
                if row>column:
                    model.addConstr(G[0][row,column]==0)                
            if row==column:
                model.addConstr(G[0][row,column]>=d_min/s.weight[row])
    J_area.add(-d_min*T*s.n*beta)
    for row in range(s.n):
        for t in range(T):
            J_area.add(-G[t][row,row]*s.weight[row])
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
    # set objective
    J_area=LinExpr()
    d_min=model.addVar(lb=0.0001)
    beta=10**2 # Weight of infinity norm
    model.update()
    for row in range(s.n):
        for column in range(s.n):
            if coin<0.1:
                if row<column:
                    model.addConstr(G[0][row,column]==0)
            elif coin>0.9:
                if row>column:
                    model.addConstr(G[0][row,column]==0)                
            if row==column:
                model.addConstr(G[0][row,column]>=d_min/s.weight[row])
    J_area.add(-d_min*T*s.n*beta)
    for row in range(s.n):
        for t in range(T):
            J_area.add(-G[t][row,row]*s.weight[row])
    # Terminal Constraint
    terminal_constraint(s,x,G,T,model,state_end)
    # Starting Point
    i_start=find_mode(s,x0)
    for i in s.modes:
        model.addConstr(z[0,i]==int(i==i_start))
    model.setParam('OutputFlag',False)
    if alpha_start==-1:
        x_delta={}
        for row in range(s.n):
            x_delta[row]=model.addVar(lb=-eps/s.weight[row],ub=eps/s.weight[row])
        model.update()
#        print("only area!")
        for row in range(s.n):
            model.addConstr(x[0][row,0]==x0[row,0]+x_delta[row])
        model.setObjective(J_area)
        model.optimize()
    else:          
        model.setObjective(J_area)
        model.optimize()
    if model.Status!=2 and model.Status!=11:
        flag=False
#        print("*"*20,"Flag is False and Status is",model.Status)
        print("*"*20,"False flag",model.Status)
        return (x,u,G,theta,z,flag)
    else:
        flag=True
#        print("*"*20,"Flag is True and Status is",model.Status)
        x_n=valuation(x)
        u_n=valuation(u)
        G_n=valuation(G)
        theta_n=valuation(theta)
        z_n=mode_sequence(s,z)
        if abs(np.linalg.det(G_n[0]))<10**-5:
            flag=False
        return (x_n,u_n,G_n,theta_n,z_n,flag)   