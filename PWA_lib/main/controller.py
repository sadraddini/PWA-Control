#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 10:51:52 2018

@author: sadra
"""
import numpy as np
from gurobipy import Model,GRB,LinExpr,QuadExpr

from main.tree_locator import tree_locator_time
from main.auxilary_methods import find_mode,valuation,find_mode_control
from main.trajectory import subset_MILP

def control_vanilla(s,x):
    x_state=tree_locator_time(s,x)
    s.value_traj.append(x_state.time_to_go)
    if x_state==s.goal:
        print("*"*40,"YUHOOO! REACHED THE GOAL","*"*40)
        return [[None]*s.m] 
    u_0=x_state.successor[1]
    u_f=np.dot(x_state.successor[2],np.dot(x_state.G_eps_inv,x-x_state.x))
    print("u_0=",u_0.T)
    print("u_f=",u_f.T)
    u_control=u_0+u_f
    u_control=saturate_control(s,x,u_control)
    print("u_control=",u_control)
    s.control_traj.append(u_control)
    return u_control

def saturate_control(s,x,u_control):
    i=find_mode_control(s,x)
    scale=np.divide(np.dot(s.F[i],u_control),s.f[i])
    r=np.amax(scale)
    if r>1:
        return u_control/r
    else:
        return u_control
    
def control_convex(s,x):
    x_state=tree_locator_time(s,x)
    state_target=x_state.successor[0]
    s.value_traj.append(x_state.time_to_go)
    if state_target==s.goal:
        print("*"*40,"YUHOOO! REACHED THE GOAL","*"*40)
        return [[None]*s.m]        
    model=Model("controller")
    p=np.empty((s.n,1),dtype='object')
    delta=np.empty((s.n,1),dtype='object')
    u=np.empty((s.m,1),dtype='object')
    for row in range(s.n):
        p[row,0]=model.addVar(lb=-1,ub=1)
        delta[row,0]=model.addVar(lb=-100,ub=100)
    for row in range(s.m):
        u[row,0]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
    model.update()
    model.setParam('OutputFlag',False)
    i=find_mode_control(s,x)
    Ax=np.dot(s.A[i],x)
    for row in range(s.n):
        Bu=LinExpr()
        Gp=LinExpr()
        for k in range(s.m):
            Bu.add(s.B[i][row,k]*u[k,0])
        for k in range(s.n):
            Gp.add(state_target.G[row,k]*p[k,0])
        model.addConstr(Ax[row,0]+Bu+s.c[i][row,0]+delta[row,0]==state_target.x[row,0]+Gp)
    J=QuadExpr()
    #print (x_state.successor_flag)
    for row in range(s.n):
        J.add(delta[row,0]*delta[row,0]*s.weight[row]*s.weight[row])
    for row in range(s.F[i].shape[0]):
        Fu=LinExpr()
        for column in range(s.m):
            Fu.add(s.F[i][row,column]*u[column,0])
        model.addConstr(Fu<=s.f[i][row,0])
    model.setObjective(J)
    model.optimize()
    print(valuation(delta))
    print(valuation(p))
#    print("predicted_state is",state_target.x+np.dot(state_target.G,valuation(p))-valuation(delta))
    return valuation(u) 
  

def control_MPC(s,x0,H):
    x_state=tree_locator_time(s,x0)
    s.value_traj.append(x_state.time_to_go)
    state_target=x_state
    for j in range(H):
        state_target=state_target.successor[0]
    if state_target==s.goal:
        print("*"*40,"YUHOOO! REACHED THE GOAL","*"*40)
        return [[None]*s.m]        
    model=Model("controller_MPC")
    p=np.empty((s.n,1),dtype='object')
    delta=np.empty((s.n,1),dtype='object')
    delta_0=np.empty((s.n,1),dtype='object')
    x={}
    u={}
    z={}
    for t in range(H+1):
        x[t]=np.empty((s.n,1),dtype='object')
        u[t]=np.empty((s.m,1),dtype='object')
    for row in range(s.n):
        p[row,0]=model.addVar(lb=-1,ub=1)
        delta[row,0]=model.addVar(lb=-100,ub=100)
        delta_0[row,0]=model.addVar(lb=-100,ub=100)
    for t in range(H+1):
        for row in range(s.m):
            u[t][row,0]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
        for row in range(s.n):
            x[t][row,0]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
    for t in range(H):
        for i in s.modes:
            z[t,i]=model.addVar(vtype=GRB.BINARY)
    model.update()
    model.setParam('OutputFlag',True)
    # Trajectory Constraints:
    bigM=100
    for t in range(H):
        for i in s.modes:
            for row in range(s.n):
                Ax=LinExpr()
                for k in range(s.n):
                    Ax.add(s.A[i][row,k]*x[t][k,0])
                for k in range(s.m):
                    Ax.add(s.B[i][row,k]*u[t][k,0])
                model.addConstr(x[t+1][row,0]<=Ax+s.c[i][row]+bigM-bigM*z[t,i])
                model.addConstr(x[t+1][row,0]>=Ax+s.c[i][row]-bigM+bigM*z[t,i])
    # Constraints of modes:
    for t in range(H):
        sum_z=LinExpr()
        for i in s.modes:
            sum_z.add(z[t,i])
        model.addConstr(sum_z==1)
    # Constraints of mode subsets
    for t in range(H):
        for i in s.modes:
            n_h=s.H[i].shape[0]
            for row in range(n_h):
                s_left=LinExpr()
                for k in range(s.n):
                    s_left.add(s.H[i][row,k]*x[t][k,0])
                model.addConstr(s_left<=s.h[i][row,0]+bigM-bigM*z[t,i]) 
            n_h=s.F[i].shape[0]
            for row in range(n_h):
                s_left=LinExpr()
                for k in range(s.m):
                    s_left.add(s.F[i][row,k]*u[t][k,0])
                model.addConstr(s_left<=s.f[i][row,0]+bigM-bigM*z[t,i]) 
    # set objective
    J=QuadExpr()
    for row in range(s.n):
        model.addConstr(x[0][row,0]==x0[row,0]+delta_0[row,0])
    #print (x_state.successor_flag)
    for row in range(s.n):
        J.add(delta[row,0]*delta[row,0]*s.weight[row]*s.weight[row])
        J.add(1000*delta_0[row,0]*delta_0[row,0]*s.weight[row]*s.weight[row])
    for row in range(s.n):
        Gp=LinExpr()
        for k in range(s.n):
            Gp.add(state_target.G[row,k]*p[k,0])
        model.addConstr(x[H][row,0]+delta[row,0]==state_target.x[row,0]+Gp)  
    model.setObjective(J)
    model.optimize()
    u_n=valuation(u)
    print(valuation(delta).T)
    print(valuation(delta_0).T)
    return u_n[0]