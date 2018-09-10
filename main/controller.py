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

def control_vanilla(s,x):
    x_state=tree_locator_time(s,x)
    s.value_traj.append(x_state.time_to_go)
    if x_state==s.goal:
        print("*"*40,"YUHOOO! REACHED THE GOAL","*"*40)
        return "END"
    u_0=x_state.successor[1]
    u_f=np.dot(x_state.successor[2],np.dot(x_state.Ginv,x-x_state.x))
    print("u_0=",u_0.T)
    print("u_f=",u_f.T)
    u_control=u_0+u_f
    print("u_control=",u_control)
    s.control_traj.append(u_control)
    return u_control

def control_convex(s,x):
    x_state=tree_locator_time(s,x)
    state_target=x_state.successor[0]
    s.value_traj.append(x_state.time_to_go)
    if state_target==s.goal:
        print("*"*40,"YUHOOO! REACHED THE GOAL","*"*40)
        return "END"        
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
    print("predicted_state is",state_target.x+np.dot(state_target.G,valuation(p))-valuation(delta))
    return valuation(u)

def control_convex_old(s,x):
    # WARNING: To be removed later
    state_target=tree_locator_time(s,x).successor[0]
    if state_target==s.goal:
        print("*"*40,"YUHOOO! REACHED THE GOAL","*"*40)
        return "END"        
    model=Model("controller")
    p={}
    delta={}
    u={}
    u={}
    for row in range(s.n):
        p[row]=model.addVar(lb=-1,ub=1)
        delta[row]=model.addVar(lb=-100,ub=100)
    for row in range(s.m):
        u[row]=model.addVar(lb=-4,ub=4)
    model.update()
    i=find_mode(s,x)
    Ax=np.dot(s.A[i],x)
    for row in range(s.n):
        Bu=LinExpr()
        Gp=LinExpr()
        for k in range(s.m):
            Bu.add(s.B[i][row,k]*u[k])
        for k in range(s.n):
            Gp.add(state_target.G[row,k]*p[k])
        model.addConstr(Ax[row,0]+Bu+s.c[i][row,0]==state_target.x[row,0]+delta[row]+Gp)
    J=QuadExpr()
    for row in range(s.n):
        J.add(delta[row]*delta[row]/s.weight[row]+p[row]*p[row]*0.00)
    model.setObjective(J)
    model.optimize()
    return np.array(u[0].X).reshape(1,1)   