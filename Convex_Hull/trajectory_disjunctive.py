#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 14:38:17 2018

@author: sadra
"""



# Primary imports
import numpy as np
from gurobipy import Model,GRB,LinExpr,QuadExpr
from random import choice as rchoice
from random import random
from time import time

### Internal imports
import sys
sys.path.append('..')

# Secondary imports
from main.auxilary_methods import find_mode,valuation,mode_sequence
from main.ana_system import state
from main.trajectory import subset_MILP

def polytopic_trajectory_to_set_of_polytopes(s,x0,T,list_of_polytopes,eps=0,method="bigM",timelimit=60):
    (model,x,u,G,theta,z)=s.library[T]
    n_vars=len(model.getVars())
    n_constraints=len(model.getConstrs())
    J_area=LinExpr()
    d_min=model.addVar(lb=0.0001,name="new var")
    beta=10**2 # Weight of infinity norm
    model.update()
    coin=random()
    for row in range(s.n):
        for column in range(s.n):
            if coin<0.1:
                if row<column:
                    model.addConstr(G[0][row,column]==0,name="constrain")
            elif coin>0.9:
                if row>column:
                    model.addConstr(G[0][row,column]==0)
            if row==column:
                model.addConstr(G[0][row,column]>=d_min/s.weight[row])
    J_area.add(-d_min*T*s.n*beta)
    for row in range(s.n):
        for t in range(T):
            J_area.add(-G[t][row,row]*s.weight[row])
    # Starting Point and initial condition
    i_start=find_mode(s,x0)
    for i in s.modes:
        model.addConstr(z[0,i]==int(i==i_start))
    x_delta={}
    for row in range(s.n):
        x_delta[row]=model.addVar(lb=-eps/s.weight[row],ub=eps/s.weight[row])
    model.update()
    for row in range(s.n):
        model.addConstr(x[0][row,0]==x0[row,0]+x_delta[row])
    model.setParam('OutputFlag',True)
    model.setParam('TimeLimit', timelimit)
    # Terminal Constraint
    if method=="bigM":
        z_pol=terminal_constraint_set_bigM(s,x[T],G[T],model,list_of_polytopes)
    elif method in ["convexhull","Convex_hull","convex_hull","chull"]:
        z_pol=terminal_constraint_convex_hull(s,x[T],G[T],model,list_of_polytopes)
    else:
        raise(method," was not recognized. Enter either 'big-M' or 'Convex_hull' as method")
    model.setObjective(J_area)
    model.optimize()
    if model.SolCount>0:
        flag=True
        print("+"*20,"Flag is True and Status is",model.Status)
        x_n=valuation(x)
        u_n=valuation(u)
        G_n=valuation(G)
        theta_n=valuation(theta)
        z_n=mode_sequence(s,z)
        if abs(np.linalg.det(G_n[0]))<10**-15:
            flag=False
        for polytope in list_of_polytopes:
            if abs(z_pol[polytope].X-1)<=0.1:
                state_end=polytope
#                print(polytope,polytope.x.T,polytope.G)
        final=(x_n,u_n,G_n,theta_n,z_n,flag,state_end)
    elif model.Status!=2 and model.Status!=11:
        flag=False
        print("-"*20,"False flag",model.Status)
        final=(x,u,G,theta,z,flag,s.goal)
    else:
        pass
    print("starting removal process")
    print("start=",n_vars,"variables and ",n_constraints," constraints")
    new_n_vars=len(model.getVars())
    new_n_constraints=len(model.getConstrs())
    print("new:",new_n_vars,"variables and ",new_n_constraints," constraints")    
    time_start=time()
    for i in range(new_n_vars-n_vars):
        model.remove(model.getVars()[-i-1])
    model.update()
    for i in range(new_n_constraints-n_constraints):
        model.remove(model.getConstrs()[-i-1]) 
    model.update()
    n_vars=len(model.getVars())
    n_constraints=len(model.getConstrs())
    print("final:",n_vars,"variables and ",n_constraints," constraints")    
    print("end of removal in ",time()-time_start," seconds")
    return final

def terminal_constraint_set_bigM(s,x,G,model,list_of_polytopes):
    """
    Facts: Here we use H-rep. Therefore, G_eps is used!
    Everything is used with bigM
    """
    Lambda={}
    (n,n_g)=G.shape
    (n_p,n)=s.Pi.shape
    (n_h,n)=(2*s.n,s.n)
    z_pol={}
    bigM=2
    for polytope in list_of_polytopes:
        Lambda[polytope]=np.empty((n_h,n_p),dtype='object')
        for row in range(n_h):
            for column in range(n_p):
                Lambda[polytope][row,column]=model.addVar(lb=0)
        z_pol[polytope]=model.addVar(vtype=GRB.BINARY)
    model.update()
    z_sum=LinExpr()
    for polytope in list_of_polytopes:
        z_sum.add(z_pol[polytope])
        for row in range(n_h):
            for column in range(n_g):
                s_left=LinExpr()
                s_right=LinExpr()
                for k in range(n_p):
                    s_left.add(Lambda[polytope][row,k]*s.Pi[k,column])
                for k in range(n):
                    s_right.add(polytope.polytope.H[row,k]*G[k,column])
                model.addConstr(s_left==s_right)
        # Lambda * 1 <= H*
        for row in range(n_h):
            s_left=LinExpr()
            s_right=LinExpr()
            for k in range(n_p):
                s_left.add(Lambda[polytope][row,k])
            for k in range(n):
                s_right.add(polytope.polytope.H[row,k]*x[k,0])
            model.addConstr(s_left<=polytope.polytope.h[row,0]-s_right+bigM-bigM*z_pol[polytope])
    model.addConstr(z_sum==1)
    return z_pol

def terminal_constraint_convex_hull(s,x,G,model,list_of_states):
    """
    Facts: Here we use H-rep. Therefore, G_eps is used!
    Everything is used with bigM
    """
    print("Using Convex Hull")
    Lambda={}
    z_pol={}
    x_pol={}
    G_pol={}
    for y in list_of_states:
        Lambda[y]=np.empty((y.polytope.H.shape[0],s.Pi.shape[0]),dtype='object')
        x_pol[y]=np.empty((s.n,1),dtype='object')
        G_pol[y]=np.empty((s.n,s.n),dtype='object')
        for row in range(y.polytope.H.shape[0]):
            for column in range(s.Pi.shape[0]):
                Lambda[y][row,column]=model.addVar(lb=0)
        z_pol[y]=model.addVar(vtype=GRB.BINARY)
        for row in range(s.n):
            x_pol[y][row,0]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
        for row in range(s.n):
            for column in range(s.n):
                G_pol[y][row,column]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)                
    model.update()
    z_sum=LinExpr()
    G_sum=np.empty((s.n,s.n),dtype='object')
    x_sum=np.empty((s.n,1),dtype='object')
    for row in range(s.n):
        x_sum[row,0]=LinExpr()
        for column in range(s.n):
            G_sum[row,column]=LinExpr()
    for y in list_of_states:
        z_sum.add(z_pol[y])
        for row in range(s.n):
            x_sum[row,0].add(x_pol[y][row,0])
            for column in range(s.n):
                G_sum[row,column].add(G_pol[y][row,column])        
        for row in range(y.polytope.H.shape[0]):
            for column in range(s.n):
                s_left=LinExpr()
                s_right=LinExpr()
                for k in range(s.Pi.shape[0]):
                    s_left.add(Lambda[y][row,k]*s.Pi[k,column])
                for k in range(s.n):
                    s_right.add(y.polytope.H[row,k]*G_pol[y][k,column])
                model.addConstr(s_left==s_right)
        for row in range(y.polytope.H.shape[0]):
            s_left=LinExpr()
            s_right=LinExpr()
            for k in range(s.Pi.shape[0]):
                s_left.add(Lambda[y][row,k])
            for k in range(s.n):
                s_right.add(y.polytope.H[row,k]*x_pol[y][k,0])
            model.addConstr(s_left<=y.polytope.h[row,0]*z_pol[y]-s_right)
    model.addConstr(z_sum==1)
    for row in range(s.n):
        model.addConstr(x_sum[row,0]==x[row,0])
        for column in range(s.n):
            model.addConstr(G_sum[row,column]==G[row,column])
    return z_pol


def terminal_constraint_convex_hull_old(s,x,G,model,list_of_polytopes):
    """
    Facts: Here we use H-rep. Therefore, G_eps is used!
    Everything is used with bigM
    """
    Lambda={}
    (n,n_g)=G.shape
    (n_p,n)=s.Pi.shape
    (n_h,n)=(2*s.n,s.n)
    z_pol={}
    x_pol={}
    G_pol={}
    G_bound=100
    for polytope in list_of_polytopes:
        Lambda[polytope]=np.empty((n_h,n_p),dtype='object')
        x_pol[polytope]=np.empty((n,1),dtype='object')
        G_pol[polytope]=np.empty((n,n_g),dtype='object')
        for row in range(n_h):
            for column in range(n_p):
                Lambda[polytope][row,column]=model.addVar(lb=0)
        z_pol[polytope]=model.addVar(vtype=GRB.BINARY)
        for row in range(n):
            x_pol[polytope][row,0]=model.addVar(lb=-G_bound,ub=G_bound)
        for row in range(n):
            for column in range(n_g):
                G_pol[polytope][row,column]=model.addVar(lb=-G_bound,ub=G_bound)                
    model.update()
    z_sum=LinExpr()
    G_sum=np.empty((n,n_g),dtype='object')
    x_sum=np.empty((n,1),dtype='object')
    for row in range(n):
        x_sum[row,0]=LinExpr()
        for column in range(n):
            G_sum[row,column]=LinExpr()
    for polytope in list_of_polytopes:
        z_sum.add(z_pol[polytope])
        for row in range(n):
            x_sum[row,0].add(x_pol[polytope][row,0])
            for column in range(n_g):
                G_sum[row,column].add(G_pol[polytope][row,column])        
        for row in range(n_h):
            for column in range(n_g):
                s_left=LinExpr()
                s_right=LinExpr()
                for k in range(n_p):
                    s_left.add(Lambda[polytope][row,k]*s.Pi[k,column])
                for k in range(n):
                    s_right.add(polytope.polytope.H[row,k]*G_pol[polytope][k,column])
                model.addConstr(s_left==s_right)
        for row in range(n_h):
            s_left=LinExpr()
            s_right=LinExpr()
            for k in range(n_p):
                s_left.add(Lambda[polytope][row,k])
            for k in range(n):
                s_right.add(polytope.polytope.H[row,k]*x_pol[polytope][k,0])
            model.addConstr(s_left<=polytope.polytope.h[row,0]*z_pol[polytope]-s_right)
    model.addConstr(z_sum==1)
    for row in range(n):
        model.addConstr(x_sum[row,0]==x[row,0])
        for column in range(n_g):
            model.addConstr(G_sum[row,column]==G[row,column])
    return z_pol
    