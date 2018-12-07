# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 15:21:17 2018

@author: sadra
"""

import numpy as np
from gurobipy import Model,GRB,LinExpr

from pypolycontain.lib.containment_encodings import subset_LP_disjunctive

def point_trajectory(system,x0,list_of_goals,T,eps=[None]):
    """
    Description: point Trajectory Optimization
    Inputs:
        system: control system in the form of sPWA
        x_0: initial point
        T= trajectory length
        list_of_goals: reaching one of the goals is enough. Each goal is a zonotope
        eps= vector, box for how much freedom is given to deviate from x_0 in each direction
    """
    model=Model("Point Trajectory Optimization")
    x=model.addVars(T+1,range(system.n),lb=-GRB.INFINITY,ub=GRB.INFINITY)
    u=model.addVars(T,range(system.m),lb=-GRB.INFINITY,ub=GRB.INFINITY)
    index_x,index_u,index_delta=[],[],[]
    for t in range(T):
        for n in system.sum_PWA:
            for i in system.A[n].keys():
                index_delta.append((t,n,i))
                for j in range(system.n):
                    index_x.append((t,n,i,j))
                for j in range(system.m):
                    index_u.append((t,n,i,j))
    x_PWA=model.addVars(index_x,lb=-GRB.INFINITY,ub=GRB.INFINITY,name="x_pwa"+str(index_x))
    u_PWA=model.addVars(index_u,lb=-GRB.INFINITY,ub=GRB.INFINITY,name="u_pwa"+str(index_u))
    delta_PWA=model.addVars(index_delta,vtype=GRB.BINARY)
    model.update()
    # Initial Condition
    add_initial_condition(system,model,x,x0,eps)
    # Convexhull Dynamics
    model.addConstrs([t+1,j]==x_PWA.sum(t+1,"*","*",j) for t in range(T) for j in range(system.n))
    for t in range(T):
        for n in system.sum_PWA:
            for i in system.A[n].keys():
                model.addConstrs(x.prod(system.H[n,i],t,i,j,"*")<=system.h[n,i][j,0]*delta_PWA[t,n,i] for j in range(system.h[n,j].shape[0]))
    for t in range(T):
        for n in system.sum_PWA:
            for i in system.A[n].keys():
                model.addConstrs(x_PWA[t+1,n,i,j]==x_PWA.prod(system.A[n,i],t,i,j,"*")+u_PWA.prod(system.B[n,i],t,i,j,"*") for j in range(system.n))
    # Integer Variables
    model.addConstrs(delta_PWA.sum(t,n,"*")==1 for t in range(n) for n in system.sum_PWA)
    # Final Goal Constraints
    mu=model.addVars(list_of_goals,vtype=GRB.BINARY)
    _p=model.addVars(list_of_goals,range(system.n),lb=-1,ub=1)
#    _y=model.addVars(list_of_goals,range(system.n),lb=-GRB.INFINITY,ub=GRB.INFINITY)
    model.update()
    for j in range(system.n):
        L=LinExpr()
        for goal in list_of_goals:
            for k in range(goal.G.shape[1]):
                L.add(goal.G[j,k]*_p[goal,k])
            L.add(mu[goal]*goal.x[j,0])
        model.addConstr(L==x[T,j])
    # Cost Engineering
        
    # Optimize
    model.optimize()        
    return (x,u,delta_PWA,mu)


    
def polytopic_trajectory_given_modes(x0,list_of_polytopes,goal_polytope,eps=0.1):
    """
    Description: Polytopic Trajectory Optimization
    """
    raise NotImplemented






def polytopic_trajectory(system,x0,list_of_goal_polytopes,T,eps=0.1):
    """
    Description: Polytopic Trajectory Optimization
    """
    raise NotImplementedError    






def add_initial_condition(system,model,x,x0,eps):
    """
    eps=system
    """
    if any(eps)==None:
        model.addConstrs(x[0,i]==x0[i,0] for i in range(system.n))
    else:
        eps=model.addvars(range(system.n),lb=[-e for e in eps],ub=eps)
        model.addConstrs(x[0,i]==x0[i,0]+eps[i] for i in range(system.n))
