# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 15:21:17 2018

@author: sadra
"""

import numpy as np
import scipy as sp
from gurobipy import Model,GRB,LinExpr

import time as time

from pypolycontain.lib.containment_encodings import subset_LP,subset_zonotopes
from pypolycontain.lib.polytope import polytope
from pypolycontain.lib.zonotope import zonotope


def point_trajectory(system,x0,list_of_goals,T,eps=[None]):
    """
    Description: point Trajectory Optimization
    Inputs:
        system: control system in the form of sPWA
        x_0: initial point
        T= trajectory length
        list_of_goals: reaching one of the goals is enough. Each goal is a zonotope
        eps= vector, box for how much freedom is given to deviate from x_0 in each direction
    Method:
        Uses convexhull formulation
    """
    t_start=time.time()
    model=Model("Point Trajectory Optimization")
    x=model.addVars(range(T+1),range(system.n),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="x")
    u=model.addVars(range(T),range(system.m),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="u")
    ##
    x_PWA=model.addVars([(t,n,i,j) for t in range(T+1) for n in system.list_of_sum_indices \
                         for i in system.list_of_modes[n] for j in range(system.n)],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="x_pwa")
    u_PWA=model.addVars([(t,n,i,j) for t in range(T) for n in system.list_of_sum_indices \
                         for i in system.list_of_modes[n] for j in range(system.m)],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="u_pwa")
    delta_PWA=model.addVars([(t,n,i) for t in range(T) for n in system.list_of_sum_indices \
                             for i in system.list_of_modes[n]],vtype=GRB.BINARY,name="delta_pwa")
    model.update()
    # Initial Condition
    add_initial_condition(system,model,x,x0,eps)
    # Convexhull Dynamics
    model.addConstrs(x[t,j]==x_PWA.sum(t,n,"*",j) for t in range(T+1) for j in range(system.n)\
                     for n in system.list_of_sum_indices)
    model.addConstrs(u[t,j]==u_PWA.sum(t,n,"*",j) for t in range(T) for j in range(system.m)\
                     for n in system.list_of_sum_indices)
    for t in range(T):
        for n in system.list_of_sum_indices:
            for i in system.list_of_modes[n]:
                for j in range(system.C[n,i].H.shape[0]):
                    expr=LinExpr()
                    expr.add(LinExpr([(system.C[n,i].H[j,k],x_PWA[t,n,i,k]) for k in range(system.n)]))
                    expr.add(LinExpr([(system.C[n,i].H[j,k+system.n],u_PWA[t,n,i,k]) for k in range(system.m)])) 
                    model.addConstr(expr<=system.C[n,i].h[j,0]*delta_PWA[t,n,i])
    # Dynamics
    for t in range(T):
        for j in range(system.n):
            expr=LinExpr()
            for n in system.list_of_sum_indices:
                expr.add(LinExpr([(system.c[n,i][j,0],delta_PWA[t,n,i]) for i in system.list_of_modes[n]]))
                expr.add(LinExpr([(system.A[n,i][j,k],x_PWA[t,n,i,k]) for k in range(system.n) \
                    for i in system.list_of_modes[n]]))
                expr.add(LinExpr([(system.B[n,i][j,k],u_PWA[t,n,i,k]) for k in range(system.m) \
                    for i in system.list_of_modes[n]]))
            model.addConstr(x[t+1,j]==expr)
    # Integer Variables
    for t in range(T):
        for n in system.list_of_sum_indices:
            expr=LinExpr([(1.0,delta_PWA[t,n,i]) for i in system.list_of_modes[n]])
            model.addConstr(expr==1)
    # Final Goal Constraints
    mu=model.addVars(list_of_goals,vtype=GRB.BINARY)
    _p=model.addVars(list_of_goals,range(system.n),lb=-1,ub=1)
    model.update()
    for j in range(system.n):
        L=LinExpr()
        L.add(LinExpr([(goal.G[j,k],_p[goal,k]) for goal in list_of_goals for k in range(goal.G.shape[1])]))
        L.add(LinExpr([(goal.x[j,0],mu[goal]) for goal in list_of_goals]))
        model.addConstr(L==x[T,j])
    model.addConstr(mu.sum()==1)
    # Cost Engineering
    print "model built in",time.time()-t_start," seconds"
    # Optimize
    model.write("sadra.lp")
    model.optimize()
    u_num,x_num,delta_PWA_num,mu_num={},{},{},{}
    for t in range(T+1):
        x_num[t]=np.array([x[t,i].X for i in range(system.n)]).reshape(system.n,1)
    for t in range(T):
        u_num[t]=np.array([u[t,i].X for i in range(system.m)]).reshape(system.m,1)
    for t in range(T):
        for n in system.list_of_sum_indices:
            for i in system.list_of_modes[n]:
                delta_PWA_num[t,n,i]=delta_PWA[t,n,i].X
    for goal in list_of_goals:
        mu_num[goal]=mu[goal].X
#    for key,val in x_PWA.items():
#        print key,val.X
#    for key,val in u_PWA.items():
#        print key,val.X        
    return (x_num,u_num,delta_PWA_num,mu_num)


    
def polytopic_trajectory_given_modes(x0,list_of_cells,goal,eps=0,order=1):
    """
    Description: 
        Polytopic Trajectory Optimization with the ordered list of polytopes given
        This is a convex program as mode sequence is already given
        list_of_cells: each cell has the following attributes: A,B,c, and polytope(H,h)
    """
    
    model=Model("Fixed Mode Polytopic Trajectory")
    T=len(list_of_cells)
    n,m=list_of_cells[0].B.shape
    q=int(order*n)
    x=model.addVars(range(T+1),range(n),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="x")
    u=model.addVars(range(T),range(m),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="u")
    G=model.addVars(range(T+1),range(n),range(q),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="G")
    theta=model.addVars(range(T),range(m),range(q),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="theta")
    model.update()
    
    for j in range(n):
        model.addConstr(x[0,j]<=x0[j,0]+eps)
        model.addConstr(x[0,j]>=x0[j,0]-eps)

    for t in range(T):
        print "adding constraints of t",t
        cell=list_of_cells[t]
        A,B,c,p=cell.A,cell.B,cell.c,cell.p
        for j in range(n):
            expr_x=LinExpr([(A[j,k],x[t,k]) for k in range(n)])
            expr_u=LinExpr([(B[j,k],u[t,k]) for k in range(m)])
            model.addConstr(x[t+1,j]==expr_x+expr_u+c[j,0])
        for i in range(n):
            for j in range(q):
                expr_x=LinExpr([(A[i,k],G[t,k,j]) for k in range(n)])
                expr_u=LinExpr([(B[i,k],theta[t,k,j]) for k in range(m)])
                model.addConstr(G[t+1,i,j]==expr_x+expr_u)
        x_t=np.array([x[t,j] for j in range(n)]).reshape(n,1)
        u_t=np.array([u[t,j] for j in range(m)]).reshape(m,1)
        G_t=np.array([G[t,i,j] for i in range(n) for j in range(q)]).reshape(n,q)
        theta_t=np.array([theta[t,i,j] for i in range(m) for j in range(q)]).reshape(m,q)
        GT=sp.linalg.block_diag(G_t,theta_t)
        xu=np.vstack((x_t,u_t))
        subset_LP(model,xu,GT,Ball(2*q),p)


    x_T=np.array([x[T,j] for j in range(n)]).reshape(n,1)
    G_T=np.array([G[T,i,j] for i in range(n) for j in range(q)]).reshape(n,q)
    z=zonotope(x_T,G_T)
    subset_zonotopes(model,z,goal)
    J=LinExpr([(1/(t+1.0),G[t,i,i]) for t in range(T+1) for i in range(n)])
    model.setObjective(J)
    model.write("sadra.lp")
    model.setParam('TimeLimit', 150)
    model.optimize()
    x_num,G_num={},{}
    for t in range(T+1):
        x_num[t]=np.array([[x[t,j].X] for j in range(n)]).reshape(n,1)
        G_num[t]=np.array([[G[t,i,j].X] for i in range(n) for j in range(q)]).reshape(n,q)
    return (x_num,G_num)






def polytopic_trajectory(system,x0,list_of_goal_polytopes,T,eps=0.1):
    """
    Description: Polytopic Trajectory Optimization
    """
    raise NotImplementedError    



def Ball(n):
    H=np.vstack((np.eye(n),-np.eye(n)))
    h=np.ones((2*n,1))
    return polytope(H,h)


def add_initial_condition(system,model,x,x0,eps):
    """
    eps=system
    """
    if eps==[None]:
        print "no epsilon is given"
        model.addConstrs(x[0,i]==x0[i,0] for i in range(system.n))
    else:
        eps=model.addVars(range(system.n),lb=[-e for e in eps],ub=eps)
        model.update()
        model.addConstrs(x[0,i]==x0[i,0]+eps[i] for i in range(system.n))


def tupledict(A):
    """
    Converts matrix to dict
    """
    if len(A.shape)==1:
        return {i:A[i] for i in range(A.shape[0])}
    elif len(A.shape)==2:
        return {(i,j):A[i,j] for i in range(A.shape[0]) for j in range(A.shape[1])}
    else:
        raise NotImplementedError
        
def valuation_t(x):
    """
    Description: given a set of Gurobi variables, output a similar object with values
    Input:
        x: dictionary or a vector, each val an numpy array, each entry a Gurobi variable
        output: x_n: dictionary with the same key as, each val an numpy array, each entry a float 
    """
    if type(x)==type(dict()):
        for key,val in x.items():
            x_n[key]=ones(val.shape)
            (n_r,n_c)=val.shape
            for row in range(n_r):
                for column in range(n_c):
                    x_n[key][row,column]=x[key][row,column].X   
        return x_n
    else:
        raise("x is neither a dictionary or a numpy array")