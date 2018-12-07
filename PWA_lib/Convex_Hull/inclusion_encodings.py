#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 16:03:39 2018

@author: sadra
"""

def Q_inside_conv_P_s(s,model,x,G,Q,list_of_polytopes):
    # WARNING: INCOMPLETE CODE
    """
    Given:
        s: PWA system that is considered
        Numbers:
            Q= {HQ x \le Hq} \subset R^q
            P_i={H_i x \le H_i} \subset R^n, i=1,...,N, 
        Symbols:
            model= Gurobi model
            T= matrix in R^{n*q}: model variable
            d= vector in R^n: model variable
    
    What the function does:
        Model with Constraints required to have TQ+d in Convexhull(Ps)
    
    Output:
        None
    """
    Lambda={}
    delta={}
    for polytope in list_of_polytopes:
        (n_H,n_HQ)=(polytope.H.shape[0],Q.H.shape[0])
        delta[polytope]=model.addVar(lb=0,ub=1)
        Lambda[polytope]=np.empty((n_H,n_Q))
        
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
        H=np.dot(s.Pi,polytope.G_eps_inv)
        h=np.ones((s.Pi.shape[0],1))+np.dot(H,polytope.x)
        factorize=np.amax(abs(H),1).reshape(s.Pi.shape[0],1)
        H=np.divide(H,factorize)
        h=np.divide(h,factorize)
        for row in range(n_h):
            for column in range(n_g):
                s_left=LinExpr()
                s_right=LinExpr()
                for k in range(n_p):
                    s_left.add(Lambda[polytope][row,k]*s.Pi[k,column])
                for k in range(n):
                    s_right.add(H[row,k]*G_pol[polytope][k,column])
                model.addConstr(s_left==s_right)
        # Lambda * 1 <= H*
        for row in range(n_h):
            s_left=LinExpr()
            s_right=LinExpr()
            for k in range(n_p):
                s_left.add(Lambda[polytope][row,k])
            for k in range(n):
                s_right.add(H[row,k]*x_pol[polytope][k,0])
            model.addConstr(s_left<=h[row,0]*z_pol[polytope]-s_right)
    model.addConstr(z_sum==1)
    for row in range(n):
        model.addConstr(x_sum[row,0]==x[row,0])
        for column in range(n_g):
            model.addConstr(G_sum[row,column]==G[row,column])
    return z_pol

def  TQd_inside_one_of_the_states(model,T,Q,d,z,list_of_states):
    """
    Inputs:
        model: Gurobi model
        T= Transformation Matrix
        
    """
    pass