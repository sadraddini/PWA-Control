#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 11:19:03 2018

@author: sadra
"""

import sys
sys.path.append('..')


import numpy as np
from gurobipy import Model,GRB,LinExpr,QuadExpr

from main.polytope import polytope,TQ_to_polytope,if_subset
from main.auxilary_methods import PI,valuation


def post_polytope(s,state_X):
    """
    Inputs:
        state_X: its polytopic form with successor map
    Output:
        polytope of (AG+Btheta)P
    """
    i=state_X.mode
    theta=state_X.successor[2]
    u=state_X.successor[1]
    G=np.dot(s.A[i],state_X.G)+np.dot(s.B[i],theta)
    d=np.dot(s.A[i],state_X.x)+np.dot(s.B[i],u)+s.c[i]
    return TQ_to_polytope(G,d) 

def successor_check(s,state_X):
    p=post_polytope(s,state_X)
    return if_subset(p,state_X.successor[0].polytope)

def subset_polytope(p2,p1):
    # Return if p2 in p1
    n1=p1.H.shape[1]
    n2=p2.H.shape[1]
    if n1==n2:
        n=n1
    else:
        raise("ERROR: Two polytopes are in different dimensions: %d and %d"%(n1,n2))
    model=Model("if a polytope is a subset of another")
    Lambda_main=np.empty((p1.H.shape[0],p2.H.shape[0]),dtype='object')
    for row in range(p1.H.shape[0]):
        for column in range(p2.H.shape[0]):
            Lambda_main[row,column]=model.addVar(lb=0,ub=GRB.INFINITY)
    model.update()
    constraints_AB_eq_CD(model,Lambda_main,p2.H,p1.H,np.eye(n))
    constraints_Ab_smaller_c(model,Lambda_main,p2.h,p1.h)
    model.optimize()  
    return model.objVal

def subset_transformed_polytope(T,d,Q,p):
    # answer if TQ+d \subset p
    if T.shape[1]!=Q.H.shape[1]:
        raise("ERROR: dimension mistmatch")
    model=Model("if a transformed polytope is a subset of another")
    Lambda_main=np.empty((p.H.shape[0],Q.H.shape[0]),dtype='object')
    for row in range(p.H.shape[0]):
        for column in range(Q.H.shape[0]):
            Lambda_main[row,column]=model.addVar(lb=0,ub=GRB.INFINITY)
    model.update()
    constraints_AB_eq_CD(model,Lambda_main,Q.H,p.H,T)
    constraints_Ab_smaller_c(model,Lambda_main,Q.h,p.h-np.dot(p.H,d))
    model.optimize()  
    return model.objVal
 
def constraints_AB_eq_CD(model,A,B,C,D):
    for row in range(A.shape[0]):
        for column in range(B.shape[1]):
            lhs=LinExpr()
            rhs=LinExpr()
            for k in range(A.shape[1]):
                lhs.add(A[row,k]*B[k,column])
            for k in range(C.shape[1]):
                rhs.add(C[row,k]*D[k,column])
            model.addConstr(rhs==lhs)
            
def constraints_Ab_smaller_c(model,A,b,c):
    delta={}
    for row in range(A.shape[0]):
        lhs=LinExpr()
        delta[row]=model.addVar(lb=0,obj=1)
        model.update()
        for k in range(A.shape[1]):
            lhs.add(A[row,k]*b[k,0])
        model.addConstr(lhs<=c[row,0]+delta[row]) 
        
def verification_of_successors(s,list_of_states,option="mild"):
    if option=="heavy":
        Q=polytope(PI(s.n),np.ones((2*s.n,1)))
        for state_X in list_of_states:
            state_X.successor_flag={}
            state_X.successor_flag["Boolean-H"]=subset_polytope(post_polytope(s,state_X),state_X.successor[0].polytope)
            state_X.successor_flag["Epsilon"]=successor_check(s,state_X)
            i=state_X.mode
            T=np.dot(s.A[i],state_X.G)+np.dot(s.B[i],state_X.successor[2])
            d=np.dot(s.A[i],state_X.x)+np.dot(s.B[i],state_X.successor[1])+s.c[i]
            state_X.successor_flag["Boolean-Td"]=subset_transformed_polytope(T,d,Q,state_X.successor[0].polytope)
            state_X.successor_flag["Td in H"]=subset_transformed_polytope(state_X.G,state_X.x,Q,state_X.polytope)
    elif option=="mild":
        Q=polytope(PI(s.n),np.ones((2*s.n,1)))
        for state_X in list_of_states:
            state_X.successor_flag={}
            i=state_X.mode
            T=np.dot(s.A[i],state_X.G)+np.dot(s.B[i],state_X.successor[2])
            d=np.dot(s.A[i],state_X.x)+np.dot(s.B[i],state_X.successor[1])+s.c[i]
            state_X.successor_flag["Boolean-Td"]=subset_transformed_polytope(T,d,Q,state_X.successor[0].polytope)
            state_X.successor_flag["Td in H"]=subset_transformed_polytope(state_X.G,state_X.x,Q,state_X.polytope)       


def print_flags(s):
    verification_of_successors(s,[y for y in s.X if y!=s.goal])
    for y in [y for y in s.X if y!=s.goal]:
        print(y,"d:",y.polytope.dimensions,y.successor[0],"d:",y.successor[0].polytope.dimensions,"\t",y.successor_flag)
        
def re_verification(s,state_X,state_end):
    """
    This is a method to deal with numerical inaccuracies in high dimensioanl MILPs that are terminated early.
    Requires solving a linear program. If feasible, the computed control startegy is recomputed. 
    If infeasible, report back, raise a flag, and abort the tree extention
    Inputs: state_X, state_end
    Output: flag, u, theta
    """
    model=Model("verifying tree extention")
    Lambda_main=np.empty((state_end.polytope.H.shape[0],s.Pi.shape[0]),dtype='object')
    for row in range(state_end.polytope.H.shape[0]):
        for column in range(s.Pi.shape[0]):
            Lambda_main[row,column]=model.addVar(lb=0,ub=GRB.INFINITY)
    G_dynamic=np.empty((s.n,s.n),dtype='object')
    theta=np.empty((s.m,s.n),dtype='object')
    u=np.empty((s.m,1),dtype='object')
    Hx_center=np.empty((state_end.polytope.H.shape[0],1),dtype='object') #H(Ax+Bu)
    for row in range(s.n):
        for column in range(s.n):
            G_dynamic[row,column]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
    for row in range(s.m):
        u[row,0]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
        for column in range(s.n):
            theta[row,column]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
    for row in range(state_end.polytope.H.shape[0]):
        Hx_center[row,0]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
    model.update()
    i=state_X.mode
    for row in range(s.n):
        for column in range(s.n):
            AG=LinExpr()
            for k in range(s.n):
                AG.add(s.A[i][row,k]*state_X.G[row,column])
            for k in range(s.m):
                AG.add(s.B[i][row,k]*theta[k,column])
            model.addConstr(G_dynamic[row,column]==AG)
    for row in range(state_end.polytope.H.shape[0]):
        HAx=LinExpr()
        for k in range(s.n):
            for j in range(s.m):
                HAx.add(state_end.polytope.H[row,k]*s.A[i][k,j]*u[j,0]+np.dot(state_end.polytope.H,np.dot(s.A[i],state_end.x))[row,0])
        model.addConstr(HAx==Hx_center[row,0])
    constraints_AB_eq_CD(model,Lambda_main,s.Pi,state_end.polytope.H,G_dynamic)
    constraints_Ab_smaller_c(model,Lambda_main,np.ones((2*s.n,1)),state_end.polytope.h-Hx_center)
    subset(model,theta,s.Pi,s.F[i],s.f[i],u)
    model.optimize() 
    if model.Status==2:
        print("Tree extention verified succesuflly!")
        return (valuation(u),valuation(theta))
    else:
        print("\n-"*10,"Tree extention failed!","\n-"*10)
        return None
    
def subset(model,G,Pi,H,h,x):
    """
    Description: Add Farkas lemma constraints for subset inclusion of x+GP in {e|H.<h}
    Inputs: 
        model: Gurobi optimization model
        G: n * n_g generator matrix
        Pi:  n_pi*n matrix where {x|Pi x< 1} is the primitive polytope
        {x| Hx<h} is the set constraint
        x is the point
    Output:
        no direct output. Adds constraints to the model. 
        FUTURE: we may add lambda positive matrix as an output for debugging
    """
    (n,n_g)=G.shape
    (n_p,n)=Pi.shape
    (n_h,n)=H.shape
    Lambda=np.empty((n_h,n_p),dtype='object')
    for row in range(n_h):
        for column in range(n_p):
            Lambda[row,column]=model.addVar(lb=0)
    model.update()
    # Lambda * Pi = H * G
    for row in range(n_h):
        for column in range(n_g):
            s_left=LinExpr()
            s_right=LinExpr()
            for k in range(n_p):
                s_left.add(Lambda[row,k]*Pi[k,column])
            for k in range(n):
                s_right.add(H[row,k]*G[k,column])
            model.addConstr(s_left==s_right)
    # Lambda * 1 <= H*
    for row in range(n_h):
        s_left=LinExpr()
        s_right=LinExpr()
        for k in range(n_p):
            s_left.add(Lambda[row,k])
        for k in range(n):
            s_right.add(H[row,k]*x[k,0])
        model.addConstr(s_left<=h[row,0]-s_right)