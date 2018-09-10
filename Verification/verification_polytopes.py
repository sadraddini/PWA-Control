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
from main.auxilary_methods import PI


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
        
def verification_of_successors(s,list_of_states):
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

def print_flags(s):
    verification_of_successors(s,[y for y in s.X if y!=s.goal])
    for y in [y for y in s.X if y!=s.goal]:
        print(y,"d:",y.polytope.dimensions,y.successor[0],"d:",y.successor[0].polytope.dimensions,"\t",y.successor_flag)