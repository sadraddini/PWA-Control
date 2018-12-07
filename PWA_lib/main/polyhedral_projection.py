#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 18:12:19 2018

@author: sadra
"""

import numpy as np
from gurobipy import Model, GRB, LinExpr

import sys
sys.path.append('..')

from main.auxilary_methods import rank_matrix

def canonical_polytope(H,h,flag=None,atol=10**-8):
    """
    Given a polytope in form {H x <= h}, provide canonical polytope by finding and removing redundant rows
    Also scale H such that its largest absolute element is 1
    """
    print(H.shape)
    row=0
    while row<H.shape[0]:
        if check_redundancy_row(H,h,row,atol)==True: # The row should be removed
            (H,h)=remove_row(H,h,row)
            row=row
        else:
            row+=1
    return normalize(H,h)
    
def remove_row(H,h,row):
    """
    Given {x| Hx <= h}, remove the row'th row of H and h
    """
    remaining_rows=list(range(0,row))+list(range(row+1,H.shape[0]))
    H_new=H[remaining_rows,:]
    h_new=h[remaining_rows,:]
    assert H_new.shape[0]==H.shape[0]-1
    return (H_new,h_new)
    
                  

def check_redundancy_row(H,h,ROW,atol=10**-8):
    model=Model("Row Redundancy Check")
    n=H.shape[1]
    x=np.empty((n,1),dtype='object')
    for row in range(n):
        x[row,0]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
    model.update()
    for row in [r for r in range(H.shape[0]) if r!=ROW]:
        Hx=LinExpr()
        for column in range(n):
            Hx.add(H[row,column]*x[column,0])
        model.addConstr(Hx<=h[row,0])
    J=LinExpr()
    for column in range(n):
        J.add(H[ROW,column]*x[column,0])
    model.setObjective(J, GRB.MAXIMIZE)
    model.setParam('OutputFlag',False)
    model.optimize()
    if model.Status==2:
        if J.getValue()>h[ROW,0]+atol:
            return False # It is NOT redundant
        else:
            return True # It is redudant
    else:
        return False
    
def normalize(H,h):
#    H_max=np.amax(abs(H),axis=1)
    H_max=np.linalg.norm(H,axis=1)
    H=np.dot(np.diag(1/H_max),H)
    h=np.dot(np.diag(1/H_max),h)
    return (H,h)
        
def fourier_motzkin_eliminate_single(var_index,A,b,C=None,d=None,atol=10**-8):
    if type(C)==type(np.array([1])):
        A=np.vstack((A,C,-C))
        b=np.vstack((b,d,-d))
        return fourier_motzkin_eliminate_single(var_index,A,b,None,None,atol)
    else:
        phi_positive=[i for i in range(A.shape[0]) if A[i,var_index]>=atol] # list of positive var entries
        phi_negative=[i for i in range(A.shape[0]) if A[i,var_index]<=-atol]  # list of negative var entries
        phi_core=[i for i in range(A.shape[0]) if abs(A[i,var_index])<atol]  # list of zero var entries
        s_smaller=np.diag(1/A[phi_positive,var_index]) # positive
        s_larger=np.diag(1/A[phi_negative,var_index]) # negative
        A_positive=np.dot(s_smaller,A[phi_positive,:]) # A of postives scaled by var entries
        b_positive=np.dot(s_smaller,b[phi_positive,:])
        A_negative=np.dot(s_larger,A[phi_negative,:])
        b_negative=np.dot(s_larger,b[phi_negative,:]) 
        """ We have A_positive x_other + x_r <= b_positive
        --> We have A_negative x_other + x_r >= b_negative
        --> We have b_postive - b_negative >= (A_neg - A _pos) * x_other (all combinations)
        """
        A_new=np.empty((0,A.shape[1]-1))
        b_new=np.empty((0,1))
        other=list(range(0,var_index))+list(range(var_index+1,A.shape[1]))
        for i in range(len(phi_positive)):
            for j in range(len(phi_negative)):
                alpha=(-A_negative[j,other]+A_positive[i,other]).reshape(1,len(other))
                beta=b_positive[i,:]-b_negative[j,:]
                A_new=np.vstack((A_new,alpha))
                b_new=np.vstack((b_new,beta))
        if phi_core!=[]:
            A_new=np.vstack((A_new,A[phi_core,:][:,other]))
            b_new=np.vstack((b_new,b[phi_core,:]))
        print("F-M elimination: ",A.shape," to ",A_new.shape)
        return canonical_polytope(A_new,b_new)

def project(T,d_translate,A,b,C=None,d=None,atol=10**-8):
    """
    Finds the H-representation of T{Ax<=b, Cx=d}+d_translate
    """
    print("projecting polytope!")
    (m,n)=T.shape # m: y, n: x, y=Tx
    print("rank of T is",rank_matrix(T))
    if m==n and rank_matrix(T)==n:
        Tinv=np.linalg.inv(T)
        A=np.dot(A,Tinv)
        b=b+np.dot(A,d_translate)
        return normalize(A,b)
    elif C==None and d==None:
        A=np.hstack((np.zeros((A.shape[0],m)),A))
        b=b
        C=np.hstack((-np.eye(m),T))
        d=-d_translate
        (A,b)=fourier_motzkin_eliminate_single(n+m-1,A,b,C,d,atol)
        print("Elimination 1")
        for j in range(n-1):
            print("Elimination %d"%(j+2))
            (A,b)=fourier_motzkin_eliminate_single(A.shape[1]-1,A,b,None,None,atol)
        print("end of projecting polytope")
        return (A,b)
    else:
        raise("I don't know how to handle this")
        
def check_copy_rows(H,h):
    """
    Remove dedundancies that come from sources 
    """
    for row in range(H.rows[0]):
        pass
    