#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 18:12:19 2018

@author: sadra
"""

import numpy as np
from gurobipy import Model, GRB, LinExpr

def canonical_polytope_old(H,h,atol=10**-8):
    # Given H, h, provide canonical polytope by finding and removing redundant rows
    # Also scale H such that its largest absolute element is 1
    # REMOVE LATER
    n=H.shape[1]
    H_final=np.empty((0,n))
    H_max=np.amax(abs(H),axis=1)
    H_final=np.empty((0,n))
    h_final=np.empty((0,1))
    print("before removing redundancy",H,h)
    for ROW in range(H.shape[0]):
        if check_redundancy_row(H,h,ROW)==False:
            H_scale=np.asscalar(H_max[ROW])
            H_final=np.vstack((H_final,H[ROW,:]/H_scale))
            h_final=np.vstack((h_final,h[ROW,:]/H_scale))
        else:
            pass
    return (H_final,h_final)

def canonical_polytope(H,h,flag=None,atol=10**-8):
    """
    Given a polytope in form {H x <= h}, provide canonical polytope by finding and removing redundant rows
    Also scale H such that its largest absolute element is 1
    """
    if flag==None: # Construct flag for each row
        flag={} 
        for ROW in range(H.shape[0]):
            if check_redundancy_row(H,h,ROW)==True:
                flag[ROW]=False # It "may" be removed
            else:
                flag[ROW]=True # It must remain
        return canonical_polytope(H,h,flag) 
    elif [row for row in range(H.shape[0]) if flag[row]==False]==[]: # No row is needed to be removed
        return normalize(H,h)
    elif len([row for row in range(H.shape[0]) if flag[row]==False])==1: # Remove the only remaining false row
        a_false_row=[row for row in range(H.shape[0]) if flag[row]==False][0]
        removed_ROW=list(range(0,a_false_row))+list(range(a_false_row+1,H.shape[0]))
        H=H[removed_ROW,:]
        h=h[removed_ROW,:]
        return normalize(H,h)         
    else: # Let's remove one of the False rows
        a_false_row=[row for row in range(H.shape[0]) if flag[row]==False][0]
        rows_one_removed=list(range(0,a_false_row))+list(range(a_false_row+1,H.shape[0]))
        flag_new={}
        row=0
        for ROW in rows_one_removed:
            flag_new[row]=flag[ROW]
            row+=1
        H=H[rows_one_removed,:]
        h=h[rows_one_removed,:]
        for row in range(H.shape[0]):
            if flag_new[row]==False:
                if check_redundancy_row(H,h,row)==False:
                    flag_new[row]=True
        return canonical_polytope(H,h,flag_new,atol) 
                    

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
    H_max=np.amax(abs(H),axis=1)
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
        return canonical_polytope(A_new,b_new)

def project(T,d_translate,A,b,C=None,d=None,atol=10**-8):
    """
    Finds the H-representation of T{Ax<=b, Cx=d}+d_translate
    """
    (m,n)=T.shape # m: y, n: x, y=Tx
    if C==None and d==None:
        A=np.hstack((np.zeros((A.shape[0],m)),A))
        b=b
        C=np.hstack((-np.eye(m),T))
        d=-d_translate
        (A,b)=fourier_motzkin_eliminate_single(n+m-1,A,b,C,d,atol)
        for j in range(n-1):
            (A,b)=fourier_motzkin_eliminate_single(A.shape[1]-1,A,b,None,None,atol)
        return (A,b)