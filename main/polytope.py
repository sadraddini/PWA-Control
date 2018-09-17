#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 11:13:56 2018

@author: sadra
"""

import sys
sys.path.append('..')


from numpy.linalg import svd
import numpy as np
from gurobipy import Model,GRB,LinExpr,QuadExpr
from random import random as rand

from main.auxilary_methods import PI,valuation
from main.polyhedral_projection import project


class polytope:
    """
    A polytope is an object. H-rep is used: {x | H x \le h}
    """
    def __init__(self,H,h,dimensions="full"):
        self.H=H
        self.h=h
        self.dimensions=dimensions
        
    def __repr__(self):
        return "polytope in R^%d"%self.H.shape[1]

def state_to_polytope(state_X,method,tol=10**-4):
    """
    Input:
        a state in the system that is a paralleltope
    Output:
        polytope for the system
    """                  
    return TQ_to_polytope(state_X.G,state_X.x,method,tol)
    
def TQ_to_polytope(T,d,method,tol=10**-4):
    """
    Input:
        s: system
        T: matrix in R^n * n_Q
        d: off-set in d
    Output:
        H-rep for TP+d, where P is the base polytope in s, P \in R^n
    """
    if method=="eps":
        (u,s,v)=svd(T)
        T=np.dot(u,np.dot(np.diag(s+tol*(s<tol)),v))
    elif method=="FM":
        pass
    else:
        raise("Unknown method for projection")
    n=T.shape[0]
    (H,h)=project(T,d,PI(n),np.ones((2*n,1)))
    rank=rank_matrix(T)          
    return polytope(H,h,rank)     
        


def rank_matrix(A, atol=1e-13, rtol=0):
    """Estimate the rank (i.e. the dimension of the nullspace) of a matrix.
    
    The algorithm used by this function is based on the singular value
    decomposition of `A`.

    Parameters
    ----------
    A : ndarray
        A should be at most 2-D.  A 1-D array with length n will be treated
        as a 2-D with shape (1, n)
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    Return value
    ------------
    r : int
        The estimated rank of the matrix.

    See also
    --------
    numpy.linalg.matrix_rank
        matrix_rank is basically the same as this function, but it does not
        provide the option of the absolute tolerance.
    """

    A = np.atleast_2d(A)
    s = svd(A, compute_uv=False)
    tol = max(atol, rtol * s[0])
    rank = int((s >= tol).sum())
    return rank


def nullspace(A, atol=1e-13, rtol=0):
    """Compute an approximate basis for the nullspace of A.

    The algorithm used by this function is based on the singular value
    decomposition of `A`.

    Parameters
    ----------
    A : ndarray
        A should be at most 2-D.  A 1-D array with length k will be treated
        as a 2-D with shape (1, k)
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    Return value
    ------------
    ns : ndarray
        If `A` is an array with shape (m, k), then `ns` will be an array
        with shape (k, n), where n is the estimated dimension of the
        nullspace of `A`.  The columns of `ns` are a basis for the
        nullspace; each element in numpy.dot(A, ns) will be approximately
        zero.
    """

    A = np.atleast_2d(A)
    u, s, vh = svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns

def sample_from_polytope(polytope):
    """
        A random point in H,h
    """
    model=Model("Polytope Sampling")
    n=polytope.H.shape[1]
    alpha=np.random.random((n,1))-0.5
    theta=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
    model.update()
    Hx_Anchor=np.dot(polytope.H,polytope.anchor)
    H_alpha=np.dot(polytope.H,alpha)
    for row in range(polytope.H.shape[0]):
        model.addConstr(H_alpha[row,0]*theta+Hx_Anchor[row,0]<=polytope.h[row])
    model.setObjective(theta,GRB.MINIMIZE)
    model.setParam('OutputFlag',False)
    model.optimize()
    theta_min=theta.X
    model.reset()
    model.setObjective(theta,GRB.MAXIMIZE)
    model.optimize()
    theta_max=theta.X
    c=rand()
    x_sample=(polytope.anchor+alpha*theta_min)*c+(polytope.anchor+alpha*theta_max)*(1-c)
    polytope.anchor=x_sample
    return x_sample

def anchor_point(polytope):
    """
        A point in H,h
    """
    model=Model("Polytope Sampling")
    n=polytope.H.shape[1]
    x=np.empty((n,1),dtype="object")
    rho=np.empty((polytope.H.shape[0],1),dtype="object")
    for row in range(n):
        x[row,0]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
    for row in range(polytope.H.shape[0]):
        rho[row,0]=model.addVar(lb=0,ub=GRB.INFINITY)
    model.update()
    J=QuadExpr(0)
    for row in range(polytope.H.shape[0]):
        a=LinExpr()
        for column in range(polytope.H.shape[1]):
            a.add(polytope.H[row,column]*x[column,0])
        model.addConstr(a+rho[row,0]==polytope.h[row])
        J.add(rho[row,0]*rho[row,0])
    model.setParam('OutputFlag',False)
    model.setObjective(J)
    model.optimize()
    return valuation(x)

def if_subset(p2,p1):
    # Answer if p2 is in p1
    return directed_distance_polytope_to_polytope(p1,p2,"uniform",None)


def hausdorff_distance_polytope_to_polytope(p1,p2,norm="uniform",ball_approximation=2):
    d12=directed_distance_polytope_to_polytope(p1,p2,norm,ball_approximation)
    d21=directed_distance_polytope_to_polytope(p2,p1,norm,ball_approximation)
    return max(d12,d21)
    
def directed_distance_polytope_to_polytope(p1,p2,norm,ball_approximation):
    # Distance between two polytopes
    # Inputs: polytope 1, polytope 2
    # Output: distance from polytope 2 to polytope 1: 
        # minimum epsilon such that polytope 2 \subset in polytope 1+ epsilon
    n1=p1.H.shape[1]
    n2=p2.H.shape[1]
    if n1==n2:
        n=n1
    else:
        raise("ERROR: Two polytopes are in different dimensions: %d and %d"%(n1,n2))
    model=Model("Distance between two polytopes")
    p_ball=ball_polytope(n,norm)
    Lambda_main=np.empty((p1.H.shape[0],p2.H.shape[0]),dtype='object')
    Lambda_ball=np.empty((p_ball.H.shape[0],p2.H.shape[0]),dtype='object')
    T_main=np.empty((n,n),dtype='object')
    T_ball=np.empty((n,n),dtype='object')
    d_main=np.empty((n,1),dtype='object')
    d_ball=np.empty((n,1),dtype='object')
    epsilon=model.addVar(lb=0,obj=1,ub=10)
    for row in range(p1.H.shape[0]):
        for column in range(p2.H.shape[0]):
            Lambda_main[row,column]=model.addVar(lb=0,ub=GRB.INFINITY)
    for row in range(p_ball.H.shape[0]):
        for column in range(p2.H.shape[0]):
            Lambda_ball[row,column]=model.addVar(lb=0,ub=GRB.INFINITY)  
    for row in range(n):
        for column in range(n):
            T_main[row,column]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
            T_ball[row,column]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY) 
        d_main[row,0]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY) 
        d_ball[row,0]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY) 
    model.update()
    constraints_AB_eq_CD(model,Lambda_main,p2.H,p1.H,T_main)
    constraints_AB_eq_CD(model,Lambda_ball,p2.H,p_ball.H,T_ball)
    constraints_AB_smaller_c_H_d(model,Lambda_main,p2.h,p1.h,1,p2.H,d_main)
    constraints_AB_smaller_c_H_d(model,Lambda_ball,p2.h,p_ball.h,epsilon,p_ball.H,d_ball)
    for row in range(n):
        for column in range(n):
            model.addConstr(T_main[row,column]+T_ball[row,column]==int(row==column))
        model.addConstr(d_ball[row,0]+d_main[row,0]==0)
    model.optimize()  
    return epsilon.X

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
            
def constraints_AB_smaller_c(model,A,b,c,epsilon):
    for row in range(A.shape[0]):
        lhs=LinExpr()
        for k in range(A.shape[1]):
            lhs.add(A[row,k]*b[k,0])
            model.addConstr(lhs<=c[row,0]*epsilon) 
            
def constraints_AB_smaller_c_H_d(model,A,b,c,epsilon,H,d):
    # A*b <= c*epsilon - H*d
    for row in range(A.shape[0]):
        lhs=LinExpr()
        Hd=LinExpr()
        for k in range(A.shape[1]):
            lhs.add(A[row,k]*b[k,0])
        for k in range(H.shape[1]):
            Hd.add(H[row,k]*d[k,0])
        model.addConstr(lhs<=c[row,0]*epsilon-Hd) 
            
def ball_polytope(n,norm):
    if norm=="uniform":
        return polytope(PI(n),np.ones((2*n,1)))
            