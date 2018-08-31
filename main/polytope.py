#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 11:13:56 2018

@author: sadra
"""
from numpy.linalg import svd
import numpy as np
from gurobipy import Model,GRB,LinExpr

from main.auxilary_methods import PI


class polytope:
    """
    A polytope is an object. H-rep is used: {x | H x \le h}
    """
    def __init__(self,H,h):
        self.H=H
        self.h=h
        
    def __repr__(self):
        return "polytope in R^%d"%self.H.shape[1]

def state_to_polytope(state_X,tol=10**-8):
    """
    Input:
        a state in the system that is a paralleltope
    Output:
        polytope for the system
    """                  
    return TQ_to_polytope(state_X.G,state_X.x,tol)
    
def TQ_to_polytope(T,d,tol=10**-8):
    """
    Input:
        s: system
        T: matrix in R^n * n_Q
        d: off-set in d
    Output:
        H-rep for TP+d, where P is the base polytope in s, P \in R^n
    """
    n=d.shape[0]
    u,sigma,v = svd(T)
    rank = int((sigma >= tol).sum())
    sigma_reduced=sigma[0:rank]
    sigma_plus=np.zeros(n)
    sigma_plus[0:rank]=1/sigma_reduced
    T_pinv=np.dot(v.T,np.dot(np.diag(sigma_plus),u.T))
    Pi=PI(n)
    H=np.dot(Pi,T_pinv)
    h=np.ones((Pi.shape[0],1))+np.dot(H,d)
    H_other=u[:,rank:n].T
    h_other=np.dot(u[:,rank:n].T,d)
    H=np.vstack((H,H_other,-H_other))
    h=np.vstack((h,h_other,-h_other))
    (H_final,h_final)=canonical_polytope(H,h)              
    return polytope(H_final,h_final)

def canonical_polytope(H,h):
    # Given H, h, provide canonical polytope by finding and removing redundant rows
    # Also scale H such that its largest absolute element is 1
    n=H.shape[1]
    H_final=np.empty((0,n))
    H_max=np.amax(abs(H),axis=1)
    H_final=np.empty((0,n))
    h_final=np.empty((0,1))
    for ROW in range(H.shape[0]):
        if check_redundancy_row(H,h,ROW)==False:
            H_scale=np.asscalar(H_max[ROW])
            H_final=np.vstack((H_final,H[ROW,:]/H_scale))
            h_final=np.vstack((h_final,h[ROW,:]/H_scale))
    if H_final.shape[0]!=2*n:
        raise("i is", H_final.shape,"polytope H-rep redundancy exists while it should not. check the code") 
    return (H_final,h_final)
                    

def check_redundancy_row(H,h,ROW):
    model=Model("Row Redundancy Check")
    n=H.shape[1]
    x=np.empty((n,1),dtype='object')
    for row in range(n):
        x[row,0]=model.addVar(lb=-100,ub=100)
    model.update()
    for row in [r for r in range(H.shape[0]) if r!=ROW]:
        Hx=LinExpr()
        for column in range(n):
            Hx.add(H[row,column]*x[column,0])
        model.addConstr(Hx<=h[row,0])
    J=LinExpr()
    for column in range(n):
        J.add(x[column,0]*H[ROW,column])
    model.setObjective(J, GRB.MAXIMIZE)
    model.setParam('OutputFlag',False)
    model.optimize()
    if J.getValue()>h[ROW,0]:
        return False # It is NOT redundant
    else:
        return True # It is redudncat
    
        


def rank(A, atol=1e-13, rtol=0):
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
        