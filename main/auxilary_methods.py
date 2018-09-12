#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 14:55:26 2018

@author: sadra
"""
import numpy as np
from numpy.linalg import svd


def find_mode(s,x,eps=10**-9):
    """
    Description: Given system s and state x
    """
    for mode in s.modes:
        if np.amin(s.h[mode]-np.dot(s.H[mode],x))>=-eps:
            return mode
    raise(ValueError("Error! out of state space:",x.T))

def find_mode_control(s,x):
    """
    Description: Given system s and state x, find the closest mode
    """
    d_min=1000
    best_mode=None
    for mode in s.modes:
        d=-np.amin(s.h[mode]-np.dot(s.H[mode],x))
        if d<d_min:
            best_mode=mode
            d_min=d
    return best_mode


def valuation(x):
    """
    Description: given a set of Gurobi variables, output a similar object with values
    Input:
        x: dictionary or a vector, each val an numpy array, each entry a Gurobi variable
        output: x_n: dictionary with the same key as, each val an numpy array, each entry a float 
    """
    if type(x)==type(dict()):
        x_n={}
        for key,val in x.items():
            x_n[key]=np.ones(val.shape)
            (n_r,n_c)=val.shape
            for row in range(n_r):
                for column in range(n_c):
                    x_n[key][row,column]=x[key][row,column].X   
        return x_n
    elif type(x)==type(np.array([])):
        x_n=np.empty(x.shape)
        for row in range(x.shape[0]):
            for column in range(x.shape[1]):
                x_n[row,column]=x[row,column].X
        return x_n
    else:
        raise("x is neither a dictionary or a numpy array")

def mode_sequence(s,z):
    seq={}
    for key,val in z.items():
        t=key[0]
        if round(val.X)==1:
            if t in seq.keys():
                raise(ValueError("Two modes are assigned to time %d: %d and %d"%(t,seq[t],key[1])))
            seq[t]=key[1]
    return seq

def vertices_cube(T):
    """
    Description: 2**n * n array of vectors of vertices in unit cube in R^n
    """
    from itertools import product 
    v=list(product(*zip([-1]*T,[1]*T)))
    return np.array(v)

def sample(l,u):
    """
    Description: sample from box of l:lower corner, u: upper corner, with uniform distribution
    """
    return l+(u-l)*np.random.random_sample(l.shape)

def inside_X(s,x,eps=10**-9):
    """
    Description: x inside state-space: True, Otherwise: False
    """
    for mode in s.modes:
        if np.amin(s.h[mode]-np.dot(s.H[mode],x))>=-eps:
            return True
    return False

def PI(n):
    return np.vstack((np.eye(n),-np.eye(n)))   

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