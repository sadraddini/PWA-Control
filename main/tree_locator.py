#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 10:38:54 2018

@author: sadra
"""

import numpy as np

from main.auxilary_methods import find_mode,find_mode_control


def array_tree(s):
    N=len(s.X)
    s.x_vector=np.empty((N,s.n,1))
    s.G_eps_vector=np.empty((N,s.n,s.n))
    s.G_eps_inv_vector=np.empty((N,s.n,s.n))
    for index in range(N):
        s.x_vector[index,:,:]=s.X[index].x
        s.G_eps_vector[index,:,:]=s.X[index].G_eps
        s.G_eps_inv_vector[index,:,:]=s.X[index].G_eps_inv

def d_tree(s,x,bigM_mode=100):
    """
    Note: array_tree should be constructed before.
    """
    i=find_mode_control(s,x)
    mode_distance=np.array([y.mode!=i for y in s.X]).reshape(len(s.X),1)*bigM_mode
    p_eps=np.matmul(s.G_eps_inv_vector,x-s.x_vector)
    p_eps_star=np.maximum(np.minimum(p_eps,1),-1)
    d=np.amax(abs(np.matmul(s.G_eps_vector,p_eps_star-p_eps)),axis=1)
    return d+mode_distance

def sorted_distance_states(s,x):
    N=len(s.X)
    d=d_tree(s,x).reshape(N)
    indices=list(np.argsort(d))
    return [s.X[i] for i in indices]

def zero_distance_states(s,x):
    return [s.X[index] for index in range(len(s.X)) if d_tree(s,x)[index][0]<10**-6]

def tree_locator_time(s,x):
    STATES=zero_distance_states(s,x)
    if STATES==[]:
        print("---- out of tree")
        return sorted_distance_states(s,x)[0]
    else:
        print("++++ inside the tree")
        T=1000
        for X in STATES:
            if X.time_to_go<T:
                T=X.time_to_go
                X_best=X
        return X_best
            

def inside_tree(s,x):
    if zero_distance_states(s,x)==[]:
        return False
    else:
        return True
    
def all_vertices_out_of_tree(s,x,G):
    if inside_tree(s,x)==True:
        return False
    for i in range(2**s.n):
        if (s,x+np.dot(G,s.vertices[i,:].T))==True:
            return False
    return True

def inside_tree_iterations(s,x,it):
    zero_distance_states=[s.X[index] for index in range(s.tree_size[it]) if d_tree(s,x)[index][0]<10**-6]
    if zero_distance_states==[]:
        return False
    else:
        return True