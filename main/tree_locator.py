#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 10:38:54 2018

@author: sadra
"""

import numpy as np

from main.auxilary_methods import find_mode


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
    i=find_mode(s,x)
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
    return [s.X[index] for index in range(len(s.X)) if d_tree(s,x)[index][0]<10**-8]
    
def tree_locator(s,x):
    pass

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
    
"""
OLD MATERIAL: Going to be removed
"""    

def array(s,target):
    """
    Description: 
    """
    I_left=np.empty((len(target),s.Pi.shape[0],s.n))
    I_right=np.empty((len(target),s.Pi.shape[0],1))
    volume=np.empty((len(target),1))
    x_vector=np.empty((len(target),s.n,1))
    G_vector=np.empty((len(target),s.n,s.n))
    Ginv_vector=np.empty((len(target),s.n,s.n))
    G_VolDiff_vector=np.empty((len(target),s.n,s.n))
    for index in range(len(target)):
        I_left[index,:,:]=np.dot(s.Pi,target[index].Ginv)
        I_right[index,:,:]=np.ones((s.Pi.shape[0],1))+np.dot(I_left[index,:,:],target[index].x)
        volume[index,:]=target[index].volume_flag
        x_vector[index,:,:]=target[index].x 
        G_vector[index,:,:]=target[index].G
        Ginv_vector[index,:,:]=target[index].Ginv
        G_VolDiff_vector[index,:,:]=np.dot(target[index].G,target[index].Ginv)-np.eye(s.n)
    s.I_left=I_left
    s.I_right=I_right
    s.volume=volume
    s.x_vector=x_vector
    s.G_vector=G_vector
    s.Ginv_vector=Ginv_vector
    s.G_VolDiff_vector=G_VolDiff_vector
    
def state_closest_distance(s,x,target):
    d=np.amax(abs(s.x_vector-x),1) # infinity norm distance between centers
    #d=(s.volume)*np.amax(d)+d # 0 volume: d, 1 volume: d+d_max
    i=find_mode(s,x)
    mode_distance=np.array([y.mode!=i for y in target]).reshape(len(target),1) # 0 if same mode, 1 otherwise
    bigM=10**8
    d=d+mode_distance*bigM
    return (target[np.argmin(d)],np.amin(d))

def state_closest_volume(s,x,target):
    i=find_mode(s,x)
    mode_distance=np.array([y.mode!=i for y in target]).reshape(len(target),1) # 0 if same mode, 1 otherwise
    p=np.matmul(s.Ginv_vector,x-s.x_vector)
    p_box=np.maximum(np.minimum(p,1),-1)
    q=abs(np.matmul(s.G_vector,p_box)+s.x_vector-x)
    q=np.amax(q,1)
    bigM=10**8
    q=q+mode_distance*bigM
    return (target[np.argmin(q)],np.amin(q))

def state_in_graph_set_volumes(s,x,target):
    q=s.I_right-np.matmul(s.I_left, x) #Positive: in, negative: out
    q=0.5*(abs(q)-q) # Rectifier, zero if positive, non-zero if out
    q=np.amax(q,axis=1) # Produces the max
    thinness=np.matmul(s.G_VolDiff_vector,x-s.x_vector)
    thinness=np.amax(abs(thinness),axis=1)<=10**-8 # Zero= in the range-good, 1: bad
    q=(1-thinness)*np.max(q)+q # zero: good, does not effect
#    q=(1-s.volume)*np.max(q)+q # 0 volume: q+q_max, 1 volume: q
    return [target[index] for index in range(len(q)) if q[index]<=10**-8]        

def state_in_graph_set_distances(s,x,target):
    d=np.amax(abs(s.x_vector-x),1) # infinity norm distance between centers
    return [target[index] for index in range(len(d)) if d[index]<=10**-8]


def tree_locator_old(s,x):
    """
        WARNING: to be removed in the next commit
    """
    target=s.X
    exact=state_in_graph_set_volumes(s,x,target)+state_in_graph_set_distances(s,x,target)
    if exact!=[]:
        print("+"*10,"exact machine state, %d states available"%len(exact))
        return exact[np.argmin(np.array([state.time_to_go for state in exact]))]
    else:
        (state_d,d)=state_closest_distance(s,x,target)
        (state_v,v)=state_closest_volume(s,x,target)
        if d>v:
            print("-"*10,"Warning: out of finite-state system: volume based:",v)
            return state_v
        else:
            print("-"*10,"Warning: out of finite-state system: distance based:",d)
            return state_d 
        
def inside_tree_old(s,x,eps):
    array(s,s.X)
    exact=state_in_graph_set_volumes(s,x,s.X)+state_in_graph_set_distances(s,x,s.X)
    if exact!=[]:
        return True
    else:
        (state_d,d)=state_closest_distance(s,x,s.X)
        (state_v,v)=state_closest_volume(s,x,s.X)
        return min(v,d)<=eps
    return state_d

def tree_K_nearest(s,x,K):
    array(s,s.X)
    d=np.amax(abs(s.x_vector-x),1) # infinity norm distance between centers
    #d=(s.volume)*np.amax(d)+d # 0 volume: d, 1 volume: d+d_max
    i=find_mode(s,x)
    mode_distance=np.array([y.mode!=i for y in target]).reshape(len(target),1) # 0 if same mode, 1 otherwise
    bigM=10**8
    d=d+mode_distance*bigM
    return (target[np.argmin(d)],np.amin(d))

    i=find_mode(s,x)
    mode_distance=np.array([y.mode!=i for y in target]).reshape(len(target),1) # 0 if same mode, 1 otherwise
    p=np.matmul(s.Ginv_vector,x-s.x_vector)
    p_box=np.maximum(np.minimum(p,1),-1)
    q=abs(np.matmul(s.G_vector,p_box)+s.x_vector-x)
    q=np.amax(q,1)
    bigM=10**8
    q=q+mode_distance*bigM
    
    
    exact=state_in_graph_set_volumes(s,x,s.X)+state_in_graph_set_distances(s,x,s.X)
    (state_d,d)=state_closest_distance(s,x,s.X)
    (state_v,v)=state_closest_volume(s,x,s.X)
    return state_d