#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 15:37:31 2018

@author: sadra
"""

from gurobipy import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from random import choice as rchoice

from main.auxilary_methods import vertices_cube


class state:
    def __init__(self,x,G,mode,ID,t,character,epsilon=10**-8):
        self.name="s"+str(ID)+"-"+str(mode)+"-"+str(t)+"-"+str(character)
        self.x=x
        self.G=G
        self.Ginv=np.linalg.pinv(G)
        if np.linalg.det(G)<10**-6:
            (U,s,V)=np.linalg.svd(G)
            self.G_eps=np.dot(U,np.dot(np.diag(s+epsilon),V))
            self.G_eps_inv=np.linalg.pinv(self.G_eps)
        else: # Large volume:
            self.G_eps=self.G
            self.G_eps_inv=self.Ginv
        if G.shape[0]==G.shape[1]:
            self.volume=abs(np.linalg.det(G))
        self.volume_flag=self.volume>10**-8
        self.mode=mode
        self.ID=ID
        self.t=t
        self.character=character # 0 for goal, 1 for ordinary paralleltope, 2 for 
        v=vertices_cube(G.shape[1])
        self.vertices=(np.dot(G,v.T)).T
        self.backward_zero=-1 # 0 for not computed, -1 for computed but not useful, 1 some large regions leading it have been computed!
        self.cost_to_go=0
        self.time_to_go=0
        self.cost_to_child=0
        self.parent=[]
            
    def __repr__(self):
        return self.name
    
class system:
    def __init__(self,n=1,m=1,name="PWA System"):
        self.n=n   
        self.m=m
        self.modes=[1]
        self.A={}
        self.B={}
        self.c={}
        self.H={}
        self.h={}
        self.F={}
        self.f={}
        self.R={}
        self.r={}
        self.dt=0.1
        self.Pi=0 # Determine this as soon as possible!
        self.name=name
        # Box for each label!
        self.l={}
        self.u={}
        # Finite Abstraction
        self.X=[]
        self.streams=[]
        self.leafs=[]
        self.branches=[] # These are for tree!
        # List of no transitions
        self.blocked_transitions=[]
        # cost functions
        self.tree_iterations=0
        self.tree_size={}
        # ID
        self.ID=0
        
    def __repr__(self):
        return self.name+" with "+str(len(self.modes))+" modes"
    
def states_time_order(s):
    indices=list(np.argsort(np.array([x.time_to_go for x in s.X])))
    return [s.X[i] for i in indices]

def states_cost_order(s):
    indices=list(np.argsort(np.array([x.cost_to_go for x in s.X])))
    return [s.X[i] for i in indices]    

class tree:
    def __init__(self,goal):
        self.goal=goal
        self.nodes=[goal]
        self.edges=[]
        self.weight={}
        self.value_function={}
        self.successor={}
        
def cost_state_old(s,state_considered,L,Q,gamma):
    """
    Asscoiate each state and its child transition a cost
    Assumption: paralleltopes
    """
    if s==s.goal:
        return 0
    model=Model("trajectory of polytopes")
    p={}
    for row in range(s.n):
        p[row]=model.addVar(lb=-1,ub=1)
    model.update()
    GLG=np.dot(state_considered.G.T,np.dot(L,state_considered.G))
    theta=state_considered.successor[2]
    u=state_considered.successor[1]
    i=state_considered.mode
    theta_Q_theta=np.dot(theta.T,np.dot(Q,theta))
    J=QuadExpr()
    for row in range(s.n):
        for k in range(s.n):
            J.add(p[row]*p[k]*GLG[row,k]+p[row]*p[k]*theta_Q_theta[row,k])
    model.setParam('OutputFlag',False)
    model.setObjective(J)
    model.optimize()
    return model.ObjVal+np.asscalar(np.dot(state_considered.x.T,np.dot(L,state_considered.x))+np.dot(u.T,np.dot(Q,u))+gamma)  

def cost_state(s,state_considered,L,Q,gamma):
    """
    Asscoiate each state and its child transition a cost
    Assumption: paralleltopes
    """
    if s==s.goal:
        return 0
    GLG=np.dot(state_considered.G.T,np.dot(L,state_considered.G))
    theta=state_considered.successor[2]
    u=state_considered.successor[1]
    theta_Q_theta=np.dot(theta.T,np.dot(Q,theta))
    v=vertices_cube(s.n)
    J={}
    for index in range(v.shape[0]):
        p=v[index,:].reshape(s.n,1)
        J[index]=0
        for row in range(s.n):
            for k in range(s.n):
                J[index]+=np.asscalar(p[row]*p[k]*GLG[row,k]+p[row]*p[k]*theta_Q_theta[row,k])
    return max(J.values())+np.asscalar(np.dot(state_considered.x.T,np.dot(L,state_considered.x))+np.dot(u.T,np.dot(Q,u))+gamma)                