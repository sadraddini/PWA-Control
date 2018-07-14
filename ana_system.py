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

from auxilary_methods import vertices_cube


class state:
    def __init__(self,x,G,mode,ID,t,character):
        self.name="s"+str(ID)+"-"+str(mode)+"-"+str(t)+"-"+str(character)
        self.x=x
        self.G=G
        self.Ginv=np.linalg.pinv(G)
        if G.shape[0]==G.shape[1]:
            self.volume=abs(np.linalg.det(G))
            self.volume_flag=abs(np.linalg.det(G))>10**-9
        self.mode=mode
        self.ID=ID
        self.t=t
        self.character=character # 1 for in, 2 for out, 0 for string, -1 for self-loop, 3 for weaving, 4 for backward Funnel, 6: free end funnel, 7: state_end funnel
        v=vertices_cube(G.shape[1])
        self.vertices=(np.dot(G,v.T)).T
        self.backward_zero=-1 # 0 for not computed, -1 for computed but not useful, 1 some large regions leading it have been computed!
        self.cost_to_go=0
        self.time_to_go=0
            
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
        
    def __repr__(self):
        return self.name+" with "+str(len(self.modes))+" modes"

class tree:
    def __init__(self,goal):
        self.goal=goal
        self.nodes=[goal]
        self.edges=[]
        self.weight={}
        self.value_function={}
        self.successor={}