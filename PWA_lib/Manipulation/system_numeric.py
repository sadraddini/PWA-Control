#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 16:06:32 2019

@author: sadra
"""
import numpy as np
import scipy as sp
from gurobipy import Model,GRB,LinExpr,QuadExpr,tuplelist

class system_hard_contact_PWA_numeric:
    def __init__(self,symbolic_system):
        self.A={}
        self.B_u={}
        self.B_lambda={}
        self.c={}
        
        self.Eta=[]
        self.E={}
        self.e={}
        
        self.symbolic_system=symbolic_system
        self.list_of_contact_points=symbolic_system.list_of_contact_points
        self.n,self.m_u,self.m_lambda=symbolic_system.n,symbolic_system.m_u,symbolic_system.m_lambda

       
    def add_environment(self,Eta,epsilon_max=[0],epsilon_min=[0]):
        self.Eta.append(Eta)
        self.E[Eta],self.e[Eta]={},{}
        H,h={},{}
        dynamical_matrices=self.symbolic_system._evaluate_dynamical_matrices(Eta.dict_of_values)
        for contact_point in self.list_of_contact_points:
            H[contact_point],h[contact_point]=contact_point.Evaluate_polyhedral_matrices(dynamical_matrices,Eta.dict_of_values,epsilon_max,epsilon_min)
            self.E[Eta][contact_point],self.e[Eta][contact_point]={},{}
            for sigma in contact_point.Sigma:
                self.E[Eta][contact_point][sigma]=H[contact_point][sigma]
                self.e[Eta][contact_point][sigma]=h[contact_point][sigma]
        self.A[Eta]=dynamical_matrices.A
        self.B_u[Eta]=dynamical_matrices.B_u
        self.B_lambda[Eta]=dynamical_matrices.B_lambda
        self.c[Eta]=dynamical_matrices.c
        
    

class environment:
    def __init__(self,name=0):
        self.dict_of_values={}
        self.name="Env "+str(name)
        
    def __repr__(self):
        return self.name
    
    
def merge_timed_vectors_glue(list_of_timed_dicts):
    X={}
    t0=0
    for D in list_of_timed_dicts:
        T=max(D.keys())
        for t in range(T+1):
            X[t0+t]=D[t]
        t0+=T
    return X

def merge_timed_vectors_skip(list_of_timed_dicts):
    X={}
    t0=0
    for D in list_of_timed_dicts:
        T=max(D.keys())
        for t in range(T+1):
            X[t0+t]=D[t]
        t0+=T+1
    return X

def simulate_system_numeric(sys,x,u,Eta):
    """
    System sys is numeric
    """
    Exu={}
    for i in sys.list_of_contact_points:
        Exu[i]={}
        for sigma in sys.E[Eta][i].keys():
            Exu[i][sigma]=sys.e[Eta][i][sigma]-np.dot(sys.E[Eta][i][sigma][:,0:sys.n+sys.m_u],np.vstack((x,u)))
    print Exu
    raise NotImplementedError        
        