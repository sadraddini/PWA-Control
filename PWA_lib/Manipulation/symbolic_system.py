#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 13:52:25 2019

@author: sadra
"""

from sympy import Symbol,pi,sin,cos,Function,diff,Matrix,symbols,lambdify
import numpy as np

class symbolic_system:
    def __init__(self,name):
        self.name=name
        self.t=Symbol('t') # Symbol time, used for derivatives
        self.x=[] # list of system states
        self.u=[] # list of system control inputs
        self.contact_points=[] # list of contact points
        
    def __repr__(self):
        return self.name+" with %d states, %d controls, and %d contact points"\
            %(len(self.x),len(self.u),len(self.contact_points))
        
    def sub_derivative(phi):
        raise NotImplementedError
        
    def sys_lambdify(self,list_of_expressions):
        return [lambdify(self.x+self.u,expression,"numpy") for expression in list_of_expressions]
    
    def evaluate_handles(self,D,x,u):
        """
        D: a list of handles
        x: N*n list of states
        u: N*m list of corresponding control inputs
        """
        assert x.shape[1]==len(self.x)
        assert u.shape[1]==len(self.u)        
        xu=np.hstack((x,u))
        G,i={},0
        for handle in D:
            G[i]=handle(*[xu[:,k] for k in range(len(self.x)+len(self.u))])
            i+=1
        return G
    
    def _extract(self,q):
        n=len(self.x)
        J_n=np.array(q[0:n]).reshape(n,1)
        J_t=np.array(q[n:2*n]).reshape(n,1)
        J_n_x=np.array(q[2*n:3*n]).reshape(n,n)
        J_t_x=np.array(q[3*n:4*n]).reshape(n,n)
        return J_n,J_t,J_n_x,J_t_x