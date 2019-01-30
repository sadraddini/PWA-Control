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
    
    def get_contact_free_dynamics(self,x_sample,u_sample):

        N=x_sample.shape[0]
        n,m,nC=len(self.x),len(self.u),len(self.contact_points)
        E,i=[0]*n,0
        for f_c in self.f:
#            print i
            f_x=[diff(f_c,var) for var in self.x]
#            print f_x
            D=self.sys_lambdify(f_x)
            E[i]=self.evaluate_handles(D,x_sample,u_sample)
#            print E
            i+=1
        print E
        A={}
        for k in range(N):
            A[k]=np.zeros((n,n))
            for i in range(n):
                for j in range(n):
                    if type(E[i][j]) in [type(0.0),type(0),type(np.float64(0))]:
                        A[k][i,j]=E[i][j]
                    elif type(E[i][j]) ==type(np.array(0)):
                        A[k][i,j]=E[i][j][k]
                    else:
                        raise NotImplementedError
        c_lambda=self.sys_lambdify(self.f)
        c=self.evaluate_handles(c_lambda,x_sample,u_sample)
        print c
        return A,c
