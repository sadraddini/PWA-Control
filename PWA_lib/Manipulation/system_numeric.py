#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 16:06:32 2019

@author: sadra
"""

class system_HACTS:
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

       
    def add_environment(self,Eta,epsilon):
        self.Eta.append(Eta)
        self.E[Eta],self.e[Eta]={},{}
        H,h={},{}
        dynamical_matrices=self.symbolic_system._evaluate_dynamical_matrices(Eta.dict_of_values)
        for contact_point in self.list_of_contact_points:
            H[contact_point],h[contact_point]=contact_point.Evaluate_polyhedral_matrices(dynamical_matrices,Eta.dict_of_values,epsilon)
            for sigma in contact_point.Sigma:
                self.E[Eta][contact_point,sigma]=H[contact_point][sigma]
                self.e[Eta][contact_point,sigma]=h[contact_point][sigma]
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