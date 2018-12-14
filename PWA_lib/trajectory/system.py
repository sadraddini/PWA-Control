# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 17:21:45 2018

@author: sadra
"""

import numpy as np

#import pickle

class system:
    def __init__(self):
        print "A new system is created"
        self.name="system"
        self.A={} # Matrix A
        self.B={} # Matrix B
        self.c={} # vector c
        self.C={} # Cell polytopes
        self.goal=None
    
    def __repr__(self):
        return self.name
        
    def build(self,T=50):
#        self.C_dict={(t,n,i,j,k):self.C[n,i].H[j,k] for t in range(T)
#        for n in self.list_of_sum_indices for i in self.list_of_modes[n]
#        for j in range(self.C[n,i].H.shape[0]) for k in range(self.n)}
        pass
    

class linear_cell:
    def __init__(self,A,B,c,p):
        self.A=A
        self.B=B
        self.c=c
        self.p=p