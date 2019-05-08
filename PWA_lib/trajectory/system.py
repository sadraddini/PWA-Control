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
    
    
    def build_cells(self):
        from pypolycontain.utils.redundancy_reduction import canonical_polytope
        from pypolycontain.lib.polytope import polytope
        from itertools import product
        all_modes=list(product(*[self.list_of_modes[n] for n in self.list_of_sum_indices]))
        self.cell={}
        print all_modes
        for mode in all_modes:
            delta_PWA={}
            for n in self.list_of_sum_indices:
                for j in range(len(self.list_of_sum_indices)):
                    if self.list_of_sum_indices[j]==n:
                        break
                for i in self.list_of_modes[n]:
#                    delta_PWA[n,i]=int(mode[n]==i)
                    delta_PWA[n,i]=int(mode[j]==i)
            print delta_PWA
            A=sum([self.A[n,i]*delta_PWA[n,i] for n in self.list_of_sum_indices for i in self.list_of_modes[n]])
            B=sum([self.B[n,i]*delta_PWA[n,i] for n in self.list_of_sum_indices for i in self.list_of_modes[n]])
            c=sum([self.c[n,i]*delta_PWA[n,i] for n in self.list_of_sum_indices for i in self.list_of_modes[n]])
            H=np.vstack([self.C[n,i].H for n in self.list_of_sum_indices for i in self.list_of_modes[n] if delta_PWA[n,i]==1])
            h=np.vstack([self.C[n,i].h for n in self.list_of_sum_indices for i in self.list_of_modes[n] if delta_PWA[n,i]==1])
#            H,h=canonical_polytope(H,h)
            cell=linear_cell(A,B,c,polytope(H,h))    
            self.cell[mode]=(cell)
            
    def simulate_forward():
        pass

class linear_cell:
    """
    This is a linear cell
    """
    def __init__(self,A,B,c,p):
        self.A=A
        self.B=B
        self.c=c
        self.p=p
        
    def __repr__(self):
        return "Linear cells with %d states and %d controls"%(self.B.shape[0],self.B.shape[1])
        
class disturbed_linear_cell:
    """
    This is a linear cell
    """
    def __init__(self,A,B,w,p_x,p_u):
        self.A=A
        self.B=B
        self.w=w # Must be a Zonotope
        self.p_x=p_x # A polytope
        self.p_u=p_u # A polytope
