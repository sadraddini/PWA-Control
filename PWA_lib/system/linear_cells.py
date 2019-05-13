#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 12:37:25 2019

@author: sadra
"""

"""
The class of linear cells is used for local constrained linear systems
"""

class linear_cell:
    """
    This is a linear cell
    """
    def __init__(self,A,B,c,p):
        self.A=A
        self.B=B
        self.c=c
        self.p=p
        self.name=None
        
    def __repr__(self):
        if self.name==None:
            return "A linear cell with %d states and %d controls"%(self.B.shape[0],self.B.shape[1])
        else:
            return "The linear cell of %s"%str(self.name)
        
class disturbed_linear_cell:
    """
    This is a linear cell for systems with disturbances
    """
    def __init__(self,A,B,w,p_x,p_u):
        self.A=A
        self.B=B
        self.w=w # Must be a Zonotope
        self.p_x=p_x # A polytope for state
        self.p_u=p_u # A polytope for controls
        
    def __repr__(self):
        return "A disturbed linear cell with %d states and %d controls"%(self.B.shape[0],self.B.shape[1])