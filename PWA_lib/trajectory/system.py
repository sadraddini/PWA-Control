# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 17:21:45 2018

@author: sadra
"""

import numpy as np

import pickle

def system():
    def __init__(self,goal=None,name="system"):
        self.name=name
        self.A={} # Matrix A
        self.B={} # Matrix B
        self.c={} # vector c
        self.H={} # vector d
        self.goal=goal
    
    def __repr__(self):
        return self.name