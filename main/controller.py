#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 10:51:52 2018

@author: sadra
"""
import numpy as np
from main.tree_locator import tree_locator

def control_vanilla(s,x):
    x_state=tree_locator(s,x)
    print("tree node is",x_state)
    if x_state==s.goal:
        print("*"*40,"YUHOOO! REACHED THE GOAL","*"*40)
        return False
    u_0=x_state.successor[1]
    u_f=np.dot(x_state.successor[2],np.dot(x_state.Ginv,x-x_state.x))
    print("u_0=",u_0.T)
    print("u_f=",u_f.T)
    return u_0+u_f