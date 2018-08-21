#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 15:54:23 2018

@author: sadra
"""
# Primary imports
import numpy as np
from random import choice as rchoice
from random import random
from time import time

### Internal imports
import sys
sys.path.append('..')

# Secondary imports
from main.auxilary_methods import find_mode,valuation,mode_sequence
from main.ana_system import state
from main.trajectory import subset_MILP


def reshape_tree(s):
    for y in s.X:
        if len(y.parent)>1:
            Group={}
            for i in s.modes:
                Group[i]=[state_y for state_y in y.parent if state_y.mode==i]
                if len(Group[i]>1):
                    pass

def regroup_states(s,list_of_states):
    for y in list_of_states:
        s.X.remove(y)
    new_state=state()