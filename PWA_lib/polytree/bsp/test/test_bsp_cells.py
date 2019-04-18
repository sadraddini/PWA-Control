#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 15:49:25 2019

@author: sadra
"""

import numpy as np

from pypolycontain.lib.AH_polytope import AH_polytope
from pypolycontain.lib.polytope import polytope,Box,translate
from pypolycontain.lib.zonotope import zonotope

from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes_ax as visZ
from pypolycontain.visualization.visualize_2D import visualize_2D_ax as vis


from PWA_lib.polytree.bsp.bsp import BSP_tree_cells

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


n=2
q=4
B=Box(q)
N=200
X=translate(Box(2,12),np.array([10,10]).reshape(2,1))
#list_of_zonotopes=[zonotope(np.random.random((n,1))*10+np.array([0.1*i,0.5*i]).reshape(2,1),np.random.random((n,q))*1) for i in range(N)]
list_of_zonotopes=[zonotope(np.random.random((n,1))*20*np.array([1,2]).reshape(2,1),np.random.random((n,q))*2-1) for i in range(N)]

fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
visZ(ax,list_of_zonotopes)


mytree=BSP_tree_cells(X,list_of_zonotopes)
mytree.construct_tree(D=10  ,N=30)