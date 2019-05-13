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
from pypolycontain.visualization.visualize_2D import visualize_2D as vis


from PWA_lib.polytree.bsp.bsp import BSP_tree_cells

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


n=2
q=3
B=Box(q)
N=200
X_1=translate(Box(n,12),np.ones((n,1))*10)
H=np.array([[1,0],[0,1],[-1,0],[0,-1],[1,1]]).reshape(5,2)
h=np.array([42,5,-2,-2,35]).reshape(5,1)
X_1=polytope(H,h)
#X_2=translate(Box(n,5),np.ones((n,1))*15)
#X_3=translate(Box(n,4),np.array([[15,5]]).reshape(2,1))
#X_3=translate(Box(n,4),np.array([[15,5]]).reshape(2,1))

s=2
#list_of_zonotopes=[zonotope(np.random.random((n,1))*10+np.array([0.1*i,0.5*i]).reshape(2,1),np.random.random((n,q))*1) for i in range(N)]
list_of_zonotopes=[zonotope(np.random.random((n,1))*np.array([40,3]).reshape(2,1),np.random.random((n,q))*s-s/2) for i in range(N)]
for Z in list_of_zonotopes:
    Z.J=np.random.random()

fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
visZ(ax,list_of_zonotopes)


mytree=BSP_tree_cells([X_1],list_of_zonotopes)
mytree.construct_tree(D=10,N=5)
mytree.draw_cells(alpha=0.8)

x=np.array([8,3]).reshape(2,1)
(Z,D)=mytree.query_closest(x)
C=mytree._query_find_cell(x)
mytree.draw_leaf(C)