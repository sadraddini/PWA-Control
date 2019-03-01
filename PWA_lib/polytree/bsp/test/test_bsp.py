#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 16:47:31 2019

@author: sadra
"""

import numpy as np

from pypolycontain.lib.AH_polytope import AH_polytope
from pypolycontain.lib.polytope import polytope,Box,translate
from pypolycontain.lib.zonotope import zonotope

from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes_ax as visZ

from PWA_lib.polytree.bsp.bsp import create_hyperplane,BSP_tree

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle





n=2
q=3
B=Box(q)
N=100
X=translate(Box(2,N/2.0+2),np.array([N/2,N/2]).reshape(2,1))
list_of_polytopes=[AH_polytope(np.random.random((n,q))*1,np.random.random((n,1))*10+1*i,B) for i in range(N)]
zonotopes={poly:zonotope(poly.t,poly.T) for poly in list_of_polytopes}
fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
visZ(ax,zonotopes.values())

mytree=BSP_tree(list_of_polytopes)

mytree.construct_tree(6,N=3)

mytree._get_polyhedrons(X)
P_list=mytree.build_Plist()

for leaf in mytree.leafs:
    print leaf,len(P_list[leaf])