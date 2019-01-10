#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 13:05:19 2019

@author: sadra
"""

# External Imports
import numpy as np

# Internal Imports
from pypolycontain.lib.zonotope import zonotope
from pypolycontain.utils.utils import vertices_cube as vcube 
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ

from PWA_lib.trajectory.poly_trajectory import point_trajectory,polytopic_trajectory_given_modes

class state:
    def __init__(self,p,branch=0,t=0):
        self.name="state branch=%d t=%d"%(branch,t)
        self.p=p
        self.x=p.x
        self.G=p.G+np.random.random(p.G.shape)*0.0001
        self.p=zonotope(self.x,self.G,color=p.color)
        self.t=t
        v=vcube(p.G.shape[1])
        self.vertices=p.x.T+(np.dot(p.G,v.T)).T
        self.cost_to_go=0
        self.time_to_go=0
        self.cost_to_child=0
        self.parent=[]
        self.child=None
        self.control=None
            
    def __repr__(self):
        return self.name

    
class tree():
    # Tree is defined for a system
    def __init__(self,sys):
        self.system=sys
        self.goal=state(sys.goal)
        self.states=[state(sys.goal)]
        self.branches=[]
        self.fill="discrete"
        
    def __repr__(self):
        return "** system: "+ self.system.name + " **: %d states and %d branches" %(len(self.states),len(self.branches))
        
        
    def visualize(self,axis_limit=[None]):
        visZ([x.p for x in self.states],axis_limit=axis_limit)
    
    def add_branch(self,x,u,G,theta,goal):
        """
        add x,u,G,theta
        """
        T=len(u)
        Z={}
        for t in range(T):
            if self.fill=="discrete":
                p=zonotope(x[t],G[t],color=(t/(T+0.01), 1-t/(T+0.01), 0.5))
                Z[t]=state(p,len(self.branches),t)
            elif self.fill=="continous":
                G_c=np.hstack( ( (G[t]+G[t+1])/2.0, (x[t+1]-x[t])/2.0 ) )
                p_c=zonotope((x[t]+x[t+1])/2.0,G_c,color=(t/(T+0.01), 1-t/(T+0.01), 0.5))
                Z[t]=state(p_c,len(self.branches),t)
            else:
                raise self.fill, ": not understood, please enter 'discrete' or 'continous' for tree.fill"
        Z[T+1]=state(goal)
        for t in range(T-1):
            Z[t].child=Z[t+1]
            Z[t].control=(u[t],theta[t])
            Z[t+1].parent=Z[t]
        self.states.extend(Z.values())
        self.branches.append(Z)
        
    def extend_RRT(self,x_sample,T,eps=1):
        sys=self.system
        (x_nom,u_nom,delta_PWA,mu_nom)=point_trajectory(self.system,x_sample,[x.p for x in self.states],T,eps=[None])
        list_of_cells=[]
        for t in range(T):
            mode=tuple([i for n in sys.list_of_sum_indices for i in sys.list_of_modes[n] if delta_PWA[t,n,i]==1])    
            list_of_cells.append(sys.cell[mode])
        goal=None
        for x in self.states:
            if mu_nom[x.p]==1:
                goal=x.p
        assert goal!=None #Something is found!
        (x,u,G,theta)=polytopic_trajectory_given_modes(x_sample,list_of_cells,sys.goal,eps,order=1,scale=sys.scale)
        self.add_branch(x,u,G,theta,goal)