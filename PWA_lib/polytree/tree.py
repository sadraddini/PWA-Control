#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 13:05:19 2019

@author: sadra
"""

# External Imports
import numpy as np

# Internal Imports
from pypolycontain.lib.zonotope import zonotope,zonotope_distance_point
from pypolycontain.utils.utils import vertices_cube as vcube 
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes_ax as visZ,plt

from PWA_lib.trajectory.poly_trajectory import point_trajectory_sPWA,polytopic_trajectory_given_modes


class state:
    def __init__(self,p,branch=0,t=0):
        self.name="state branch=%d t=%d"%(branch,t)
        self.p=p
        self.x=p.x
        self.G=p.G+np.random.random(p.G.shape)*0.0001
        self.p=zonotope(self.x,self.G,color=p.color,name="zonotope"+self.name)
        self.t=t
        v=vcube(p.G.shape[1])
        self.vertices=p.x.T+(np.dot(p.G,v.T)).T
        self.cost_to_go=0
        self.time_to_go=0
        self.cost_to_child=0
        self.parents=[]
        self.child=None
        self.control=None
            
    def __repr__(self):
        return "%s - %d parents" %(self.name,len(self.parents))

    
class tree():
    # Tree is defined for a system
    def __init__(self,sys):
        self.system=sys
        self.goal=state(sys.goal)
        self.states=[]
        self.branches=[]
        self.fill="discrete"
        self.name=sys.name
        self.fig,self.ax=plt.subplots()
        
    def __repr__(self):
        return "** system: "+ self.system.name + " **: %d states and %d branches" %(len(self.states),len(self.branches))
        
    def restart_fig(self):
        self.fig,self.ax=plt.subplots()
        
    def visualize(self,axis_limit=[True]):
        visZ(self.ax,[x.p for x in self.states],axis_limit=axis_limit,title="%s - %d zonotopes %d branches" %(self.name,len(self.states),len(self.branches)))
        self.fig.savefig("figures/tree_%d"%len(self.branches))
            
    def add_branch(self,x,u,G,theta,goal):
        """
        add x,u,G,theta
        """
        T=len(u)
        Z={}
        import random as rand
        c=(rand.random(),rand.random(),rand.random()) #over-ride previous one
        for t in range(T):
#            c=(t/(T+0.01), 1-t/(T+0.01), 0.5)
            if self.fill=="discrete":
                p=zonotope(x[t],G[t],color=c)
                Z[t]=state(p,len(self.branches),t)
            elif self.fill=="continous":
                G_c=np.hstack( ( G[t], (x[t+1]-x[t])/1.0 ) )
                p_c=zonotope((x[t]+x[t+1])/2.0,G_c,color=c)
                Z[t]=state(p_c,len(self.branches),t)
            else:
                raise self.fill, ": not understood, please enter 'discrete' or 'continous' for tree.fill"
        Z[T+1]=state(goal)
        for t in range(T-1):
            Z[t].child=Z[t+1]
            Z[t].control=(u[t],theta[t])
            Z[t+1].parents.append(Z[t])
        self.states.extend([Z[t] for t in range(T)])
        self.branches.append(Z)
        
    def extend_RRT(self,x_sample,T,eps_point=1,eps_poly=1):
        sys=self.system
        if self.states!=[]:
            list_of_goals=[x.p for x in self.states]
        else:
            print "WARNING: no gaol. Goal fixed to system goal. Is this the initizaliztion?"
            list_of_goals=[self.system.goal]
        (x_nom,u_nom,delta_PWA,mu_nom,flag)=point_trajectory_sPWA(self.system,x_sample,list_of_goals,T,eps=eps_point)
        if flag==False:
            return False
        else:
            list_of_cells=[]
            for t in range(T):
                mode=tuple([i for n in sys.list_of_sum_indices for i in sys.list_of_modes[n] if delta_PWA[t,n,i]==1])    
                list_of_cells.append(sys.cell[mode])
    #            print t,"mode",mode,x_nom[t].T
    #        verify(x_nom,u_nom,delta_PWA,list_of_cells)
            goal=None
            for _goal in list_of_goals:
                if mu_nom[_goal]==1:
                    goal=_goal
            print "goal is",goal,goal.x.T,"G=",goal.G,"x_nom[T]",x_nom[T].T
    #        print "mu=",mu_nom
            assert goal!=None #Something is found!
            (x,u,G,theta)=polytopic_trajectory_given_modes(x_nom[0],list_of_cells,goal,eps=eps_poly,scale=sys.scale)
            self.add_branch(x,u,G,theta,goal)
            return True



    def construct_tree(self,list_of_samples,T_start=10,T_max=50,eps_point=0.1,eps_poly=0.1):
        for x_sample in list_of_samples:
            if self.inside_tree(x_sample):
                print "sample already inside the tree"
                continue
            else:
                T=T_start
                while self.extend_RRT(x_sample,T,eps_point=eps_point,eps_poly=eps_poly)==False and T<T_max:
                    print "T=",T
                    T+=1
        print "*** End of Construction of Tree"
        
    def inside_tree(self,x,eps=10**-4):
        for state in self.states:
            d=zonotope_distance_point(state.p,x)
            if d<eps:
                return True
        return False

def tree_locate(tree,x):
    d={}
    for state in tree.states:
        d[state]=zonotope_distance_point(state.p,x)
    return d

    
            
            
        
        
"""
Not really needed functions:
"""        



def verify(x,u,delta_PWA,list_of_cells):
    T=len(u)
    for t in range(T):
        p=list_of_cells[t].p
        A=list_of_cells[t].A
        B=list_of_cells[t].B
        c=list_of_cells[t].c
        e_p=p.h-np.dot(p.H,np.vstack((x[t],u[t])))
        e_x=abs(x[t+1]-np.dot(A,x[t])-np.dot(B,u[t])-c)
        print t,np.amin(e_p)<0
        print t,np.amax(e_x)>10**-8
        