#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 19:52:33 2019

@author: sadra
"""

import numpy as np
import random

from scipy.optimize import linprog as lp
from gurobipy import Model,LinExpr,QuadExpr,GRB

from pypolycontain.utils.utils import valuation
from pypolycontain.lib.zonotope import zonotope
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes_ax as visZ
import matplotlib.pyplot as plt


class node:
    def __init__(self,list_of_polytopes,hyperplane=[0,0,0],parent_node=False):
        """
        parent_node is None for root
        hyperplane is a tuple (c,g,sign)
        """
        self.list_of_polytopes=list_of_polytopes
        self.hyperplane=hyperplane
        self.parent_node=parent_node
        
    def __repr__(self):
        return "node with %d polytopes"%len(self.list_of_polytopes)

    def division(self):
        c,g,mu_p,mu_n=create_hyperplane(self.list_of_polytopes)
        node_positive=node([poly for poly in self.list_of_polytopes if mu_p[poly]==1 or mu_n[poly]+mu_p[poly]==0],hyperplane=(c,g,"+"),parent_node=self)
        node_negative=node([poly for poly in self.list_of_polytopes if mu_n[poly]==1 or mu_p[poly]+mu_n[poly]==0],hyperplane=(c,g,"-"),parent_node=self)
        return (node_positive,node_negative)
    
    
    

        
class BSP_tree:
    def __init__(self,list_of_polytopes):
        self.list_of_polytopes=list_of_polytopes
        self.nodes=[node(list_of_polytopes)]
        self.leafs=list(self.nodes)
        
    def __repr__(self):
        return "BSP_tree of AH-polytopes"
    
    def deepen_tree_one_step(self):
        new_leafs=[]
        for leaf in self.leafs:
            node_positive,node_negative=leaf.division()
            new_leafs.extend([node_positive,node_negative])
            self.nodes.extend([node_positive,node_negative])
        self.leafs=list(new_leafs)
        
    def visualize(self):
        """
        Assumption: all AH-polytopes are zonotopes, and the problem is in 2D
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
        zonotopes=[]
        for leaf in self.leafs:
            color=(random.random(),random.random(),random.random())
            zonotopes.extend([zonotope(poly.t,poly.T,color=color) for poly in leaf.list_of_polytopes])
        for leaf in self.nodes:
            if leaf.parent_node==False:
                continue
            elif leaf.parent_node.hyperplane[2]==0:
                x_min,x_max=-200,200
                c,g=leaf.hyperplane[0:2]
                ax.plot([x_min,x_max],[(g-x_min*c[0,0])/(c[1,0]+0.001),(g-x_max*c[0,0])/(c[1,0]+0.001)],color=(0.2,0,0), linewidth=2,linestyle='dashed')
            else:              
                if leaf.parent_node.hyperplane[2]=="+":
#                    x_p,y_p=find_the_intersection_of_hyperplanes(leaf.hyperplane,leaf.parent_node.hyperplane)
                    x_min,x_max=-200,200
                    c,g=leaf.hyperplane[0:2]
#                    _p=np.array([-200,(g+200*c[0,0])/(c[1,0]+0.001)]).reshape(2,1)
#                    _q=np.array([200,(g-200*c[0,0])/(c[1,0]+0.001)]).reshape(2,1)
#                    if np.asscalar(np.dot(_q.T,c))>=g:
#                        x_min=x_p
#                    else:
#                        x_max=x_p
                    ax.plot([x_min,x_max],[(g-x_min*c[0,0])/(c[1,0]+0.001),(g-x_max*c[0,0])/(c[1,0]+0.001)],color=(0.0,0.5,0), linewidth=1,linestyle='dashed')  
                elif leaf.parent_node.hyperplane[2]=="-":
#                    x_p,y_p=find_the_intersection_of_hyperplanes(leaf.hyperplane,leaf.parent_node.hyperplane)
                    x_min,x_max=-200,200
                    c,g=leaf.hyperplane[0:2]
#                    _p=np.array([-200,(g+200*c[0,0])/(c[1,0]+0.001)]).reshape(2,1)
#                    _q=np.array([200,(g-200*c[0,0])/(c[1,0]+0.001)]).reshape(2,1)
#                    if np.asscalar(np.dot(_p.T,c))<=g:
#                        x_max=x_p
#                    else:
#                        x_min=x_p
                    ax.plot([x_min,x_max],[(g-x_min*c[0,0])/(c[1,0]+0.001),(g-x_max*c[0,0])/(c[1,0]+0.001)],color=(0.0,0,0.5), linewidth=1,linestyle='dashed') 
        visZ(ax,zonotopes,title="Random Zonotopes")


def find_the_intersection_of_hyperplanes(hyperplane1,hyperplane2):
    c1,g1,c2,g2=hyperplane1[0:2]+hyperplane2[0:2]
    C=np.vstack((c1.T,c2.T))
    g=np.array([g1,g2]).reshape(2,1)
    p=np.dot(np.linalg.inv(C),g).reshape(2)
    return p[0],p[1]
    
        
        
            
        

def create_hyperplane(list_of_polytopes):
    """
    Given a list of polytopes, create a hyperplane that roughly cuts the objects half in group
    J=(mu_pos-N/2)**2+(mu_neg-N/2)**2
    """
    model=Model("Hyperplane generation")
    n=list_of_polytopes[0].t.shape[0]
    c=tupledict_to_array(model.addVars(range(n),[0],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="c"))
    g=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
    u={}
    epsilon={}
    mu={}
    bigM=100
    eps_tol=0.001
    for poly in list_of_polytopes:
        for sign in ["+","-"]:
            u[poly,sign]=tupledict_to_array(model.addVars(range(poly.P.h.shape[0]),[0],lb=0,ub=GRB.INFINITY,name="u"))
            epsilon[poly,sign]=model.addVar(lb=-bigM/2,ub=bigM/2,name="eps")
            mu[poly,sign]=model.addVar(vtype=GRB.BINARY,name="mu(-)")
    mu["sum","+"]=model.addVar(lb=0,name="sum+")
    mu["sum","-"]=model.addVar(lb=0,name="sum-")
    model.update()
    for poly in list_of_polytopes:
        constraints_list_of_tuples(model,[(poly.P.h.T,u[poly,"+"]),
                                          (np.eye(1),np.array([g]).reshape(1,1)),
                                          (np.eye(1),np.array([epsilon[poly,"+"]]).reshape(1,1)),
                                          (-poly.t.T,c)],sign="=")
        constraints_list_of_tuples(model,[(poly.P.h.T,u[poly,"-"]),
                                          (-np.eye(1),np.array([g]).reshape(1,1)),
                                          (np.eye(1),np.array([epsilon[poly,"-"]]).reshape(1,1)),
                                          (poly.t.T,c)],sign="=")    
        constraints_list_of_tuples(model,[(poly.P.H.T, u[poly,"+"]),(poly.T.T,c)],sign="=")     
        constraints_list_of_tuples(model,[(poly.P.H.T, u[poly,"-"]),(-poly.T.T,c)],sign="=")     
        # bigM_formulation
        for sign in ["+","-"]:
            model.addConstr(bigM*mu[poly,sign]>=epsilon[poly,sign]-eps_tol)
            model.addConstr(bigM-bigM*mu[poly,sign]>=-epsilon[poly,sign]+eps_tol)
    # sum of mu values
    model.addConstr(mu["sum","+"]==LinExpr([(1,mu[poly,"+"]) for poly in list_of_polytopes]))
    model.addConstr(mu["sum","-"]==LinExpr([(1,mu[poly,"-"]) for poly in list_of_polytopes]))
    # The first value of c is one
    model.addConstr(c[0,0]==1)
    # Cost objective
    J=QuadExpr()
    N=len(list_of_polytopes)+0.0
    J.add(mu["sum","+"]*mu["sum","+"]+mu["sum","-"]*mu["sum","-"]-N*mu["sum","+"]-N*mu["sum","-"]+N**2/2.0-mu["sum","+"]-mu["sum","-"])
    J=J*100/(N**2)
    model.setObjective(J,GRB.MINIMIZE)
    model.optimize()
#    for poly in list_of_polytopes:
#        print poly
#        print poly.t.T,"mu+:",mu[poly,"+"].X,"mu-:",mu[poly,"-"].X,"eps+:",epsilon[poly,"+"].X,"eps-:",epsilon[poly,"-"].X
#        print "u+=",valuation(u[poly,"+"]).T,"u-=",valuation(u[poly,"-"]).T
#        print "+",np.dot(valuation(u[poly,"+"]).T,poly.P.h)+np.dot(valuation(c).T,poly.t)-g.X#-epsilon[poly,"+"].X
#        print "-",np.dot(valuation(u[poly,"-"]).T,poly.P.h)-np.dot(valuation(c).T,poly.t)+g.X#-epsilon[poly,"+"].X
    mu_n_positive={poly:np.round(mu[poly,"+"].X) for poly in list_of_polytopes}
    mu_n_negative={poly:np.round(mu[poly,"-"].X) for poly in list_of_polytopes}
    print "sums",mu["sum","-"].X,mu["sum","+"].X
    return valuation(c),g.X,mu_n_positive,mu_n_negative

    
        
    
    
    
    
    
"""
Auxilary Gurobi Shortcut Functions
used for convenience
"""

def tupledict_to_array(mytupledict):
    # It should be 2D
    n,m=max(mytupledict.keys())
    n+=1
    m+=1
    array=np.empty((n,m),dtype="object")
    for i in range(n):
        for j in range(m):
            array[i,j]=mytupledict[i,j]
    return array

def constraints_list_of_tuples(model,mylist,sign="="):
    term_0=mylist[0]
    ROWS,COLUMNS=term_0[0].shape[0],term_0[1].shape[1]
    for row in range(ROWS):
        for column in range(COLUMNS):
            expr=LinExpr()
            for term in mylist:
                q,qp=term[0].shape[1],term[1].shape[0]
                if q!=qp:
                    raise ValueError(term,"q=%d qp=%d"%(q,qp))
                if type(term[1][0,column])==type(model.addVar()):
                    expr.add(LinExpr([(term[0][row,k],term[1][k,column]) for k in range(q)]))
                elif type(term[0][row,0])==type(model.addVar()):
                    expr.add(LinExpr([(term[1][k,column],term[0][row,k]) for k in range(q)]))
                else:
                    expr.addConstant(sum([term[1][k,column]*term[0][row,k] for k in range(q)]))
            if sign=="<":
                model.addConstr(expr<=0)
            elif sign=="=":
                model.addConstr(expr==0)
            elif sign==">=":
                model.addConstr(expr>=0)
            else:
                raise "sign indefinite"