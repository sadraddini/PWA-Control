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
from pypolycontain.lib.polytope import polytope

from pypolycontain.lib.AH_polytope import minimum_distance
from pypolycontain.lib.hausdorff.hausdorff import Hausdorff_directed


from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes_ax as visZ
from pypolycontain.visualization.visualize_2D import visualize_2D_ax as vis

import matplotlib.pyplot as plt

import time




class node:
    def __init__(self,list_of_polytopes,hyperplane=[0,0,0],parent_node=False,depth=0):
        """
        parent_node is None for root
        hyperplane is a tuple (c,g,sign)
        """
        self.list_of_polytopes=list_of_polytopes
        self.hyperplane=hyperplane
        self.parent_node=parent_node
        self.depth=depth
        
    def __repr__(self):
        return "node with %d polytopes at depth %d"%(len(self.list_of_polytopes),self.depth)

    def division(self):
        """
        Divide the node into two nodes based on polytopes inside
        """
        c,g,mu_p,mu_n=create_hyperplane(self.list_of_polytopes)
        node_positive=node([poly for poly in self.list_of_polytopes if mu_p[poly]==1 or mu_n[poly]+mu_p[poly]==0],hyperplane=(c,g,"+"),parent_node=self,depth=self.depth+1)
        node_negative=node([poly for poly in self.list_of_polytopes if mu_n[poly]==1 or mu_p[poly]+mu_n[poly]==0],hyperplane=(c,g,"-"),parent_node=self,depth=self.depth+1)
        return (node_positive,node_negative)




class BSP_tree:
    def __init__(self,list_of_polytopes):
        self.list_of_polytopes=list_of_polytopes
        self.nodes=[node(list_of_polytopes)]
        self.leafs=list(self.nodes)
        
    def __repr__(self):
        return "BSP_tree of AH-polytopes"
    
    def get_dimensions(self):
        return self.list_of_polytopes[0].T.shape[0]
    
    def deepen_tree_one_step(self,N=3):
        new_leafs=[]
        for leaf in self.leafs:
            if len(leaf.list_of_polytopes)<N:
                continue
            else:
                node_positive,node_negative=leaf.division()
                new_leafs.extend([node_positive,node_negative])
                self.nodes.extend([node_positive,node_negative])
        self.leafs=list(new_leafs)
        
    def construct_tree(self,D,N=5,visualize=True):
        for d in range(D):
            self.deepen_tree_one_step(N)
            if visualize:
                self.visualize()
                
    def _get_polyhedrons(self,X=True):
        for leaf in self.leafs:
            leaf.polytope=_get_polyhedron(self,leaf,X)
    
    def _get_minimums(self):
        start=time.time()
        _=minimum_distance(self.list_of_polytopes[0],self.leafs[0].polytope)
        T=time.time()-start
        print "*"*100
        print "Getting Minimums, estimated time: %0.02f seconds"%(T*len(self.list_of_polytopes)*len(self.leafs))
        print "*"*100 
        return {(leaf,P): minimum_distance(leaf.polytope,P) for leaf in self.leafs for P in self.list_of_polytopes}
    
    def _get_maximums(self):
        start=time.time()
        _=Hausdorff_directed(self.list_of_polytopes[0],self.leafs[0].polytope)
        T=time.time()-start
        print "*"*100
        print "Getting Maximums, estimated time: %0.02f seconds"%(T*len(self.list_of_polytopes)*len(self.leafs))
        print "*"*100
        return {(leaf,P): Hausdorff_directed(leaf.polytope,P) for leaf in self.leafs for P in self.list_of_polytopes}
    
    def build_Plist(self):
        D_min=self._get_minimums()
        D_max=self._get_maximums()
        P_list={}
        for leaf in self.leafs:
            d_max_min=min([D_max[leaf,P] for P in self.list_of_polytopes])
            print "dmax_min is",d_max_min
            P_list[leaf]=[poly for poly in self.list_of_polytopes if D_min[(leaf,poly)]<=d_max_min]
        return P_list
        
        
    def visualize(self):
        """
        Assumption: all AH-polytopes are zonotopes, and the problem is in 2D
        """
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

    

class cell:
    def __init__(self,P,list_of_polytopes=[],depth=0):
        self.polytope=P
        self.list_of_polytopes=list_of_polytopes
        self.depth=depth
    
    def __repr__(self):
        return "cell at depth %d with %d polytopes"%(self.depth,len(self.list_of_polytopes))
        
    def cut_half(self):
        """
        Divide the node into two nodes based on the polytope itself
        """
        X=np.hstack([z.x for z in self.list_of_polytopes])
        u,s,v=np.linalg.svd(X)
        c_u=u[random.randint(1,len(u))-1,:].reshape(len(u),1)
        c=c_u+np.random.normal(size=(len(u),1))
        C1,C2,hyperplane=_cut_half(c,self.polytope)
        S={}
        for C in [C1,C2]:
            D_max={P: Hausdorff_directed(C,P) for P in self.list_of_polytopes} 
            d_max_min=min(D_max.values())
            D_min={P: minimum_distance(C,P) for P in self.list_of_polytopes} 
            list_of_polytopes=[P for P in self.list_of_polytopes if D_min[P]<=d_max_min]
            S[C]=cell(C,list_of_polytopes,depth=self.depth+1)
        return S[C1],S[C2],hyperplane
        

class BSP_tree_cells:
    def __init__(self,X,list_of_polytopes):
        self.list_of_polytopes=list_of_polytopes
        self.cells=[cell(X,list_of_polytopes)]
        self.leafs=list(self.cells)
        self.hyperplanes=[]
        
    def __repr__(self):
        return "BSP_tree of AH-polytopes"
    
    def _get_dimensions(self):
        return self.list_of_polytopes[0].T.shape[0]    

    def _deepen_tree_one_step(self,N):
        new_leafs=[]
        for leaf in self.leafs:
            if len(leaf.list_of_polytopes)<=N:
                new_leafs.append(leaf)
                continue
            else:
                cell_positive,cell_negative,hyperplane=leaf.cut_half()
                new_leafs.extend([cell_positive,cell_negative])
                self.cells.extend([cell_positive,cell_negative])
                self.hyperplanes.append(hyperplane)
        self.leafs=list(new_leafs)
        
    def construct_tree(self,D,N=5,visualize=True):
        for d in range(D):
            print "depth:",d
            print [len(leaf.list_of_polytopes) for leaf in self.leafs]
            self._deepen_tree_one_step(N)
            if visualize:
                self.visualize()

    def visualize(self):
        """
        Assumption: all AH-polytopes are zonotopes, and the problem is in 2D
        """
        fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
        zonotopes=[]
        for leaf in self.leafs:
            color=(random.random(),random.random(),random.random())
            zonotopes.extend([zonotope(poly.x,poly.G,color=color) for poly in leaf.list_of_polytopes])
        for hyperplane in self.hyperplanes:
            x_min,x_max=-200,200
            c,g=hyperplane[0:2]
            g=np.asscalar(g)
#            print c,c.shape,g
            ax.plot([x_min,x_max],[(g-x_min*c[0,0])/(c[1,0]+0.001),(g-x_max*c[0,0])/(c[1,0]+0.001)],color=(0.2,0,0), linewidth=1,linestyle='dashed')
        visZ(ax,zonotopes,title="BSP_tree")

    def draw_leaf(self,leaf):
        """
        Assumption: all AH-polytopes are zonotopes, and the problem is in 2D
        """
#        assert leaf in self.cells
        fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
        for P in self.list_of_polytopes:
            if P in leaf.list_of_polytopes:
                P.color=(0,0,1)
            else:
                P.color=(1,0,0)
        vis(ax,[leaf.polytope])
        visZ(ax,self.list_of_polytopes,title="BSP_tree for %s"%leaf.__repr__()) 
        
    def draw_cells(self):
        """
        Assumption: all AH-polytopes are zonotopes, and the problem is in 2D
        """
        fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
        vis(ax,[cell.polytope for cell in self.cells])
    

def _get_polyhedron(tree,leaf,X=True):
    assert leaf in tree.leafs
    n=tree.get_dimensions()
    if X==True:
        H,h=np.empty((0,n)),np.empty((0,1))
    else:
        H,h=X.H,X.h
    u=leaf
    while u.parent_node!=False:
        c,g,s=u.hyperplane[0:3]
        if s=="-":
            H=np.vstack((H,c.T))
            h=np.vstack((h,g))
        if s=="+":
            H=np.vstack((H,-c.T))
            h=np.vstack((h,-g))
        u=u.parent_node
    return polytope(H,h)
        



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

    
        
def _cut_half(c,P):
    """
    Given polytopen and direction c, cut it into half
    """
    n=c.shape[0]
    assert n==P.H.shape[1]
    model=Model("n")
    x=tupledict_to_array(model.addVars(range(n),[0],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="x"))
    model.update()
    constraints_list_of_tuples(model,[(P.H,x),(-np.eye(P.h.shape[0]),P.h)],sign="<")
    J=LinExpr([(c[i,0],x[i,0]) for i in range(n)])
    model.setParam('OutputFlag', False)
    model.setObjective(J,GRB.MINIMIZE)
    model.optimize()
    g_min=model.ObjVal
    model.setObjective(J,GRB.MAXIMIZE)
    # Max
    model.reset()
    model.optimize()
    g_max=model.ObjVal
    g=np.array([(g_max+g_min)/2.0]).reshape(1,1)
    P1=polytope(np.vstack((P.H,c.T)),np.vstack((P.h,g)))
    P2=polytope(np.vstack((P.H,-c.T)),np.vstack((P.h,-g)))
    return P1,P2,(c,g)
    
            
    
    
    
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