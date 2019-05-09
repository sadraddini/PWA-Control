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

from pypolycontain.lib.AH_polytope import minimum_distance,is_inside,distance_point
from pypolycontain.lib.hausdorff.hausdorff import Hausdorff_directed


from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes_ax as visZ
from pypolycontain.visualization.visualize_2D import visualize_2D_ax as vis

import matplotlib.pyplot as plt

import time


    

class cell:
    """
    Each cell is a polytope in the space Partition
    """
    def __init__(self,P,list_of_polytopes=[],depth=0,hyperplane=None,parent=None):
        self.polytope=P
        self.list_of_polytopes=list_of_polytopes
        self.depth=depth
        self.hyperplane=hyperplane
        self.parent=parent
        self.child_left=None # Positive c and g
        self.child_right=None # Negative c and g
    
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
        P1,P2,(c,g)=_cut_half(c,self.polytope)
        S={}
        self.hyperplane=(c,g)
        for C in [P1,P2]:
            D_max={P: Hausdorff_directed(C,P) for P in self.list_of_polytopes} 
            d_max_min=min(D_max.values())
            D_min={P: minimum_distance(C,P) for P in self.list_of_polytopes} 
            list_of_polytopes=[P for P in self.list_of_polytopes if D_min[P]<=d_max_min]
            S[C]=cell(C,list_of_polytopes,depth=self.depth+1,parent=self)
        self.child_left,self.child_right=S[P1],S[P2]
        return S[P1],S[P2],(c,g)
    
    def shorten_list_by_J(self,tol=10**-4):
        list_of_inside_polytopes=[]
        for P in self.list_of_polytopes:
            if Hausdorff_directed(self.polytope,P)<=tol:
                list_of_inside_polytopes.append(P)
        if len(list_of_inside_polytopes)<=1:
            return # List is empty or singular. Just return
        J_min=10*5
        P_best=None
        for P in list_of_inside_polytopes:
            if J_min>P.J:
                J_min=P.J
                P_best=P
        new_list=[P_best]+[P for P in self.list_of_polytopes if P not in list_of_inside_polytopes]
        self.list_of_polytopes=new_list
        

class BSP_tree_cells:
    def __init__(self,list_of_X,list_of_polytopes):
        self.list_of_polytopes=list_of_polytopes
        self.cells=[cell(X,list_of_polytopes) for X in list_of_X]
        self.leafs=list(self.cells)
        self.hyperplanes=[]
        self.list_of_X=list_of_X
        
    def __repr__(self):
        return "BSP_tree of AH-polytopes with %d-depth, %d leafs and %d cells"%(max([leaf.depth for leaf in self.leafs]) , len(self.leafs), len(self.cells))
    
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
                # Shorten Lists
                cell_positive.shorten_list_by_J()
                cell_negative.shorten_list_by_J()
                # End of Shortening Lists
                new_leafs.extend([cell_positive,cell_negative])
                self.cells.extend([cell_positive,cell_negative])
                self.hyperplanes.append(hyperplane)
        self.leafs=list(new_leafs)
        
    def construct_tree(self,D,N,visualize=False):
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
        
    def draw_cells(self,alpha=0.6):
        """
        Assumption: all AH-polytopes are zonotopes, and the problem is in 2D
        """
        fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
        fig.set_size_inches(10, 10)
        fig.gca().set_aspect('equal')
        vis(ax,[cell.polytope for cell in self.cells],alpha=alpha)
        
        
    def insert_polytope(self):
        """
        Insert a new polytope
        """
        raise NotImplementedError
        
    def query_closest(self,x,norm="L2"):
        d_star=10**12
        P_star=None
        list_to_check=self._query_find_cell(x).list_of_polytopes
        print "size of list to check is %d"%len(list_to_check)
        for P in list_to_check:
            d=distance_point(P,x,norm=norm)
            if d<d_star:
                d_star=d
                P_star=P
        return P_star,d_star
    
    def _query_find_cell(self,x):
        # First find the X
        ana_cell=None
        for k in range(len(self.list_of_X)):
            if is_inside(self.cells[k].polytope,x):
                ana_cell=self.cells[k]
                break
        assert ana_cell!=None
        # Now go down the tree!
        while ana_cell.hyperplane!=None:
            (c,g)=ana_cell.hyperplane
            if all(np.dot(c.T,x)<=g):
                ana_cell=ana_cell.child_left
            else:
                ana_cell=ana_cell.child_right
        return ana_cell
        
    

#def _get_polyhedron(tree,leaf,X=True):
#    assert leaf in tree.leafs
#    n=tree.get_dimensions()
#    if X==True:
#        H,h=np.empty((0,n)),np.empty((0,1))
#    else:
#        H,h=X.H,X.h
#    u=leaf
#    while u.parent_node!=False:
#        c,g,s=u.hyperplane[0:3]
#        if s=="-":
#            H=np.vstack((H,c.T))
#            h=np.vstack((h,g))
#        if s=="+":
#            H=np.vstack((H,-c.T))
#            h=np.vstack((h,-g))
#        u=u.parent_node
#    return polytope(H,h)
#        
#
#
#
#def find_the_intersection_of_hyperplanes(hyperplane1,hyperplane2):
#    c1,g1,c2,g2=hyperplane1[0:2]+hyperplane2[0:2]
#    C=np.vstack((c1.T,c2.T))
#    g=np.array([g1,g2]).reshape(2,1)
#    p=np.dot(np.linalg.inv(C),g).reshape(2)
#    return p[0],p[1]
    
        
        
            
        

#def create_hyperplane(list_of_polytopes):
#    """
#    Given a list of polytopes, create a hyperplane that roughly cuts the objects half in group
#    J=(mu_pos-N/2)**2+(mu_neg-N/2)**2
#    """
#    model=Model("Hyperplane generation")
#    n=list_of_polytopes[0].t.shape[0]
#    c=tupledict_to_array(model.addVars(range(n),[0],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="c"))
#    g=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
#    u={}
#    epsilon={}
#    mu={}
#    bigM=100
#    eps_tol=0.001
#    for poly in list_of_polytopes:
#        for sign in ["+","-"]:
#            u[poly,sign]=tupledict_to_array(model.addVars(range(poly.P.h.shape[0]),[0],lb=0,ub=GRB.INFINITY,name="u"))
#            epsilon[poly,sign]=model.addVar(lb=-bigM/2,ub=bigM/2,name="eps")
#            mu[poly,sign]=model.addVar(vtype=GRB.BINARY,name="mu(-)")
#    mu["sum","+"]=model.addVar(lb=0,name="sum+")
#    mu["sum","-"]=model.addVar(lb=0,name="sum-")
#    model.update()
#    for poly in list_of_polytopes:
#        constraints_list_of_tuples(model,[(poly.P.h.T,u[poly,"+"]),
#                                          (np.eye(1),np.array([g]).reshape(1,1)),
#                                          (np.eye(1),np.array([epsilon[poly,"+"]]).reshape(1,1)),
#                                          (-poly.t.T,c)],sign="=")
#        constraints_list_of_tuples(model,[(poly.P.h.T,u[poly,"-"]),
#                                          (-np.eye(1),np.array([g]).reshape(1,1)),
#                                          (np.eye(1),np.array([epsilon[poly,"-"]]).reshape(1,1)),
#                                          (poly.t.T,c)],sign="=")    
#        constraints_list_of_tuples(model,[(poly.P.H.T, u[poly,"+"]),(poly.T.T,c)],sign="=")     
#        constraints_list_of_tuples(model,[(poly.P.H.T, u[poly,"-"]),(-poly.T.T,c)],sign="=")     
#        # bigM_formulation
#        for sign in ["+","-"]:
#            model.addConstr(bigM*mu[poly,sign]>=epsilon[poly,sign]-eps_tol)
#            model.addConstr(bigM-bigM*mu[poly,sign]>=-epsilon[poly,sign]+eps_tol)
#    # sum of mu values
#    model.addConstr(mu["sum","+"]==LinExpr([(1,mu[poly,"+"]) for poly in list_of_polytopes]))
#    model.addConstr(mu["sum","-"]==LinExpr([(1,mu[poly,"-"]) for poly in list_of_polytopes]))
#    # The first value of c is one
#    model.addConstr(c[0,0]==1)
#    # Cost objective
#    J=QuadExpr()
#    N=len(list_of_polytopes)+0.0
#    J.add(mu["sum","+"]*mu["sum","+"]+mu["sum","-"]*mu["sum","-"]-N*mu["sum","+"]-N*mu["sum","-"]+N**2/2.0-mu["sum","+"]-mu["sum","-"])
#    J=J*100/(N**2)
#    model.setObjective(J,GRB.MINIMIZE)
#    model.optimize()
##    for poly in list_of_polytopes:
##        print poly
##        print poly.t.T,"mu+:",mu[poly,"+"].X,"mu-:",mu[poly,"-"].X,"eps+:",epsilon[poly,"+"].X,"eps-:",epsilon[poly,"-"].X
##        print "u+=",valuation(u[poly,"+"]).T,"u-=",valuation(u[poly,"-"]).T
##        print "+",np.dot(valuation(u[poly,"+"]).T,poly.P.h)+np.dot(valuation(c).T,poly.t)-g.X#-epsilon[poly,"+"].X
##        print "-",np.dot(valuation(u[poly,"-"]).T,poly.P.h)-np.dot(valuation(c).T,poly.t)+g.X#-epsilon[poly,"+"].X
#    mu_n_positive={poly:np.round(mu[poly,"+"].X) for poly in list_of_polytopes}
#    mu_n_negative={poly:np.round(mu[poly,"-"].X) for poly in list_of_polytopes}
#    print "sums",mu["sum","-"].X,mu["sum","+"].X
#    return valuation(c),g.X,mu_n_positive,mu_n_negative

    
        
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