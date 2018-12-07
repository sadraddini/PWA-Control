# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 18:11:28 2018

@author: sadra
"""

from gurobipy import Model,GRB,LinExpr
import time as time

N=100
T=4

model=Model("1")
t=time.time()
x=model.addVars(range(N),range(T),lb=-2,ub=2,obj=-1)
model.update()
model.addConstrs(x.sum(j,"*")==0 for j in range(N))
print "1:",time.time()-t
model.optimize()

#model=Model("2")
#t2=time.time()
#x={}
#for n in range(N):
#    for t in range(T):
#        x=model.addVars(range(N),range(T),lb=-2,ub=2,obj=-1)
#model.update()
#print "2:",time.time()-t2
