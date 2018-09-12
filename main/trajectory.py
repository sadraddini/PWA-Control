"""
Created on Mon June 25 2018

@author:    Sadra Sadraddini
            CSAIL, MIT 
            32 Vassar st, Cambridge, MA 02139
            s a d r a @ mit.edu
"""


# Primary imports
import numpy as np
from gurobipy import Model,GRB,LinExpr,QuadExpr
from random import choice as rchoice
from random import random

from time import time

import sys
sys.path.append('..')

# Secondary imports
from main.auxilary_methods import find_mode,valuation,mode_sequence
from main.ana_system import state,cost_state

from Verification.verification_polytopes import re_verification

def polytope_trajectory(s,x0,state_end,T,alpha_start,eps=0.1,coin=random()):
    (model,x,u,G,theta,z)=s.library[T]
    n_vars=len(model.getVars())
    n_constraints=len(model.getConstrs())
    new_var_count=0
    new_constraint_count=0
    J_area=LinExpr()
    d_min=model.addVar(lb=0.0001,name="new var %d"%new_var_count)
    new_var_count+=1
    beta=10**2 # Weight of infinity norm
    model.update()
    print("coin=",coin)
    for row in range(s.n):
        for column in range(s.n):
            if coin<0.1:
                if row<column:
                    model.addConstr(G[0][row,column]==0,name="constraint %d"%new_constraint_count)
                    new_constraint_count+=1
            elif coin>0.9:
                if row>column:
                    model.addConstr(G[0][row,column]==0,name="constraint %d"%new_constraint_count)
                    new_constraint_count+=1                
            if row==column:
                model.addConstr(G[0][row,column]>=d_min/s.weight[row],name="constraint %d"%new_constraint_count)
                new_constraint_count+=1
    J_area.add(-d_min*T*s.n*beta)
    for row in range(s.n):
        for t in range(T):
            J_area.add(-G[t][row,row]*s.weight[row])
    # Terminal Constraint
    terminal_constraint(s,x,G,T,model,state_end)
    # Starting Point
    i_start=find_mode(s,x0)
    for i in s.modes:
        model.addConstr(z[0,i]==int(i==i_start))
    # model.setParam('OutputFlag',False)
    if alpha_start==-1:
        x_delta={}
        for row in range(s.n):
            x_delta[row]=model.addVar(lb=-eps/s.weight[row],ub=eps/s.weight[row])
        model.update()
        for row in range(s.n):
            model.addConstr(x[0][row,0]==x0[row,0]+x_delta[row])
        model.setObjective(J_area)
        model.optimize()
    else:          
        model.setObjective(J_area)
        model.optimize()
    if model.Status!=2 and model.Status!=11:
        flag=False
        print("*"*20,"False flag",model.Status)
        final=(x,u,G,theta,z,flag)
    else:
        flag=True
        x_n=valuation(x)
        u_n=valuation(u)
        G_n=valuation(G)
        theta_n=valuation(theta)
        z_n=mode_sequence(s,z)
        if abs(np.linalg.det(G_n[0]))<10**-15:
            flag=False
        final=(x_n,u_n,G_n,theta_n,z_n,flag)
    print("starting removal process")
    print("start=",n_vars,"variables and ",n_constraints," constraints")
    new_n_vars=len(model.getVars())
    new_n_constraints=len(model.getConstrs())
    print("new:",new_n_vars,"variables and ",new_n_constraints," constraints")    
    time_start=time()
    for i in range(new_n_vars-n_vars):
        model.remove(model.getVars()[-i-1])
    model.update()
    for i in range(new_n_constraints-n_constraints):
        model.remove(model.getConstrs()[-i-1]) 
    model.update()
    n_vars=len(model.getVars())
    n_constraints=len(model.getConstrs())
    print("final:",n_vars,"variables and ",n_constraints," constraints")    
    print("end of removal in ",time()-time_start," seconds")
    return final
    
def make_state_trajectory_state_end(s,x,u,seq,G,theta,T,state_end):
    flag=re_verification(s,state(x[T-1],G[T-1],seq[T-1],-1,-1,-1),state_end)
    if flag==None:
        return
    else:
        u[T-1]=flag[0]
        theta[T-1]=flag[1]
        x_temp={}
        ID=s.ID
        s.ID+=1
        for t in range(0,T):
            x_temp[t]=state(x[t],G[t],seq[t],ID,t,1)
            x_temp[t].time_to_go=T-t+state_end.time_to_go
        # Build Transitons
        for t in range(0,T-1):
            x_temp[t].successor=(x_temp[t+1],u[t],theta[t])
        x_temp[T-1].successor=(state_end,u[T-1],theta[T-1])
        for t in range(1,T):
            x_temp[t].parent.append(x_temp[t-1])
        state_end.parent.append(x_temp[T-1])
        s.X.extend(list(x_temp.values())[::-1])
    

def subset_MILP(model,G,Pi,H,h,x,z_time_mode):
    """
    Description: Add Farkas lemma constraints for subset inclusion of x+GP in {e|H.<h}
    Inputs: 
        model: Gurobi optimization model
        G: n * n_g generator matrix
        Pi:  n_pi*n matrix where {x|Pi x< 1} is the primitive polytope
        {x| Hx<h} is the set constraint
        x is the point
    Output:
        no direct output. Adds constraints to the model. 
        FUTURE: we may add lambda positive matrix as an output for debugging
    """
    bigM=10
    (n,n_g)=G.shape
    (n_p,n)=Pi.shape
    (n_h,n)=H.shape
    Lambda=np.empty((n_h,n_p),dtype='object')
    for row in range(n_h):
        for column in range(n_p):
            Lambda[row,column]=model.addVar(lb=0)
    model.update()
    # Lambda * Pi = H * G
    for row in range(n_h):
        for column in range(n_g):
            s_left=LinExpr()
            s_right=LinExpr()
            for k in range(n_p):
                s_left.add(Lambda[row,k]*Pi[k,column])
            for k in range(n):
                s_right.add(H[row,k]*G[k,column])
            model.addConstr(s_left==s_right)
    # Lambda * 1 <= H*
    for row in range(n_h):
        s_left=LinExpr()
        s_right=LinExpr()
        for k in range(n_p):
            s_left.add(Lambda[row,k])
        for k in range(n):
            s_right.add(H[row,k]*x[k,0])
        model.addConstr(s_left<=h[row,0]-s_right+bigM-bigM*z_time_mode) 

def subset(model,G,Pi,H,h,x):
    """
    Description: Add Farkas lemma constraints for subset inclusion of x+GP in {e|H.<h}
    Inputs: 
        model: Gurobi optimization model
        G: n * n_g generator matrix
        Pi:  n_pi*n matrix where {x|Pi x< 1} is the primitive polytope
        {x| Hx<h} is the set constraint
        x is the point
    Output:
        no direct output. Adds constraints to the model. 
        FUTURE: we may add lambda positive matrix as an output for debugging
    """
    (n,n_g)=G.shape
    (n_p,n)=Pi.shape
    (n_h,n)=H.shape
    Lambda=np.empty((n_h,n_p),dtype='object')
    for row in range(n_h):
        for column in range(n_p):
            Lambda[row,column]=model.addVar(lb=0)
    model.update()
    # Lambda * Pi = H * G
    for row in range(n_h):
        for column in range(n_g):
            s_left=LinExpr()
            s_right=LinExpr()
            for k in range(n_p):
                s_left.add(Lambda[row,k]*Pi[k,column])
            for k in range(n):
                s_right.add(H[row,k]*G[k,column])
            model.addConstr(s_left==s_right)
    # Lambda * 1 <= H*
    for row in range(n_h):
        s_left=LinExpr()
        s_right=LinExpr()
        for k in range(n_p):
            s_left.add(Lambda[row,k])
        for k in range(n):
            s_right.add(H[row,k]*x[k,0])
        model.addConstr(s_left<=h[row,0]-s_right) 

def terminal_constraint(s,x,G,T,model,state):
    """
        Terminal Constraint 
    """
    subset(model,G[T],s.Pi,state.polytope.H,state.polytope.h,x[T])
            
def terminal_constraint_old(s,x,G,T,model,state):
        # Terminal Constraint
    if state.character!=-1 and state.volume_flag==True:
        print("taking G approach")
        H=np.dot(s.Pi,state.Ginv)
        h=np.ones((s.Pi.shape[0],1))+np.dot(H,state.x)
        #print("Nonzero Volume!","Ginv=",state2.Ginv,"H=",H,"h=",h)
        subset(model,G[T],s.Pi,H,h,x[T])
    else:
        print("taking Lambda approach")
        # The vertices by G[i,T]
        Lambda={}
        for alpha in range(2**s.n):
            for beta in range(2**state.G.shape[1]):
                Lambda[alpha,beta]=model.addVar(lb=0)
        model.update()
        for alpha in range(2**s.n):
            for row in range(s.n):
                exp_left=LinExpr()
                exp_right=LinExpr()
                for k in range(s.n):
                    exp_left.add(G[T][row,k]*s.vertices[alpha,:].reshape(s.n,1)[k,0])
                for beta in range(2**state.G.shape[1]):
                    exp_right.add(Lambda[alpha,beta]*state.vertices[beta,:].reshape(s.n,1)[row,0])
                model.addConstr(x[T][row,0]+exp_left==exp_right+state.x[row,0])
            lambda_sum=LinExpr()
            for beta in range(2**state.G.shape[1]):
                lambda_sum.add(Lambda[alpha,beta])
            model.addConstr(lambda_sum<=1)
            
def trajectory_cost(s,x,u,seq,G,theta,T):
    J=0
    ID=1001
    state_considered={}
    for t in range(T):
        state_considered[t]=state(x[t],G[t],seq[t],ID,t,7)
    for t in range(T-1):
        state_considered[t].successor=(state_considered[t+1],u[t],theta[t])
        J+=cost_state(s,state_considered[t],s.Q,s.R,s.g)
    return J

def state_trajectory(s,x0,state_end,T):
    model=Model("trajectory of polytopes")
    x={}
    u={}
    theta={}
    z={}
    G_bound=100
    # Mode 1:
    for t in range(T):
        x[t]=np.empty((s.n,1),dtype='object') # n*1
        u[t]=np.empty((s.m,1),dtype='object') # m*1
        theta[t]=np.empty((s.m,s.n),dtype='object') # n*m
        for row in range(s.n):
            x[t][row,0]=model.addVar(lb=-G_bound,ub=G_bound)
        for row in range(s.m):
            u[t][row,0]=model.addVar(lb=-G_bound,ub=G_bound)
        for row in range(s.m):
            for column in range(s.n):
                theta[t][row,column]=0   
    for t in range(T+1):
        for i in s.modes:
            z[t,i]=model.addVar(vtype=GRB.BINARY)
    x[T]=np.empty((s.n,1),dtype='object') # Final state in Mode i
    for row in range(s.n):
        x[T][row,0]=model.addVar(lb=-G_bound,ub=G_bound)
    G={}
    for t in range(T+1):
        G[t]=np.empty((s.n,s.n),dtype='object')
        for row in range(s.n):
            for column in range(s.n):
                G[t][row,column]=0
    model.update()
    # Trajectory Constraints:
    # Mode i:
    bigM=100
    for i in s.modes:
        for t in range(T):
            for row in range(s.n):
                Ax=LinExpr()
                for k in range(s.n):
                    Ax.add(s.A[i][row,k]*x[t][k,0])
                for k in range(s.m):
                    Ax.add(s.B[i][row,k]*u[t][k,0])
                model.addConstr(x[t+1][row,0]<=Ax+s.c[i][row]+bigM-bigM*z[t,i])
                model.addConstr(x[t+1][row,0]>=Ax+s.c[i][row]-bigM+bigM*z[t,i])
    # Generator Dynamics Constraints:
    for i in s.modes:
        for t in range(T):
            for row in range(s.n):
                for column in range(s.n):
                    AG=LinExpr()
                    for k in range(s.n):
                        AG.add(s.A[i][row,k]*G[t][k,column])
                    for k in range(s.m):
                        AG.add(s.B[i][row,k]*theta[t][k,column])
                    model.addConstr(G[t+1][row,column]<=AG+bigM-bigM*z[t,i])
                    model.addConstr(G[t+1][row,column]>=AG-bigM+bigM*z[t,i])
    # Constraints of modes:
    for t in range(T+1):
        sum_z=LinExpr()
        for i in s.modes:
            sum_z.add(z[t,i])
        model.addConstr(sum_z==1)
    # Constraints of mode subsets
    for t in range(T):
        for i in s.modes:
            subset_MILP(model,G[t],s.Pi,s.H[i],s.h[i],x[t],z[t,i])
            subset_MILP(model,theta[t],s.Pi,s.F[i],s.f[i],u[t],z[t,i])   
    # set objective
    model.update()
    # Terminal Constraint
    terminal_constraint(s,x,G,T,model,state_end)
    # Starting Point
    i_start=find_mode(s,x0)
    for i in s.modes:
        model.addConstr(z[0,i]==int(i==i_start))
    model.setParam('OutputFlag',False)
    for row in range(s.n):
        model.addConstr(x[0][row,0]==x0[row,0])
    model.optimize()
    if model.Status!=2 and model.Status!=11:
        flag=False
#        print("*"*20,"Flag is False and Status is",model.Status)
        return (x,u,z,flag)
    else:
        flag=True
#        print("*"*20,"Flag is True and Status is",model.Status)
        x_n=valuation(x)
        u_n=valuation(u)
        z_n=mode_sequence(s,z)
        return (x_n,u_n,z_n,flag)