#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 16:06:32 2019

@author: sadra
"""
import numpy as np
import scipy as sp
from gurobipy import Model,GRB,LinExpr,QuadExpr,tuplelist

from itertools import product as pr

from pypolycontain.lib.polytope import polytope
from pypolycontain.lib.AH_polytope import AH_polytope,is_nonempty
from PWA_lib.system.linear_cells import linear_cell as LC

class system_hard_contact_PWA_numeric:
    def __init__(self,symbolic_system):
        self.A={}
        self.B_u={}
        self.B_lambda={}
        self.c={}
        
        self.Eta=[]
        self.E={}
        self.e={}
        
        self.symbolic_system=symbolic_system
        self.list_of_contact_points=symbolic_system.list_of_contact_points
        self.n,self.m_u,self.m_lambda=symbolic_system.n,symbolic_system.m_u,symbolic_system.m_lambda
        self.Sigma=list(pr(*[i. Sigma for i in self.list_of_contact_points]))


       
    def add_environment(self,Eta,epsilon_max=[0],epsilon_min=[0]):
        self.Eta.append(Eta)
        self.E[Eta],self.e[Eta]={},{}
        H,h={},{}
        dynamical_matrices=self.symbolic_system._evaluate_dynamical_matrices(Eta.dict_of_values)
        for contact_point in self.list_of_contact_points:
            H[contact_point],h[contact_point]=contact_point.Evaluate_polyhedral_matrices(dynamical_matrices,Eta.dict_of_values,epsilon_max,epsilon_min)
            self.E[Eta][contact_point],self.e[Eta][contact_point]={},{}
            for sigma in contact_point.Sigma:
                self.E[Eta][contact_point][sigma]=H[contact_point][sigma]
                self.e[Eta][contact_point][sigma]=h[contact_point][sigma]
        self.A[Eta]=dynamical_matrices.A
        self.B_u[Eta]=dynamical_matrices.B_u
        self.B_lambda[Eta]=dynamical_matrices.B_lambda
        self.c[Eta]=dynamical_matrices.c
        
class environment:
    def __init__(self,name=0):
        self.dict_of_values={}
        self.name="Env "+str(name)
        
    def __repr__(self):
        return self.name
    
    
def merge_timed_vectors_glue(list_of_timed_dicts):
    X={}
    t0=0
    for D in list_of_timed_dicts:
        T=max(D.keys())
        for t in range(T+1):
            X[t0+t]=D[t]
        t0+=T
    return X

def merge_timed_vectors_skip(list_of_timed_dicts):
    X={}
    t0=0
    for D in list_of_timed_dicts:
        T=max(D.keys())
        for t in range(T+1):
            X[t0+t]=D[t]
        t0+=T+1
    return X


def enumerate_modes_E(sys,Eta):
    """
    System sys is numeric
    """
    E,e=sys.E[Eta],sys.e[Eta]
    E_new,e_new={},{}
    for sys_sigma in sys.Sigma:
        E_new[sys_sigma]=np.vstack(([E[sys.list_of_contact_points[j]][sys_sigma[j]] for j in range(len(sys.list_of_contact_points))]))
        e_new[sys_sigma]=np.vstack(([e[sys.list_of_contact_points[j]][sys_sigma[j]] for j in range(len(sys.list_of_contact_points))]))
    return E_new ,e_new
    
    
def simulate_system_numeric(sys,x,u,Eta):
    """
    System sys is numeric
    """
    E,e=enumerate_modes_E(sys,Eta)
    for sys_sigma in sys.Sigma:
        e[sys_sigma]-np.dot(E[sys_sigma][:,0:sys.n+sys.m_u],np.vstack((x,u)))
     
    raise NotImplementedError
    
def trajectory_to_list_of_linear_cells(sys,Eta,x_traj,u_traj,lambda_traj,mode_traj):
    """
    Assumption: single environment
    """
    E,e=enumerate_modes_E(sys,Eta)
    T=max(mode_traj.keys())+1
    list_of_cells=[]
    if Eta not in sys.Eta:
        raise ValueError("Enviornement not defined / added to the system yet")
    for t in range(T):
        p=polytope(E[mode_traj[t]],e[mode_traj[t]])
        new_cell=LC(sys.A[Eta],np.hstack((sys.B_u[Eta],sys.B_lambda[Eta])),sys.c[Eta],p)
        list_of_cells.append(new_cell)
        # Verify correctness:
        v=np.hstack((x_traj[t],u_traj[t],lambda_traj[t])).reshape(p.n,1)
        if p.if_inside(v,tol=10**-5)==False:
            raise ValueError("Nominal Trajectory not inside polytope! Maybe an issue of numerical tolerance")
    return list_of_cells 

def environment_from_state(symbolic_sys,x,h):
    """
    Takes a symbolic system, state, control, contact forces and time-step, build an environment
    """
    Eta_0=environment(0)
    Eta_0.dict_of_values={}
    for i in range(symbolic_sys.n):
        Eta_0.dict_of_values[symbolic_sys.x[i]]=x[i]
    for u in symbolic_sys.u:
        Eta_0.dict_of_values[u]=0
    for u_lambda in symbolic_sys.u_lambda:
        Eta_0.dict_of_values[u_lambda]=0
    Eta_0.dict_of_values[symbolic_sys.h]=h
    return Eta_0

def PWA_cells_from_state(symbolic_sys,x,h,epsilon_min,epsilon_max):
    sys=system_hard_contact_PWA_numeric(symbolic_sys)
    Eta_now=environment_from_state(symbolic_sys,x,h)
    sys.add_environment(Eta_now,epsilon_max,epsilon_min)
    E,e=enumerate_modes_E(sys,Eta_now)
    list_of_linear_cells=[]
    for mode in sys.Sigma:
        p=polytope(E[mode],e[mode])
        if is_nonempty(p):
            new_cell=LC(sys.A[Eta_now],np.hstack((sys.B_u[Eta_now],sys.B_lambda[Eta_now])),sys.c[Eta_now],p)
            new_cell.name=mode
            list_of_linear_cells.append(new_cell)
    return list_of_linear_cells

def hybrid_reachable_sets_from_state(symbolic_sys,x,h,epsilon_min,epsilon_max):
    sys=system_hard_contact_PWA_numeric(symbolic_sys)
    Eta_now=environment_from_state(symbolic_sys,x,h)
    sys.add_environment(Eta_now,epsilon_max,epsilon_min)
    E,e=enumerate_modes_E(sys,Eta_now)
    list_of_sets=[]
    for mode in sys.Sigma:
        H=E[mode][:,symbolic_sys.n:]
        q=H.shape[0]
        h=e[mode].reshape(q,1)-np.dot(E[mode][:,:symbolic_sys.n],x.reshape(symbolic_sys.n,1))
        p=polytope(H,h)
        if is_nonempty(p):
            B=np.hstack((sys.B_u[Eta_now],sys.B_lambda[Eta_now]))
            t=np.dot(sys.A[Eta_now],x)+sys.c[Eta_now]
            new_set=AH_polytope(T=B,t=t,P=p)
            new_set.name=mode
            list_of_sets.append(new_set)
    return list_of_sets

def trajectory_to_list_of_linear_cells_full_linearization(symbolic_sys,x_traj,u_traj,lambda_traj,mode_traj,h,epsilon_min,epsilon_max):
    """
    Linearize all points
    """
    list_of_cells=[]
    T=max(mode_traj.keys())+1
    for t in range(T):
        sys=system_hard_contact_PWA_numeric(symbolic_sys)
        Eta_now=environment_from_state(symbolic_sys,x_traj[t],h)
        sys.add_environment(Eta_now,epsilon_max,epsilon_min)
        E,e=enumerate_modes_E(sys,Eta_now)
        p=polytope(E[mode_traj[t]],e[mode_traj[t]])
        new_cell=LC(sys.A[Eta_now],np.hstack((sys.B_u[Eta_now],sys.B_lambda[Eta_now])),sys.c[Eta_now],p)
        list_of_cells.append(new_cell)
        # Verify correctness:
        v=np.hstack((x_traj[t],u_traj[t],lambda_traj[t])).reshape(p.n,1)
#        print p.h-np.dot(p.H,v)
#        if p.if_inside(v,tol=10**-1)==False:
#            raise ValueError("Nominal Trajectory not inside polytope! Maybe an issue of numerical tolerance")
    return list_of_cells