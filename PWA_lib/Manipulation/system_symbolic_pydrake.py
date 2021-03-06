#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 18:52:43 2019

@author: sadra
"""

import pydrake.symbolic as sym

import numpy as np

class system_symbolic:
    def __init__(self,name,d):
        self.name=name+" (%dD)"%d
        self.d=d # Dimensions
        self.list_of_contact_points=[] # list of contact points
        # Variables
        self.q_o=np.empty(0,dtype="object") # Position vector of the Object(s)
        self.q_m=np.empty(0,dtype="object") # Position vector of the Manipulator(s): those controlled by velocity commands
        self.u_torques=np.zeros(0) # Vector of external torque/forces
        # Functions
        self.M=np.empty(0,dtype="object") # Mass matrix
        self.C=np.empty(0,dtype="object") # C matrix
        self.tau_g=np.empty(0,dtype="object") # tau_g vector
        self.B=np.empty(0,dtype="object") # B matrix
        # discrete_time step
        self.h=sym.Variable(name="h")
        # Contacts:
#        self.u_lambda=np.hstack(*[contact_point.u_lambda for contact_point in self.contact_points])
        # Dynamics
        self.dynamics=dynamics(None,None,None,None)
                
    def __repr__(self):
        return self.name+" with %d states, %d controls, and %d contact points"\
            %(len(self.x),len(self.u),len(self.list_of_contact_points))
            
    def _build_state_space(self):
        """ Build the system's state and control space"""
#        self.qdot_o=np.array([sym.Variable("qdot%d"%i) for i in range(len(self.q_o))])
        self.v_o=np.array([sym.Variable(str(a.__repr__())[10:str(a.__repr__()).index(',')-1]+"_dot") for a in self.q_o])
        self.u_m=np.array([sym.Variable("u_"+str(a.__repr__())[10:str(a.__repr__()).index(',')-1]) for a in self.q_m])
#        self.u_m=np.array([sym.Variable("u_m%d"%i) for i in range(len(self.q_m))])
        if self.d==2:
            self.u_lambda=np.array([[sym.Variable("lambda_%s_n"%c.name),sym.Variable("lambda_%s_t"%c.name)] \
                                     for c in self.list_of_contact_points]).reshape(2*len(self.list_of_contact_points))
        else:
            raise NotImplementedError
        self.q=np.hstack((self.q_o,self.q_m))
        self.x=np.hstack((self.q_o,self.q_m,self.v_o))
#        self.u=np.hstack((self.u_torques,self.u_m,self.u_lambda))
        self.u=np.hstack((self.u_torques,self.u_m))
        # self.tau_c
        self.tau_c=np.dot(self.C,self.v_o)
        # The Jacobian
        self.J=np.hstack(([c.J for c in self.list_of_contact_points]))

    
    def _build_dimensions(self):
        # Dimensions
        self.n,self.m_u,self.m_lambda=len(self.x),len(self.u),len(self.u_lambda)
        
        
    def _set_mass_matrix_inverse(self,M_inv): 
        # Inverse of Mass Matrix is needed
        self.M_inv=M_inv
        
    def _polyhedrize_constraints_symbolic(self,zeta):
        """
        Given vector of symbolic constraint set 
        zeta(q,v_o,u_torques,u_m,lambdas) <= 0
        Compute the polyhedral set H(x,u,lambda) <= h
        """
        try:
            H=np.hstack((sym.Jacobian(zeta,self.q),sym.Jacobian(zeta,self.v_o),
                         sym.Jacobian(zeta,self.u_torques),sym.Jacobian(zeta,self.u_m),
                         sym.Jacobian(zeta,self.u_lambda) ))
        except:
            H=np.hstack((sym.Jacobian(zeta,self.q),sym.Jacobian(zeta,self.v_o),
                         sym.Jacobian(zeta,self.u_m), # There is no torque input
                         sym.Jacobian(zeta,self.u_lambda) ))            
        h=np.dot(H,np.hstack((self.x,self.u,self.u_lambda)))-zeta
        return (H,h)
    
    def _time_derivative(self,xi):
        """
        Time derivative of the vector of symbolic expression xi(q)
        """
        return np.dot(sym.Jacobian(xi,self.q_o),self.v_o)+np.dot(sym.Jacobian(xi,self.q_m),self.u_m)
            
    def _linearize_dynamics_symbolic(self):
        # First compute A matrix
        A_11=np.eye(self.q.shape[0])
        A_12=np.vstack((-self.h*np.eye(self.q_o.shape[0]),np.zeros((self.q_m.shape[0],self.v_o.shape[0]))))
        alpha_1=sym.Jacobian(self.tau_c,self.q)-sym.Jacobian(self.tau_g,self.q)-\
                sum([sym.Jacobian(self.B[:,i]*self.u_torques[i],self.q) for i in range(self.B.shape[1])])-\
                sum([sym.Jacobian(self.J[:,i]*self.u_lambda[i],self.q) for i in range(self.J.shape[1])])
        A_21=self.h*np.dot(self.M_inv,alpha_1)
        alpha_2=sym.Jacobian(self.tau_c,self.v_o)-sum([ 
                sym.Jacobian(self.B[:,i]*self.u_torques[i],self.v_o) for i in range(self.B.shape[1])])
        A_22=np.eye(self.v_o.shape[0])-self.h*np.dot(self.M_inv,alpha_2)
        self.dynamics.A_inv=np.vstack ( ( np.hstack((A_11,A_12))  ,  np.hstack((A_21,A_22)) ) )

        
        self.dynamics.B_u_1=np.vstack((np.zeros((self.q.shape[0],self.B.shape[1])),self.h*self.B))
        self.dynamics.B_u_2=np.vstack((np.zeros((self.q_o.shape[0],self.q_m.shape[0])),\
                                       self.h*np.eye(self.q_m.shape[0]),\
                                       np.zeros((self.v_o.shape[0],self.q_m.shape[0]))))
       
        self.dynamics.B_u=np.hstack((self.dynamics.B_u_1,self.dynamics.B_u_2))
        
        self.dynamics.B_lambda=np.vstack((np.zeros((self.q.shape[0],self.J.shape[1])),self.h*self.J))

        alpha_3=-self.tau_c+np.dot(sym.Jacobian(self.tau_c,self.q),self.q)+\
            np.dot(sym.Jacobian(self.tau_c,self.v_o),self.v_o)+\
            self.tau_g-np.dot(sym.Jacobian(self.tau_g,self.q),self.q)+\
            np.dot(self.B,self.u_torques)-\
            sum([np.dot(sym.Jacobian(self.B[:,i],self.q),self.q)*self.u_torques[i] for i in range(self.B.shape[1])])-\
            sum([np.dot(sym.Jacobian(self.B[:,i],self.v_o),self.v_o)*self.u_torques[i] for i in range(self.B.shape[1])])-\
            sum([np.dot(sym.Jacobian(self.J[:,i],self.q),self.q)*self.u_lambda[i] for i in range(self.J.shape[1])])
        
        self.dynamics.c=np.hstack((np.zeros(self.q.shape[0]),self.h*alpha_3))
    
    def _evaluate_dynamical_matrices(self,Eta):
        D=self.dynamics
        A_n=np.linalg.inv(sym.Evaluate(D.A_inv,Eta))
        B_u_n=np.dot(A_n,sym.Evaluate(D.B_u,Eta))
        B_lambda_n=np.dot(A_n,sym.Evaluate(D.B_lambda,Eta))
        c_n=np.dot(A_n,sym.Evaluate(D.c,Eta))
        return dynamics(A_n,B_u_n,B_lambda_n,c_n)
    
    def build_and_linearize(self):
        self._build_state_space()
        self._linearize_dynamics_symbolic()
        for C in self.list_of_contact_points:
            C._contact_geometery_all()
            C._contact_forces_all()
        self._build_dimensions()

            
class dynamics:
    def __init__(self,A,B_u,B_lambda,c):
        self.A=A
        self.B_u=B_u
        self.B_lambda=B_lambda
        self.c=c