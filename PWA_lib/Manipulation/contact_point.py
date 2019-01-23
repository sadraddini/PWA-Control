#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 15:50:54 2019

@author: sadra
"""
from sympy import Symbol,pi,sin,cos,Function,diff,Matrix,symbols,lambdify
import numpy as np


"""
Important Assumptions that are hard coded here: 
     * Everything here is restricted to two dimensions
     * The integration is Eulerian f(x+h)=f(x)+hf'(x)
"""
        
        
class contact_point:
    """
        Each contact point:
            * Forces as control inputs which should satisfy Friction and Complementarity Constraints
            * 2D so far: 4 Possibilites: No contact, slide left, slide right, sticking
    """
    
    """ SYMBOLIC:"""
    def __init__(self,sys,index,phi,psi,J,friction=0.5,K=100,damping=1,name="Contact Point ",contact_model="soft"):
        self.sys=sys # What system it belongs to
        self.phi=phi # Penetration position, symbolic expression
        self.psi=psi # Sliding position, symbolic expression
        self.J=J # Jacobian for the forces
        self.polytope={}
        self.index=index # Must be an integer between [0,N), N the number of contact points
        self.friction=friction # Friction
        self.contact_model=contact_model # Damper coefficient, used for soft contact
        # Only used for soft contact
        self.damping=damping # Damper coefficient, used for soft contact
        self.K=K # Spring coefficient, used for soft contact
        
        self.sys.contact_points.append(self)
        self.name=name+str(len(self.sys.contact_points))+" contact model: "+self.contact_model
        
    def __repr__(self):
        return self.name
        
    def get_determiners_symbolic(self):
        """
        The used inputs of the function are two: 
            * phi: penetration function
            * psi: sliding coordinates function
        The output of this argument is symbolic_determiner, or sigma
            * phi: phi itself 
            * psi: psi_itself
            * v_phi: phi_dot (only needed for particular contact modes)
            * v_psi: psi_dot 
            * phi_x:  derivative of phi with respect to x: vector
            * phi_u: partial derivative of phi with respect to u: vector
            * v_psi_x: derivative of v_psi with respect to x: vector
            * v_psi_u: derivative of v_psi with respect to u: vector
            What we have is a vector of 2+2+2*(n+m) vector
        """
        n=len(self.sys.x)
        m=len(self.sys.u)
        _output=[0]*(4+2*(n+m))
        t=self.sys.t
        v_phi=diff(self.phi,t) 
        v_psi=diff(self.psi,t)
        _output[0:4]=self.phi,self.psi,v_phi,v_psi
        # v_psi=self.sys.sub_derivatives(v_psi)
        phi_x=[diff(self.phi,var) for var in self.sys.x]
        phi_u=[diff(self.phi,var) for var in self.sys.u]
        v_psi_x=[diff(v_psi,var) for var in self.sys.x]
        v_psi_u=[diff(v_psi,var) for var in self.sys.u]
        _output[4:4+n]=phi_x[0:n]
        _output[4+n:4+n+m]=phi_u[0:m]
        _output[4+n+m:4+n+m+n]=v_psi_x[0:m]
        _output[4+n+m+n:4+n+m+n+m]=v_psi_u[0:m]
        return (self.sys.x,self.sys.u,self.phi,self.psi,v_phi,v_psi,np.array(phi_x),np.array(phi_u),np.array(v_psi_x),np.array(v_psi_u))

    def get_determiners_symbolic_J(self):
        """
        Given the arguments:
            contact point c: object (see above)
            mode in [no_contact,sticking,sliding_right,sliding_left]
            x_0 and u_0 are nominal points IN THE SAME MODE
        """
        J_x=[diff(self.J,var) for var in self.sys.x]
        J_u=[diff(self.J,var) for var in self.sys.u]
        return self.J,J_x,J_u 

    """ 
    *************************************************
    NUMERICAL: The arguments are numbers, not symbols
    *************************************************
    """  
    def forces_no_contact(self,determiners):
        """
        The constraints are as follows:
            * phi >= 0 --> (phi_x,phi_u)(x,u)+ phi_0 - phi_x*x0 - phi_u*u0 >= phi_epsilon
            * f_n = 0 
            * f_t = 0
        """
        assert len(determiners)==10 # This is hard coded. Wonder why I keep psi as one of them:D
        x,u,phi,psi,v_phi,v_psi,phi_x,phi_u,v_psi_x,v_psi_u=determiners
        n=len(self.sys.x)
        m=len(self.sys.u)
        nC=len(self.sys.contact_points)
        # phi>=0
        H1=-np.hstack((phi_x,phi_u))
        h1=phi-np.dot(phi_x,x)-np.dot(phi_u,u)
        # f=0
        H2=np.zeros((4,m+n))
        h2=np.zeros((4,1))
        print H2.shape,H1.shape
        H2[:,n+m-2*nC+2*self.index:n+m-2*nC+2*(self.index+1)]=np.vstack((np.eye(2),-np.eye(2)))
        H=np.vstack((H1,H2))
        h=np.vstack((h1,h2))
        return (H,h)
    
    def forces_sticking(self,determiners,v_epsilon=0.01,maximum_penetration=0.01):
        """
        The constraints are as follows:
            * phi <= 0 -> (phi_x,phi_u)(x,u)+ phi_0 - phi_x*x0 - phi_u*u0 <= 0
            * f_n >= 0 : This is a cosntraint on control
            * |f_T| < mu * f_n : This is a constraint on control
            * |v_psi| <= v_epsilon -> | (v_psi_x,v_psi_u)(x,u) +v_psi_0 - (v_psi_x,v_psi_u)(x_0,u_0) | <= v_epsilon
        """
        assert len(determiners)==10 # This is hard coded. Wonder why I keep psi as one of them:D
        x,u,phi,psi,v_phi,v_psi,phi_x,phi_u,v_psi_x,v_psi_u=determiners
        # Some Logistics
        n=len(self.sys.x)
        m=len(self.sys.u)
        nC=len(self.sys.contact_points)
        # phi <= 0
        H1=np.hstack((phi_x,phi_u))
        h1=-phi+np.dot(phi_x,x)+np.dot(phi_u,u)
        # f_n >= 0 
        H2=np.zeros((1,m+n))
        h2=np.zeros((1,1))
        H2[:,n+m-2*nC+2*self.index]=-1
        # |f_t|<=mu*f_n
        H3=np.zeros((2,m+n))
        h3=np.zeros((2,1))
        H3[:,n+m-2*nC+2*self.index:n+m-2*nC+2*(self.index+1)]=np.array([[-self.friction,-1],[-self.friction,1]])
        # |v_psi|<= v_epsilon
        H4=np.hstack((v_psi_x,v_psi_u))
        h4=v_epsilon-v_psi+np.dot(v_psi_x,x)+np.dot(v_psi_u,u)
        H5=-np.hstack((v_psi_x,v_psi_u))
        h5=v_epsilon+v_psi-np.dot(v_psi_x,x)-np.dot(v_psi_u,u)
        # Contact model
        if self.contact_model=="soft":
            # we have f_n = - K * phi - c * phi_dot --> K*phi+K*(phi_x,phi_u)(x,u)+c*phi_dot+f_n=0 
            H6=np.hstack((phi_x,phi_u))*self.K
            H6[n+m-2*nC+2*self.index]=1
            h6=-self.K*phi-self.damping*v_phi+self.K*(np.dot(phi_x,x)+np.dot(phi_u,u))
            H7=-H6
            h7=self.K*phi+self.damping*v_phi-self.K*(np.dot(phi_x,x)+np.dot(phi_u,u))
            H=np.vstack((H1,H2,H3,H4,H5,H6,H7))
            h=np.vstack((h1,h2,h3,h4,h5,h6,h7)).reshape(8,1)
            return (H,h)
        elif self.contact_model=="hard_implicit":
            # I only impose a maximum penetration constraint
            # phi>=-maximum_penetration: -(phi_x,phi_u)(x,u) + phi_x x+ phi_u u - phi <= maximum_penetration
            H6=-np.hstack((phi_x,phi_u))
            h6=phi-np.dot(phi_x,x)-np.dot(phi_u,u)+maximum_penetration
            H=np.vstack((H1,H2,H3,H4,H5,H6))
            h=np.vstack((h1,h2,h3,h4,h5,h6)).reshape(7,1)
            return (H,h)           
        else:
            raise ValueError("Unknown contact model")
        
        
        
    def force_slide_positive(self,determiners,v_epsilon=0.01,maximum_penetration=0.01):
        """
        The constraints are as as follows:
            * phi <= 0 -> (phi_x,phi_u)(x,u)+ phi_0 - phi_x*x0 - phi_u*u0 <= 0
            * f_n >= 0 : This is a cosntraint on control
            * f_T = mu * f_n : This is a constraint on control - sliding in the direction of v_psi
            * v_psi >=v_epsilon -> (v_psi_x,v_psi_u)(x,u) +v_psi_0 - (v_psi_x,v_psi_u)(x_0,u_0)  >= v_epsilon
        """
        assert len(determiners)==10 # This is hard coded. Wonder why I keep psi as one of them:D
        x,u,phi,psi,v_phi,v_psi,phi_x,phi_u,v_psi_x,v_psi_u=determiners
        # Some Logistics
        n=len(self.sys.x)
        m=len(self.sys.u)
        nC=len(self.sys.contact_points)
        # phi <= 0
        H1=np.hstack((phi_x,phi_u))
        h1=-phi+np.dot(phi_x,x)+np.dot(phi_u,u)
        # f_n >= 0 
        H2=np.zeros((1,m+n))
        h2=np.zeros((1,1))
        H2[:,n+m-2*nC+2*self.index]=-1
        # f_t==mu*f_n
        H3=np.zeros((2,m+n))
        h3=np.zeros((2,1))
        H3[:,n+m-2*nC+2*self.index:n+m-2*nC+2*(self.index+1)]=np.array([[-self.friction,1],[self.friction,-1]])
        # v_psi >= v_epsilon
        H4=-np.hstack((v_psi_x,v_psi_u))
        h4=-v_epsilon+v_psi-np.dot(v_psi_x,x)-np.dot(v_psi_u,u)
        # Contact model
        if self.contact_model=="soft":
            # we have f_n = - K * phi - c * phi_dot --> K*phi+K*(phi_x,phi_u)(x,u)+c*phi_dot+f_n=0 
            H5=np.hstack((phi_x,phi_u))*self.K
            H5[n+m-2*nC+2*self.index]=1
            h5=-self.K*phi-self.damping*v_phi+self.K*(np.dot(phi_x,x)+np.dot(phi_u,u))
            H6=-H5
            h6=self.K*phi+self.damping*v_phi-self.K*(np.dot(phi_x,x)+np.dot(phi_u,u))
            H=np.vstack((H1,H2,H3,H4,H5,H6))
            h=np.vstack((h1,h2,h3,h4,h5,h6)).reshape(7,1)
            return (H,h)
        elif self.contact_model=="hard_implicit":
            # I only impose a maximum penetration constraint
            # phi>=-maximum_penetration: -(phi_x,phi_u)(x,u) + phi_x x+ phi_u u - phi <= maximum_penetration
            H5=-np.hstack((phi_x,phi_u))
            h5=phi-np.dot(phi_x,x)-np.dot(phi_u,u)+maximum_penetration
            H=np.vstack((H1,H2,H3,H4,H5))
            h=np.vstack((h1,h2,h3,h4,h5)).reshape(6,1)
            return (H,h)           
        else:
            raise ValueError("Unknown contact model")

    def force_slide_negative(self,determiners,v_epsilon=0.01,maximum_penetration=0.01):
        """
        The constraints are as as follows:
            * phi <= 0 -> (phi_x,phi_u)(x,u)+ phi_0 - phi_x*x0 - phi_u*u0 <= 0
            * f_n >= 0 : This is a cosntraint on control
            * f_T = mu * f_n : This is a constraint on control - sliding in the direction of v_psi
            * v_psi >=v_epsilon -> (v_psi_x,v_psi_u)(x,u) +v_psi_0 - (v_psi_x,v_psi_u)(x_0,u_0)  >= v_epsilon
        """
        assert len(determiners)==10 # This is hard coded. Wonder why I keep psi as one of them:D
        x,u,phi,psi,v_phi,v_psi,phi_x,phi_u,v_psi_x,v_psi_u=determiners
        # Some Logistics
        n=len(self.sys.x)
        m=len(self.sys.u)
        nC=len(self.sys.contact_points)
        # phi <= 0
        H1=np.hstack((phi_x,phi_u))
        h1=-phi+np.dot(phi_x,x)+np.dot(phi_u,u)
        # f_n >= 0 
        H2=np.zeros((1,m+n))
        h2=np.zeros((1,1))
        H2[:,n+m-2*nC+2*self.index]=-1
        # f_t= -mu*f_n
        H3=np.zeros((2,m+n))
        h3=np.zeros((2,1))
        H3[:,n+m-2*nC+2*self.index:n+m-2*nC+2*(self.index+1)]=np.array([[self.friction,1],[-self.friction,-1]])
        # v_psi <= -v_epsilon
        H4=np.hstack((v_psi_x,v_psi_u))
        h4=-v_epsilon-v_psi+np.dot(v_psi_x,x)+np.dot(v_psi_u,u)
        # Contact model
        if self.contact_model=="soft":
            # we have f_n = - K * phi - c * phi_dot --> K*phi+K*(phi_x,phi_u)(x,u)+c*phi_dot+f_n=0 
            H5=np.hstack((phi_x,phi_u))*self.K
            H5[n+m-2*nC+2*self.index]=1
            h5=-self.K*phi-self.damping*v_phi+self.K*(np.dot(phi_x,x)+np.dot(phi_u,u))
            H6=-H5
            h6=self.K*phi+self.damping*v_phi-self.K*(np.dot(phi_x,x)+np.dot(phi_u,u))
            H=np.vstack((H1,H2,H3,H4,H5,H6))
            h=np.vstack((h1,h2,h3,h4,h5,h6)).reshape(7,1)
            return (H,h)
        elif self.contact_model=="hard_implicit":
            # I only impose a maximum penetration constraint
            # phi>=-maximum_penetration: -(phi_x,phi_u)(x,u) + phi_x x+ phi_u u - phi <= maximum_penetration
            H5=-np.hstack((phi_x,phi_u))
            h5=phi-np.dot(phi_x,x)-np.dot(phi_u,u)+maximum_penetration
            H=np.vstack((H1,H2,H3,H4,H5))
            h=np.vstack((h1,h2,h3,h4,h5)).reshape(6,1)
            return (H,h)           
        else:
            raise ValueError("Unknown contact model")
        
        


def make_mode(cp,H,h,J_x,J_u,J,x0,u0):
    """
    We have the following:
        J*u+J_u u -J_u u0
    """
    assert 1==1 # It should be H*(x0,u0)<=h
    


    
    
    
    
    

































def extract_point(determiner,index):
    """
    If iis an array: take the value index
    
        
    If just a number 
    If inside list is a number, then take it
    If list:
        if first entry a number, then take the whole list! numpy it
    """
    N=len(determiner)
    d=[0]*N
    for i in range(N):
#        print i
#        print type(determiner[i])
        if type(determiner[i])==type(np.array(0)):    
            d[i]=determiner[i][index]
        elif type(determiner[i])==type(d):  # list
            if all([type(determiner[i][j]) in [type(0.0),type(0),type(np.float64(0))] for j in range(len(determiner[i]))]):
                d[i]=determiner[i]
            elif any([type(determiner[i][j])==type(np.array(0)) for j in range(len(determiner[i]))]):
                _d=[0]*len(determiner[i])
                for var in range(len(determiner[i])):
                    if type(determiner[i][var]) in [type(0.0),type(0),type(np.float64(0))]:
                        _d[var]=determiner[i][var]
                    else:
                        _d[var]=determiner[i][var][index]
                d[i]=_d
            else:
                raise NotImplementedError
        d[i]=np.array(d[i])
    return d