#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 15:50:54 2019

@author: sadra
"""
from sympy import Symbol,pi,sin,cos,Function,diff,Matrix,symbols,lambdify
import numpy as np

from PWA_lib.trajectory.system import linear_cell
from pypolycontain.lib.polytope import polytope
from pypolycontain.utils.redundancy_reduction import canonical_polytope


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
    def __init__(self,sys,index,phi,psi,J,friction=0.4,K=100,damping=1,name="Contact Point ",contact_model="soft"):
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
        Output: 
            list of J values and J_n_x and then J_t_x
        """
        n=len(self.sys.x)
        print len(self.J)
        assert len(self.J)==2*n # This is important!
        J_x=[ [diff(J_f,var) for var in self.sys.x] for J_f in self.J]
        return self.J+J_x[0:n]+J_x[n:]

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
#        print H2.shape,H1.shape
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
        


    """
    These are the functions to construct stuff :)  
    Functions with preceding "_" indicate for inner use
    """    
    
    def build_PWA_cells(self,x_sample,u_sample,epsilon_confidence,big_eps=100):
        n,m,nC=len(self.sys.x),len(self.sys.u),len(self.sys.contact_points)
        N=x_sample.shape[0]
        assert N==u_sample.shape[0]
        D=self.get_determiners_symbolic()
        D_lambda=self.sys.sys_lambdify(D)
        E=self.get_determiners_symbolic_J()
        E_lambda=self.sys.sys_lambdify(E)
        D_n=self.sys.evaluate_handles(D_lambda,x_sample,u_sample)
        E_n=self.sys.evaluate_handles(E_lambda,x_sample,u_sample)
        print "************ u_sample is",u_sample
        # Construct Epsilon Confidence
        # Only for non affecting controls
        big_epsilon=np.array([y not in [m+n-2*nC+2*self.index,m+n-2*nC+2*self.index+1] and y>n  for y in range(n+m)]).reshape((n+m,1))*10
        epsilon_confidence=epsilon_confidence+big_epsilon
        print "************ u_sample is",u_sample        
        list_of_outputs=[None]*N
        for i in range(N):
            p={}
            z=extract_point(D_n,i)
            q=extract_point(E_n,i)
            J_n,J_t,J_n_x,J_t_x=self.sys._extract(q)
            (H1,h1)=self.forces_no_contact(z)
            (H2,h2)=self.forces_sticking(z)
            (H3,h3)=self.force_slide_positive(z)
            (H4,h4)=self.force_slide_negative(z)
            p["NC"]=polytope(H1,h1)
            p["ST"]=polytope(H2,h2)
            p["SP"]=polytope(H3,h3)
            p["SN"]=polytope(H4,h4)
            C=self.construct_PWA_cell(p,[J_n,J_t,J_n_x,J_t_x],x_sample[i,:].reshape(n,1),u_sample[i,:].reshape(m,1),epsilon_confidence)
            print "************ u_sample is",u_sample
#            list_of_outputs[n]=(p,[J_n,J_t,J_n_x,J_t_x])
            list_of_outputs[i]=C
            
        return list_of_outputs
            
            
            
            
            
    

    def construct_linear_cell(self,H,h,Jacobian_numbers,x0,u0,epsilon_confidence):
        """
        We have the following:
            A= J_t_x \lambda_n + J_n_x \lambda_t
        """
        n,m,nC=len(self.sys.x),len(self.sys.u),len(self.sys.contact_points)
        fn0,ft0=u0[m-2*nC+2*self.index,0],u0[m-2*nC+2*self.index+1,0]
        J_n,J_t,J_n_x,J_t_x=Jacobian_numbers[0:4]
#        print J_n
#        print J_t
#        print J_n_x
#        print J_t_x
        A= J_n_x*fn0 + J_t_x*ft0
        B=np.zeros((n,m))
        B[:,m-2*nC+2*self.index:m-2*nC+2*(self.index+1)]=np.hstack((J_n.reshape(n,1),J_t.reshape(n,1)))
#        print x0.shape,A.shape
        c=-np.dot(A,x0)
        # Bring epsilon in!
        H_eps=np.vstack((np.eye(n+m),-np.eye(n+m)))
        xu0=np.vstack((x0,u0))
        h_eps=np.vstack((xu0+epsilon_confidence,-xu0+epsilon_confidence))
        print "h_eps.shape=",xu0.shape,epsilon_confidence.shape,h_eps.shape
        H_new=np.vstack((H,H_eps))
        h_new=np.vstack((h,h_eps))
        (H_can,h_can)=canonical_polytope(H_new,h_new)
        assert 1==1
#        print H,h,h-np.dot(H,xu0)
#        assert all(np.dot(H,xu0)<=h)==True
#        raise NotImplementedError
        return linear_cell(A,B,c,polytope(H_can,h_can))
    
    def construct_PWA_cell(self,polytopes_dict,Jacobian_numbers,x0,u0,epsilon_confidence):
        n,m,nC=len(self.sys.x),len(self.sys.u),len(self.sys.contact_points)
        C={} # These are the cells
        fn,ft=u0[m-2*nC+2*self.index,0],u0[m-2*nC+2*self.index+1,0]
        F=self.generate_contact_forces(fn,ft)
        for mode in ["NC","ST","SP","SN"]:
            p=polytopes_dict[mode]
            H,h=p.H,p.h
            fn_new,ft_new=F[mode]
            _u=np.array(u0)
            _u[m-2*nC+2*self.index,0],_u[m-2*nC+2*self.index+1,0]=fn_new,ft_new
#            raise NotImplementedError # I have to figure out epsilon_confidence here!
            C[mode]=self.construct_linear_cell(H,h,Jacobian_numbers,x0,_u,epsilon_confidence)
#        raise NotImplementedError
        return C
    
    def generate_contact_forces(self,fn,ft,fn_dummy=1):
        """
        Given a contact force, we generate closest contact forces in all possible contact modes.
        This contact forces are later used as the nominal forces for linearization.
        Arguments:
            fn= nominal force
            ft= tangential force
            fn_dummy= if nominal force is zero (no contact), what force should be considered if contact unexpectedely occurs?
        """
        F={}
        if self.friction*fn>abs(ft): # Sticking
            F["NC"]=(0,0)
            F["ST"]=(fn,ft)
            F["SP"]=(fn,self.friction*fn)
            F["SN"]=(fn,-self.friction*fn)
        elif self.friction*fn<ft: # Sldiing highly positive!
            F["NC"]=(0,0)
            F["ST"]=(fn,self.friction*fn)
            F["SP"]=(fn,ft)
            F["SN"]=(fn,-self.friction*fn)        
        elif -self.friction*fn>-ft: # Sldiing highly negative!
            F["NC"]=(0,0)
            F["ST"]=(fn,-self.friction*fn)
            F["SP"]=(fn,self.friction*fn)
            F["SN"]=(fn,ft)    
        elif fn==0: # No Contact
            assert ft==0
            F["NC"]=(fn,ft)
            F["ST"]=(fn_dummy,0)
            F["SP"]=(fn_dummy,self.friction*fn_dummy)
            F["SN"]=(fn_dummy,-self.friction*fn_dummy) 
        return F
    
          
    


    
    
    
    
    

































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
        elif type(determiner[i]) in [type(0.0),type(0),type(np.float64(0))]:
            d[i]=determiner[i]
        else:
            raise NotImplementedError
        d[i]=np.array(d[i])
    return d