#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 19:47:22 2019

@author: sadra
"""
import numpy as np
import pydrake.symbolic as sym

from PWA_lib.trajectory.system import linear_cell
from pypolycontain.lib.polytope import polytope

"""
What a contact point is
"""

class contact_point_symbolic_2D:
    """
        Each contact point:
            * Forces as control inputs which should satisfy Friction and Complementarity Constraints
            * 2D so far: 4 Possibilites: No contact, slide left, slide right, sticking
    """
    
    """ SYMBOLIC:"""
    def __init__(self,sys,phi,psi,J,friction=0.4,K=100,damping=1,name="Contact Point ",contact_model="hard"):
        self.sys=sys # What system it belongs to
        self.phi=phi # Penetration position, symbolic expression
        self.psi=psi # Sliding position, symbolic expression
        self.J=J # Jacobian for the forces
        self.index=len(self.sys.contact_points) # Must be an integer between [0,N), N the number of contact points
#        self.u_lambda=np.array([sym.Variable("lambda%d %s"%(d,name)) for d in self.sys.d])
        # The bounds
        self.psi_max=None
        self.psi_min=None
        # Dynamical properties
        self.polytope={}
        self.friction=friction # Friction
        self.contact_model=contact_model 
        # Only used for soft contact
        self.damping=damping # Damper coefficient, used for soft contact
        self.K=K # Spring coefficient, used for soft contact
        self.sys.contact_points.append(self)
        self.name=name+str(len(self.sys.contact_points))+" contact model: "+self.contact_model
        
    def __repr__(self):
        return self.name
        
    def _contact_geometry_no_contact(self,tolerance=10**-6):
        """
        phi > 0, psi does not matter
        """
        zeta=np.array([-self.phi])
        return self.sys._polyhedrize_constraints_symbolic(zeta)
    
    def _contact_geometry_sticking(self,tolerance_phi=10**-6,tolerance_psi=10**-6):
        """
        phi < -tolerance_phi, -tolerance_psi<psi_dot<tolerance_psi
        """
        psi_dot=self.sys._time_derivative(np.array([self.psi]))[0]
        zeta=np.array([self.phi+tolerance_phi,-psi_dot-tolerance_psi,psi_dot-tolerance_psi])
        return self.sys._polyhedrize_constraints_symbolic(zeta)
    
    def _contact_geometry_sliding_positive(self,tolerance_phi=10**-6,tolerance_psi=10**-6):
        """
        phi <=-tolerance_phi psi_dot>tolerance_psi
        """
        psi_dot=self.sys._time_derivative(np.array([self.psi]))[0]
        zeta=np.array([self.phi+tolerance_phi,-psi_dot+tolerance_psi])
        return self.sys._polyhedrize_constraints_symbolic(zeta)
    
    def _contact_geometry_sliding_negative(self,tolerance_phi=10**-6,tolerance_psi=10**-6):
        """
        phi <=-tolerance_phi psi_dot<-tolerance_psi
        """
        psi_dot=self.sys._time_derivative(np.array([self.psi]))[0]
        zeta=np.array([self.phi+tolerance_phi,psi_dot+tolerance_psi])
        return self.sys._polyhedrize_constraints_symbolic(zeta)
    
    def _contact_geometery_all(self):
        H,h={},{}
        H["no_contact"],h["no_contact"]=self._contact_geometry_no_contact()
        H["sticking"],h["sticking"]=self._contact_geometry_sticking()
        H["sliding_positive"],h["sliding_positive"]=self._contact_geometry_sliding_positive()
        H["sliding_negative"],h["sliding_negative"]=self._contact_geometry_sliding_negative()
        self.H=H
        self.h=h
    
    """
    Numerical Functions
    """
    
    def time_stepping_geometrical_constraint(self,H,h,dynamical_matrices):
        """
        Arguments:
            * H, h: floats
            * dynamical_matrices: floats
        """
        A,B_u,B_lambda,c=dynamical_matrices.A,dynamical_matrices.B_u,dynamical_matrices.B_lambda,dynamical_matrices.c
        H_x,H_u,H_lambda=H[:,range(self.sys.x.shape[0])],\
                        H[:,range(self.sys.x.shape[0],self.sys.x.shape[0]+self.sys.u.shape[0])],\
                        H[:,range(self.sys.x.shape[0]+self.sys.u.shape[0],H.shape[1])]
        E_x,E_u,E_lambda,e=np.dot(H_x,A),np.dot(H_x,B_u)+H_u,np.dot(H_x,B_lambda)+H_lambda,h-np.dot(H_x,c)
        E=np.hstack((E_x,E_u,E_lambda))
        return (E,e)
    
    def Evaluate_polyhedral_matrices(self,dynamical_matrices,Eta):
        """
        H,h are dictionary of different volumes
        """
        H_n,h_n={},{}
        for contact_mode in ["no_contact","sticking","sliding_positive","sliding_negative"]:
            print contact_mode
            print self.H[contact_mode].shape,type(self.H[contact_mode])
            H=sym.Evaluate(self.H[contact_mode],Eta)
            h=sym.Evaluate(self.h[contact_mode],Eta)
            H_n[contact_mode],h_n[contact_mode]=self.time_stepping_geometrical_constraint(H,h,dynamical_matrices)
        return H_n,h_n
