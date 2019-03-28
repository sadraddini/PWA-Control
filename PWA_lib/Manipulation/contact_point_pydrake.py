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
    def __init__(self,sys,phi,psi,J,friction=0.4,K=100,damping=1,name="Contact Point",contact_model="soft"):
        self.sys=sys # What system it belongs to
        self.phi=phi # Penetration position, symbolic expression
        self.psi=psi # Sliding position, symbolic expression
        self.J=J # Jacobian for the forces
        self.index=len(self.sys.contact_points) # Must be an integer between [0,N), N the number of contact points
        self.u_lambda=np.array([sym.Variable("lambda%d %s"%(d,name)) for d in self.sys.d])
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
        
    def _contact_geometry_no_contact(self):
        """
        phi > 0, psi does not matter 
        """
        zeta=phi
        