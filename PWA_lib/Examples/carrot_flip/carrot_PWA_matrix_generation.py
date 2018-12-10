# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 16:00:33 2018

@author: sadra
"""

from sympy import Symbol,pi,sin,cos,Function,diff,Matrix
import numpy as np

import time as time
from pypolycontain.lib.polytope import polytope

from carrot_symbolic_flip import carrot_symbolic
from carrot_PWA_mode_generation import GRADIENT,SUB

 

def force_ground(symbolic_carrot):
    J=symbolic_carrot.Jacobian_ground()
    phi,v_slide=symbolic_carrot.lambda_ground_curve()
    lambda_phi=-symbolic_carrot.K_ground*phi-symbolic_carrot.c_ground*diff(phi,symbolic_carrot.t) # Positive when penetrations is negative
    lambda_slide=-lambda_phi*symbolic_carrot.mu_ground
    # Three Forces
    f0=np.dot(J,Matrix([np.zeros((2,1))]).reshape(2,1))
    f1=np.dot(J,Matrix([[lambda_phi,lambda_slide]]).reshape(2,1)) 
    f2=np.dot(J,Matrix([[lambda_phi,-lambda_slide]]).reshape(2,1))
    A0=GRADIENT(f0,symbolic_carrot.X)
    B0=GRADIENT(f0,symbolic_carrot.U)
    A1=GRADIENT(f1,symbolic_carrot.X)
    B1=GRADIENT(f1,symbolic_carrot.U)
    A2=GRADIENT(f2,symbolic_carrot.X)
    B2=GRADIENT(f2,symbolic_carrot.U)
    return (A0,B0,A1,B1,A2,B2)

def force_left(symbolic_carrot):
    J=symbolic_carrot.Jacobian_left_finger()
    phi,v_slide=symbolic_carrot.lambda_left_finger_point()
    lambda_phi=-symbolic_carrot.K_finger*phi-symbolic_carrot.c_finger*diff(phi,symbolic_carrot.t) # Positive when penetrations is negative
    lambda_slide=-lambda_phi*symbolic_carrot.mu_finger
    # Three Forces
    f0=np.dot(J,Matrix([np.zeros((2,1))]).reshape(2,1))
    f1=np.dot(J,Matrix([[lambda_phi,lambda_slide]]).reshape(2,1)) 
    f2=np.dot(J,Matrix([[lambda_phi,-lambda_slide]]).reshape(2,1))
    A0=GRADIENT(f0,symbolic_carrot.X)
    B0=GRADIENT(f0,symbolic_carrot.U)
    A1=GRADIENT(f1,symbolic_carrot.X)
    B1=GRADIENT(f1,symbolic_carrot.U)
    A2=GRADIENT(f2,symbolic_carrot.X)
    B2=GRADIENT(f2,symbolic_carrot.U)
    return (A0,B0,A1,B1,A2,B2)

def force_right(symbolic_carrot):
    J=symbolic_carrot.Jacobian_right_finger()
    phi,v_slide=symbolic_carrot.lambda_right_finger()
    lambda_phi=-symbolic_carrot.K_finger*phi-symbolic_carrot.c_finger*diff(phi,symbolic_carrot.t) # Positive when penetrations is negative
    lambda_slide=-lambda_phi*symbolic_carrot.mu_finger
    # Three Forces
    f0=np.dot(J,Matrix([np.zeros((2,1))]).reshape(2,1))
    f1=np.dot(J,Matrix([[lambda_phi,lambda_slide]]).reshape(2,1)) 
    f2=np.dot(J,Matrix([[lambda_phi,-lambda_slide]]).reshape(2,1))
    A0=GRADIENT(f0,symbolic_carrot.X)
    B0=GRADIENT(f0,symbolic_carrot.U)
    A1=GRADIENT(f1,symbolic_carrot.X)
    B1=GRADIENT(f1,symbolic_carrot.U)
    A2=GRADIENT(f2,symbolic_carrot.X)
    B2=GRADIENT(f2,symbolic_carrot.U)
    return (A0,B0,A1,B1,A2,B2)   

O=carrot_symbolic("my carrot")
(A0_right,B0_right,A1_right,B1_right,A2_right,B2_right)=force_right(O)
(A0_left,B0_left,A1_left,B1_left,A2_left,B2_left)=force_left(O)
(A0_ground,B0_ground,A1_ground,B1_ground,A2_ground,B2_ground)=force_ground(O)