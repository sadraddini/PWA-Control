#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 19:08:09 2019

@author: sadra
"""
import numpy as np
from scipy.linalg import block_diag as blk

import pydrake.symbolic as sym

from PWA_lib.Manipulation.contact_point_pydrake import contact_point_symbolic_2D
from PWA_lib.Manipulation.system_symbolic_pydrake import system_symbolic,_evaluate_polyhedral_matrices


mysystem=system_symbolic("Balacning a bar", 2)
x,y,theta=sym.Variable("x"),sym.Variable("y"),sym.Variable("theta")
x_1,y_1,x_2,y_2=sym.Variable("x_1"),sym.Variable("y_1"),sym.Variable("x_2"),sym.Variable("y_2")
mysystem.q_o=np.array([x,y,theta])
mysystem.q_m=np.array([x_1,y_1,x_2,y_2])

# Dynamics:
mysystem.C=np.zeros((3,3))
g=9.8
mysystem.tau_g=np.array([0,-g,0])
mysystem.B=np.zeros((3,0))


I=1/3.0
mysystem.M=blk(*[1,1,I])
mysystem.M_inv=np.linalg.inv(mysystem.M)


psi_1= (x_1-x)*sym.cos(theta) + (y_1-y)*sym.sin(theta)
phi_1=-( (x_1-x)*sym.sin(-theta) + (y_1-y)*sym.cos(theta)  )
J_1n=np.array([-sym.sin(theta),sym.cos(theta),psi_1]).reshape(3,1)
J_1t=np.array([sym.cos(theta),sym.sin(theta),0]).reshape(3,1)
J_1=np.hstack((J_1n,J_1t))

psi_2= (x_2-x)*sym.cos(theta) + (y_2-y)*sym.sin(theta)
phi_2=-( (x_2-x)*sym.sin(-theta) + (y_2-y)*sym.cos(theta)  )
J_2n=np.array([-sym.sin(theta),sym.cos(theta),psi_2]).reshape(3,1)
J_2t=np.array([sym.cos(theta),sym.sin(theta),0]).reshape(3,1)
J_2=np.hstack((J_2n,J_2t))

mysystem.J=np.hstack((J_1,J_2))
mysystem.u_lambda=np.array([sym.Variable("lambda_1n"),sym.Variable("lambda_1t"),sym.Variable("lambda_2n"),sym.Variable("lambda_2t")])

mysystem._build_state_space()
mysystem._linearize_dynamics_symbolic()

Eta={x:0,y:0,theta:0,x_1:1,x_2:-1,y_1:0,y_2:0,
     mysystem.v_o[0]:0,mysystem.v_o[1]:0,mysystem.v_o[2]:0,
     mysystem.u_lambda[0]:1,mysystem.u_lambda[1]:0,mysystem.u_lambda[2]:1,mysystem.u_lambda[3]:0,
     mysystem.u_m[0]:1,mysystem.u_m[1]:0,mysystem.u_m[2]:1,mysystem.u_m[3]:0,
     mysystem.h:0.01}


dynamical_matrices=mysystem._evaluate_dynamical_matrices(Eta)


C1=contact_point_symbolic_2D(mysystem,phi=phi_1,psi=psi_1,J=J_1,name="contact point")
C1._contact_geometery_all()
H,h=C1.Evaluate_polyhedral_matrices(dynamical_matrices,Eta)