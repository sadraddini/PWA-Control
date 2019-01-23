#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 13:46:55 2019

@author: sadra
"""
from sympy import Symbol,pi,sin,cos,Function,diff,Matrix,symbols,lambdify
import numpy as np

from PWA_lib.Manipulation.contact_point import contact_point,extract_point
from PWA_lib.Manipulation.symbolic_system import symbolic_system

# test
mysystem=symbolic_system("my system")
t=mysystem.t
x=Function("x")(t)
y=Function("y")(t)
ux=Function("ux")(t)
uy=Function("uy")(t)
fn_1=Function("fn_1")(t)
ft_1=Function("ft_1")(t)

mysystem.x=[x,y,diff(x,t),diff(y,t)]
mysystem.u=[ux,uy,fn_1,ft_1]
phi=y-5
psi=x
J=0

# Introduce Contact Point
c1=contact_point(mysystem,0,phi,psi,J)
D=c1.get_determiners_symbolic()
D_lambda=mysystem.sys_lambdify(D)

x_sample=np.ones((7,4))
u_sample=np.ones((7,4))

D_n=mysystem.evaluate_handles(D_lambda,x_sample,u_sample)
z=extract_point(D_n,0)
(H,h)=c1.forces_sticking(z)
raise 1
(v_phi,v_psi,phi_x,phi_u,v_psi_x,v_psi_u)=c1.make_derivatives()


assert 1==0
c1.forces_no_contact(v_phi,v_psi,phi_x,phi_u,v_psi_x,v_psi_u,phi_0,psi_0,x0,u0)

def complete(D):
    pass
    