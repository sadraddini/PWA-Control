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
J=[2*x+2*y**2,x-y,0,3*y**3,-2,1,2*x**2+3*y,-y**2]

# Introduce Contact Point
c1=contact_point(mysystem,0,phi,psi,J,contact_model="hard_implicit")
D=c1.get_determiners_symbolic()
D_lambda=mysystem.sys_lambdify(D)
E=c1.get_determiners_symbolic_J()
E_lambda=mysystem.sys_lambdify(E)


x_sample=np.ones((7,4))*3
u_sample=np.ones((7,4))

L=c1.build_PWA_cells(x_sample,u_sample,epsilon_confidence=np.ones((8,1))*0.1)
assert 1==0

#D_n=mysystem.evaluate_handles(D_lambda,x_sample,u_sample)
#E_n=mysystem.evaluate_handles(E_lambda,x_sample,u_sample)
#
#z=extract_point(D_n,0)
#q=extract_point(E_n,1)
#
#J_n,J_t,J_n_x,J_t_x=mysystem._extract(q)
#
#(H1,h1)=c1.forces_no_contact(z)
#(H2,h2)=c1.forces_sticking(z)
#(H3,h3)=c1.force_slide_positive(z)
#(H4,h4)=c1.force_slide_negative(z)