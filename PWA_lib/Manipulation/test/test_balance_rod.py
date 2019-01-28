#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 14:01:34 2019

@author: sadra
"""

from sympy import Symbol,pi,sin,cos,sqrt,Function,diff,Matrix,symbols,lambdify
import numpy as np

from PWA_lib.Manipulation.contact_point import contact_point,extract_point
from PWA_lib.Manipulation.symbolic_system import symbolic_system

# test
mysystem=symbolic_system("my system")
t=mysystem.t
x=Function("x")(t)
y=Function("y")(t)
theta=Function("theta")(t)

x_1=Function("x_1")(t)
y_1=Function("y_1")(t)

x_2=Function("x_2")(t)
y_2=Function("y_2")(t)


fn_1=Function("fn_1")(t)
ft_1=Function("ft_1")(t)

fn_2=Function("fn_2")(t)
ft_2=Function("ft_2")(t)

mysystem.x=[x,y,theta,x_1,y_1,x_2,y_2,diff(x,t),diff(y,t),diff(theta,t)]
mysystem.u=[diff(x_1,t),diff(y_1,t),diff(x_2,t),diff(y_2,t),fn_1,ft_1,fn_2,ft_2]

J1=[0,0,0,0,0,0,0,-sin(theta),cos(theta),x_1]+[0,0,0,0,0,0,0,cos(theta),sin(theta),0]
phi=y_1
psi=x_1
c1=contact_point(mysystem,0,phi,psi,J1,contact_model="hard_implicit")



J2=[0,0,0,0,0,0,0,-sin(theta),cos(theta),-x_2]+[0,0,0,0,0,0,0,cos(theta),sin(theta),0]
phi=y_2
psi=x_2
c2=contact_point(mysystem,1,phi,psi,J2,contact_model="hard_implicit")

#D=c1.get_determiners_symbolic()
#D_lambda=mysystem.sys_lambdify(D)
#E=c1.get_determiners_symbolic_J()
#E_lambda=mysystem.sys_lambdify(E)


x_sample=np.array(([0,0,0,1,0,-1,0,0,0,0])).reshape(1,10)
u_sample=np.array(([0,0,0,0,1,0,0,0])).reshape(1,8)

S1=c1.build_PWA_cells(x_sample,u_sample,epsilon_confidence=np.ones((18,1))*0.1)
S2=c2.build_PWA_cells(x_sample,u_sample,epsilon_confidence=np.ones((18,1))*0.1)


assert 1==0
h=0.01


# Introduce Contact Point 1
phi=y-3*sin(theta)+5
psi=x+2*cos(theta)+x_1
J=2*x*sin(theta)+y
c1=contact_point(mysystem,0,phi,psi,J,contact_model="hard_implicit")
# Introduce Contact Point 1
phi=y-3*sin(theta)+5
psi=x+2*cos(theta)+x_1
J=2*x*sin(theta)+y
c2=contact_point(mysystem,1,phi,psi,J,contact_model="hard_implicit")




D=c1.get_determiners_symbolic()
D_lambda=mysystem.sys_lambdify(D)

x_sample=np.ones((7,len(mysystem.x)))*0
u_sample=np.ones((7,len(mysystem.u)))

D_n=mysystem.evaluate_handles(D_lambda,x_sample,u_sample)
z=extract_point(D_n,0)
(H1,h1)=c1.forces_no_contact(z)
(H2,h2)=c1.forces_sticking(z)
(H3,h3)=c1.force_slide_positive(z)
(H4,h4)=c1.force_slide_negative(z)