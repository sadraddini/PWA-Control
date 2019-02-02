#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 15:36:51 2019

@author: sadra
"""


from sympy import Symbol,pi,sin,cos,sqrt,Function,diff,Matrix,symbols,lambdify
import numpy as np

from PWA_lib.Manipulation.contact_point import contact_point,extract_point
from PWA_lib.Manipulation.symbolic_system import symbolic_system

# test
mysystem=symbolic_system("my carrot")
t=mysystem.t
x=Function("x")(t)
y=Function("y")(t)
theta=Function("theta")(t)

D=Function("D")(t)
zeta=Function("zeta")(t)
sep=Function("sep")(t)
pen=Function("pen")(t)


fn_left=Function("fn_left")(t)
ft_left=Function("ft_left")(t)

fn_right=Function("fn_right")(t)
ft_right=Function("ft_right")(t)

fn_ground=Function("fn_ground")(t)
ft_ground=Function("ft_ground")(t)

mysystem.x=[x,y,theta,diff(x,t),diff(y,t),diff(theta,t),D,zeta,sep,pen]
mysystem.u=[diff(D,t),diff(zeta,t),diff(sep,t),diff(pen,t),\
            fn_left,ft_left,fn_right,ft_right,fn_ground,ft_ground]


R=1
p=4*R/(3*np.pi)

"""LEFT FINGER """
J1=[0,0,0,-cos(theta),sin(theta),D,0,0,0,0]+[0,0,0,-sin(theta),-cos(theta),p,0,0,0,0]
phi=pen
psi=D
c1=contact_point(mysystem,0,phi,psi,J1,contact_model="hard_implicit")

"""RIGHT FINGER """
J2=[0,0,0,-sin(zeta),cos(zeta),-p*D*sin(zeta-theta),0,0,0,0]+[0,0,0,cos(zeta),sin(zeta),R-p*D*cos(zeta-theta),0,0,0,0]
phi=sep-D*sin(zeta-theta)-R-pen*sin(zeta-theta)
psi=R*(zeta-theta)
c2=contact_point(mysystem,1,phi,psi,J2,contact_model="hard_implicit")

"""GROUND """
#D=c1.get_determiners_symbolic()
#D_lambda=mysystem.sys_lambdify(D)
#E=c1.get_determiners_symbolic_J()
#E_lambda=mysystem.sys_lambdify(E)
J3=[0,0,0,-sin(zeta),cos(zeta),-p*D*sin(zeta-theta),0,0,0,0]+[0,0,0,cos(zeta),sin(zeta),R-p*D*cos(zeta-theta),0,0,0,0]
phi=sep-D*sin(zeta-theta)-R-pen*sin(zeta-theta)
psi=R*(zeta-theta)
c3=contact_point(mysystem,2,phi,psi,J3,contact_model="hard_implicit")





gravity=9.8
h=0.01
f=[diff(x,t),diff(y,t),diff(theta,t),0,0,0,0,0,-gravity,0]
g=np.zeros((10,8))
g[3,0]=1
g[4,1]=1
g[5,2]=1
g[6,3]=1

x_sample=np.array(([0,0,0,1,0,-1,0,0,0,0])).reshape(1,10)
u_sample=np.array(([0,0,0,0,gravity/2.0,0,gravity/2.0,0])).reshape(1,8)

#x_sample=np.random.random((5,10))
#u_sample=np.random.random((5,8))

mysystem.f=f
mysystem.get_contact_free_dynamics(x_sample,u_sample)

A_dict,c_dict=mysystem.get_contact_free_dynamics(x_sample,u_sample)


h=0.01

A0=np.eye(10)+h*A_dict[0]
B0=h*g
c0=h*np.array(extract_point(c_dict,0)).reshape(10,1)
S={}
S[1]=c1.build_PWA_cells(x_sample,u_sample,epsilon_confidence=np.ones((18,1))*0.8)
#assert False
S[2]=c2.build_PWA_cells(x_sample,u_sample,epsilon_confidence=np.ones((18,1))*0.8)

import numpy as np

# My modules
from pypolycontain.lib.zonotope import zonotope
from pypolycontain.lib.polytope import polytope,Box
from pypolycontain.utils.redundancy_reduction import check_empty_polytope


# Internal imports
from PWA_lib.trajectory.system import system,linear_cell

sys=system()
sys.name="balancing rod"


sys.A[0,0]=A0
sys.B[0,0]=B0
sys.c[0,0]=c0
sys.C[0,0]=Box(18,10)

for i in [10,11,12,13]+[28,29,30,31]:
    sys.C[0,0].h[i,0]=1

for c_point in [1,2]:
    for mode in ["NC","ST","SP","SN"]:
        sys.A[c_point,mode]=h*S[c_point][0][mode].A
        sys.B[c_point,mode]=h*S[c_point][0][mode].B
        sys.c[c_point,mode]=h*S[c_point][0][mode].c
        sys.C[c_point,mode]=S[c_point][0][mode].p


sys.goal=zonotope(x_sample.reshape(10,1),0.00*np.eye(10))

sys.n=10
sys.m=8
sys.list_of_sum_indices=[0,1,2]
sys.list_of_modes={}
sys.list_of_modes[0]=[0]
sys.list_of_modes[1]=["NC","ST","SP","SN"]
sys.list_of_modes[2]=["NC","ST","SP","SN"]

sys.build()

sys.build_cells()

import pickle
pickle.dump(sys,open("sadra_bar.pkl","w"))