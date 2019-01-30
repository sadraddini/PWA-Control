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



J2=[0,0,0,0,0,0,0,-sin(theta),cos(theta),x_2]+[0,0,0,0,0,0,0,cos(theta),sin(theta),0]
phi=y_2
psi=x_2
c2=contact_point(mysystem,1,phi,psi,J2,contact_model="hard_implicit")

#D=c1.get_determiners_symbolic()
#D_lambda=mysystem.sys_lambdify(D)
#E=c1.get_determiners_symbolic_J()
#E_lambda=mysystem.sys_lambdify(E)

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
from PWA_lib.trajectory.poly_trajectory import point_trajectory,polytopic_trajectory_given_modes

sys=system()
sys.name="balancing rod"


sys.A[0,0]=A0
sys.B[0,0]=B0
sys.c[0,0]=c0
sys.C[0,0]=Box(18,10)

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

#sys.build_cells()

x0=np.array([0,0,0,1,0.1,-1,0.00,0,0,0]).reshape(10,1)
#u0=np.array(([0,0,0,0,gravity/2,0,gravity/2,0])).reshape(8,1)
u0=np.array(([0,0,0,0,0,0,gravity/2,0])).reshape(8,1)

T=40
(x_n,u,delta_PWA,mu,flag)=point_trajectory(sys,x0,[sys.goal],T,eps=0)

list_of_cells=[]
for t in range(T):
    mode=tuple([i for n in sys.list_of_sum_indices for i in sys.list_of_modes[n] if delta_PWA[t,n,i]==1]) 
    print mode
    list_of_cells.append(sys.cell[mode])



#A=A0+sys.A[1,"ST"]+sys.A[2,"ST"]
#B=B0+sys.B[1,"ST"]+sys.B[2,"ST"]
#c=c0+sys.c[1,"ST"]+sys.c[2,"ST"]
#np.dot(A,x0)
#x_future=np.dot(A,x0)+np.dot(B,u0)+c
#
#xu0=np.vstack((x0,u0))
#sys.C[1,"ST"].if_inside(xu0)


import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

plt.plot([x_n[t][2] for t in range(T)])

def animate(X,ax1):
    x,y,theta,x_1,y_1,x_2,y_2=X[0:7]
    L=1.5
    b=0.2
    # x
    x_left=x-L*np.cos(theta)
    x_left_bar=x-L*np.cos(theta)-b*np.sin(theta)
    x_right=x+L*np.cos(theta)
    x_right_bar=x+L*np.cos(theta)-b*np.sin(theta)
    # y
    y_left=y-L*np.sin(theta)
    y_left_bar=y_left+b*np.cos(theta)
    y_right=y+L*np.sin(theta)
    y_right_bar=y_right+b*np.cos(theta)
    # Good now plot
    ax1.set_xlabel("x",fontsize=20)
    ax1.set_ylabel("y",fontsize=20)
    ax1.set_xlim([-2*L,2*L])
    ax1.set_ylim([-0.2,0.5])
#    fig.gca().set_aspect('equal')
    bar=[patches.Polygon(np.array([[x_left,x_left_bar,x_right_bar,x_right],[y_left,y_left_bar,y_right_bar,y_right]]).reshape(2,4).T, True)]
    ax1.add_collection(PatchCollection(bar,color=(0.5,0.3,0),alpha=0.8,edgecolor=(0,0,0)))
    xf_1=x+x_1*np.cos(theta)+y_1*np.sin(theta)
    xf_2=x+x_2*np.cos(theta)+y_2*np.sin(theta)
    yf_1=x+x_1*np.sin(theta)-y_1*np.cos(theta)
    yf_2=x+x_2*np.sin(theta)-y_2*np.cos(theta)
    ax1.plot([xf_1,xf_2],[yf_1,yf_2],'o')
    
for t in range(T):
    fig,ax = plt.subplots()
    animate(x_n[t],ax)
    fig.savefig('figures/bar_%d.png'%t, dpi=100)