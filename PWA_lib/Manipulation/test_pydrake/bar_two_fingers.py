#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 19:08:09 2019

@author: sadra
"""
import numpy as np
from scipy.linalg import block_diag as blk

import pydrake.symbolic as sym

from pypolycontain.lib.zonotope import zonotope

from PWA_lib.Manipulation.contact_point_pydrake import contact_point_symbolic_2D
from PWA_lib.Manipulation.system_symbolic_pydrake import system_symbolic
from PWA_lib.Manipulation.system_numeric import system_hard_contact_PWA_numeric as system_numeric
from PWA_lib.Manipulation.system_numeric import environment
from PWA_lib.trajectory.poly_trajectory import point_trajectory_tishcom


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


M=1
I=M/3.0
mysystem.M=blk(*[M,M,I])
mysystem.M_inv=np.linalg.inv(mysystem.M)


psi_1= (x_1-x)*sym.cos(theta) + (y_1-y)*sym.sin(theta)
phi_1=-( (x_1-x)*sym.sin(-theta) + (y_1-y)*sym.cos(theta)  )
J_1n=np.array([-sym.sin(theta),sym.cos(theta),psi_1]).reshape(3,1)
J_1t=np.array([sym.cos(theta),sym.sin(theta),0]).reshape(3,1)
J_1=np.hstack((J_1n,J_1t))
C1=contact_point_symbolic_2D(mysystem,phi=phi_1,psi=psi_1,J=J_1,friction=0.3,name="contact point")


psi_2= (x_2-x)*sym.cos(theta) + (y_2-y)*sym.sin(theta)
phi_2=-( (x_2-x)*sym.sin(-theta) + (y_2-y)*sym.cos(theta)  )
J_2n=np.array([-sym.sin(theta),sym.cos(theta),psi_2]).reshape(3,1)
J_2t=np.array([sym.cos(theta),sym.sin(theta),0]).reshape(3,1)
J_2=np.hstack((J_2n,J_2t))
C2=contact_point_symbolic_2D(mysystem,phi=phi_2,psi=psi_2,J=J_2,friction=0.3,name="contact point")


mysystem.build_and_linearize()


Eta=environment(0)
Eta.dict_of_values={x:0,y:0,theta:0,x_1:1,x_2:-1,y_1:0,y_2:0,
     mysystem.v_o[0]:0,mysystem.v_o[1]:0,mysystem.v_o[2]:0,
     mysystem.u_lambda[0]:M*g/2,mysystem.u_lambda[1]:0,mysystem.u_lambda[2]:M*g/2,mysystem.u_lambda[3]:0,
     mysystem.u_m[0]:1,mysystem.u_m[1]:0,mysystem.u_m[2]:1,mysystem.u_m[3]:0,
     mysystem.h:0.04}

epsilon=np.array([2,2,2,1,1,1,1,5,5,5,1,1,1,1,50,1000,50,1000]).reshape(18,1)
sys=system_numeric(mysystem)
sys.add_environment(Eta,epsilon)

y_goal=0.5
#x_goal=np.array([0,y_goal,0.0,1,y_goal,-1,y_goal,0,0,0]).reshape(10,1)
x_goal=np.array([0,0,0.0,1,0,-1,0,0,0,0]).reshape(10,1)
sys.goal=zonotope(x_goal.reshape(10,1),np.diag([0,0,0,0,0,0,0,0,0,0]))
x0=np.array([0.0,0,0.0,1,0.00,-1,-0.0,1.2,0,0]).reshape(10,1)
T=20
x,u,u_lambda,x_tishcom,x_time=point_trajectory_tishcom(sys,x0,[sys.goal],T,optimize_controls_indices=[0,1,2,3])

def env_q(x):
    return {mysystem.q[i]:x[i] for i in range(len(mysystem.q))}

phi_1_n={t:phi_1.Evaluate(env_q(x[t])) for t in range(T+1)}
phi_2_n={t:phi_2.Evaluate(env_q(x[t])) for t in range(T+1)}
psi_1_n={t:psi_1.Evaluate(env_q(x[t])) for t in range(T+1)}
psi_2_n={t:psi_2.Evaluate(env_q(x[t])) for t in range(T+1)}

"""
Visualization
"""
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

def animate(X,ax1):
    x,y,theta,x_1,y_1,x_2,y_2=X[0:7]
    L=1.5
    b=0.3
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
    ax1.set_ylim([-0.7,1.8])
    fig.gca().set_aspect('equal')
    bar=[patches.Polygon(np.array([[x_left,x_left_bar,x_right_bar,x_right],[y_left,y_left_bar,y_right_bar,y_right]]).reshape(2,4).T, True)]
    ax1.add_collection(PatchCollection(bar,color=(0.8,0.3,0.4),alpha=0.8,edgecolor=(0,0,0)))
    ax1.plot([x_1,x_2],[y_1-0.07,y_2-0.07],'^',linewidth=3,markersize=10)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.5)
    
for t in range(T+1):
    fig,ax = plt.subplots()
#    ax.set_title("Balancing t=%d \n Right finger: %s \n Left Finger: %s"%(t,mode_name[list_of_modes[t][1]],mode_name[list_of_modes[t][2]]))
    animate(x[t],ax)
    if t==-1:
        ax.arrow(-1.8, 0.1, 0.3, 0.0, head_width=0.3, head_length=0.3, linewidth=10, fc='k', ec='k')
    fig.savefig('figures/bar_%d.png'%t, dpi=100)