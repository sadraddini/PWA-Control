#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 13:56:52 2019

@author: sadra
"""
# Numpy and Scipy
import numpy as np
from scipy.linalg import block_diag as blk

import pydrake.symbolic as sym

from pypolycontain.lib.zonotope import zonotope

from PWA_lib.Manipulation.contact_point_pydrake import contact_point_symbolic_2D
from PWA_lib.Manipulation.system_symbolic_pydrake import system_symbolic
from PWA_lib.Manipulation.system_numeric import system_hard_contact_PWA_numeric as system_numeric
from PWA_lib.Manipulation.system_numeric import environment
from PWA_lib.trajectory.poly_trajectory import point_trajectory_tishcom


mysystem=system_symbolic("Picking Up a Box", 2)
x,y,theta=sym.Variable("x"),sym.Variable("y"),sym.Variable("theta")
x_1,y_1,x_2,y_2=sym.Variable("x_1"),sym.Variable("y_1"),sym.Variable("x_2"),sym.Variable("y_2")
mysystem.q_o=np.array([x,y,theta])
mysystem.q_m=np.array([x_1,y_1,x_2,y_2])

a=1 # The half of the length of the vertical side
b=1 # The half the length of the horizonatl side
# Dynamics:
mysystem.C=np.zeros((3,3))
g=9.8
mysystem.tau_g=np.array([0,-g,0])
mysystem.B=np.zeros((3,0))

# Mass Matrix
M=1
I=M/3.0*a**2
mysystem.M=blk(*[M,M,I])
mysystem.M_inv=np.linalg.inv(mysystem.M)

# Contact Point 1: The ground with left corner
phi_1= y - a * sym.cos(theta) - b * sym.sin(theta)
psi_1= x + a * sym.sin(theta) - b * sym.cos(theta)
J_1n=np.array([0,1,-b*sym.cos(theta)-a*sym.sin(theta)]).reshape(3,1)
J_1t=np.array([1,0,a*sym.cos(theta)+b*sym.sin(theta)]).reshape(3,1)
J_1=np.hstack((J_1n,J_1t))
C1=contact_point_symbolic_2D(mysystem,phi=phi_1,psi=psi_1,J=J_1,friction=0.5,name="contact point")


# Contact Point 2: The ground with right corner
phi_2= y - a * sym.cos(theta) + b * sym.sin(theta)
psi_2= x + a * sym.sin(theta) + b * sym.cos(theta)
J_2n=np.array([0,1,b*sym.cos(theta)-a*sym.sin(theta)]).reshape(3,1)
J_2t=np.array([1,0,a*sym.cos(theta)-b*sym.sin(theta)]).reshape(3,1)
J_2=np.hstack((J_2n,J_2t))
C2=contact_point_symbolic_2D(mysystem,phi=phi_2,psi=psi_2,J=J_2,friction=0.5,name="contact point")


# Contact Point 3: Left finger
phi_3 = (x-x_1)*sym.cos(theta) - (y_1-y)*sym.sin(theta)- b 
psi_3 = (y_1-y)*sym.cos(theta) - (x-x_1)*sym.sin(theta)
J_3n=np.array([sym.cos(theta),sym.sin(theta),psi_3]).reshape(3,1)
J_3t=np.array([-sym.sin(theta),sym.cos(theta),-b]).reshape(3,1)
J_3=np.hstack((J_3n,J_3t))
C3=contact_point_symbolic_2D(mysystem,phi=phi_3,psi=psi_3,J=J_3,friction=0.5,name="contact point")

# Contact Point 4: Right finger
phi_4 = (x_2-x)*sym.cos(theta) + (y_2-y)*sym.sin(theta)- b 
psi_4 = (y_2-y)*sym.cos(theta) - (x-x_2)*sym.sin(theta)
J_4n=np.array([-sym.cos(theta),-sym.sin(theta),psi_4]).reshape(3,1)
J_4t=np.array([-sym.sin(theta),sym.cos(theta),b]).reshape(3,1)
J_4=np.hstack((J_4n,J_4t))
C4=contact_point_symbolic_2D(mysystem,phi=phi_4,psi=psi_4,J=J_4,friction=0.5,name="contact point")

mysystem.build_and_linearize()


Eta=environment(0)
Eta.dict_of_values={x:0,y:a,theta:0,x_1:-b,x_2:b,y_1:a,y_2:a,
     mysystem.v_o[0]:0,mysystem.v_o[1]:0,mysystem.v_o[2]:0,
     mysystem.u_lambda[0]:0,mysystem.u_lambda[1]:0,mysystem.u_lambda[2]:0,mysystem.u_lambda[3]:0,
     mysystem.u_lambda[4]:0,mysystem.u_lambda[5]:0,mysystem.u_lambda[6]:0,mysystem.u_lambda[7]:0,
     mysystem.u_m[0]:0,mysystem.u_m[1]:0,mysystem.u_m[2]:0,mysystem.u_m[3]:0,
     mysystem.h:0.1}

epsilon=np.array([2,2,2,1,1,1,1,5,5,5,1,1,1,1,50,50,50,50,50,50,50,50]).reshape(22,1)
sys=system_numeric(mysystem)
sys.add_environment(Eta,epsilon)

up_shift=0.4
#x_goal=np.array([0,y_goal,0.0,1,y_goal,-1,y_goal,0,0,0]).reshape(10,1)
x_goal=np.array([0,a+up_shift,0.0,-b,a-up_shift,b,a-up_shift,0,0,0]).reshape(10,1)
sys.goal=zonotope(x_goal.reshape(10,1),np.diag([0,0,0,1,1,1,1,0,0,0]))
x0=np.array([0,a,0.0,-b,a,b,a,0,0,0]).reshape(10,1)
T=15
x,u,u_lambda,x_tishcom,x_time=point_trajectory_tishcom(sys,x0,[sys.goal],T,optimize_controls_indices=[0,1,2,3])


"""
Visualization
"""
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

def animate(X,ax1):
    x,y,theta,x_1,y_1,x_2,y_2=X[0:7]
    a=1
    b=1
    # x
    x_left=x-b*np.cos(theta)
    x_left_bar=x-b*np.cos(theta)-b*np.sin(theta)
    x_right=x+b*np.cos(theta)
    x_right_bar=x+b*np.cos(theta)-b*np.sin(theta)
    # y
    y_left=y-a-a*np.sin(theta)
    y_left_bar=y_left+2*a*np.cos(theta)
    y_right=y-a+a*np.sin(theta)
    y_right_bar=y_right+2*a*np.cos(theta)
    # Good now plot
    ax1.set_xlabel("x",fontsize=20)
    ax1.set_ylabel("y",fontsize=20)
    ax1.set_xlim([-3*a,3*a])
    ax1.set_ylim([-0.7,4])
    fig.gca().set_aspect('equal')
    bar=[patches.Polygon(np.array([[x_left,x_left_bar,x_right_bar,x_right],[y_left,y_left_bar,y_right_bar,y_right]]).reshape(2,4).T, True)]
    ax1.add_collection(PatchCollection(bar,color=(0.8,0.3,0.4),alpha=0.8,edgecolor=(0,0,0)))
    ax1.plot([x_1,x_2],[y_1,y_2],'*',linewidth=3,markersize=10)
    ax1.plot([-3,3],[0,0],linewidth=3)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.5)
    
for t in range(T+1):
    fig,ax = plt.subplots()
#    ax.set_title("Balancing t=%d \n Right finger: %s \n Left Finger: %s"%(t,mode_name[list_of_modes[t][1]],mode_name[list_of_modes[t][2]]))
    animate(x[t],ax)
    if t==-1:
        ax.arrow(-1.8, 0.1, 0.3, 0.0, head_width=0.3, head_length=0.3, linewidth=10, fc='k', ec='k')
    fig.savefig('figures/box_%d.png'%t, dpi=100)