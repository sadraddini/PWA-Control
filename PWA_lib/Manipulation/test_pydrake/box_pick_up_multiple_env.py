#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 12:53:04 2019

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
from PWA_lib.Manipulation.system_numeric import environment,merge_timed_vectors_glue,merge_timed_vectors_skip
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
psi_1= - (x + a * sym.sin(theta) - b * sym.cos(theta) )
J_1n=np.array([0,1,-b*sym.cos(theta)-a*sym.sin(theta)]).reshape(3,1)
J_1t=np.array([1,0,a*sym.cos(theta)+b*sym.sin(theta)]).reshape(3,1)
J_1=np.hstack((J_1n,J_1t))
C1=contact_point_symbolic_2D(mysystem,phi=phi_1,psi=psi_1,J=J_1,friction=0.8,name="contact point")
C1.sliding=False
#C1.sticking=False


# Contact Point 2: The ground with right corner
phi_2= y - a * sym.cos(theta) + b * sym.sin(theta)
psi_2= - ( x + a * sym.sin(theta) + b * sym.cos(theta) )
J_2n=np.array([0,1,b*sym.cos(theta)-a*sym.sin(theta)]).reshape(3,1)
J_2t=np.array([1,0,a*sym.cos(theta)-b*sym.sin(theta)]).reshape(3,1)
J_2=np.hstack((J_2n,J_2t))
C2=contact_point_symbolic_2D(mysystem,phi=phi_2,psi=psi_2,J=J_2,friction=0.8,name="contact point")
C2.sliding=False


# Contact Point 3: Left finger
phi_3 = (x-x_1)*sym.cos(theta) - (y_1-y)*sym.sin(theta)- b 
psi_3 = (y_1-y)*sym.cos(theta) + (x-x_1)*sym.sin(theta)
J_3n=np.array([sym.cos(theta),sym.sin(theta),psi_3]).reshape(3,1)
J_3t=np.array([-sym.sin(theta),sym.cos(theta),-b]).reshape(3,1)
J_3=np.hstack((J_3n,J_3t))
C3=contact_point_symbolic_2D(mysystem,phi=phi_3,psi=psi_3,J=J_3,friction=0.5,name="contact point")
C3.sliding=False
C3.sticking=False

# Contact Point 4: Right finger
phi_4 = (x_2-x)*sym.cos(theta) + (y_2-y)*sym.sin(theta)- b 
psi_4 = (y_2-y)*sym.cos(theta) - (x_2-x)*sym.sin(theta)
J_4n=np.array([-sym.cos(theta),-sym.sin(theta),psi_4]).reshape(3,1)
J_4t=np.array([-sym.sin(theta),sym.cos(theta),b]).reshape(3,1)
J_4=np.hstack((J_4n,J_4t))
C4=contact_point_symbolic_2D(mysystem,phi=phi_4,psi=psi_4,J=J_4,friction=0.5,name="contact point")
C4.sliding=False

# Build the Symbolic MLD system
mysystem.build_and_linearize()
sys=system_numeric(mysystem)

# Add a numerical environment
Eta_0=environment(0)
Eta_0.dict_of_values={x:0,y:a,theta:0,x_1:-b,x_2:b,y_1:a,y_2:a,
     mysystem.v_o[0]:0,mysystem.v_o[1]:0,mysystem.v_o[2]:0,
     mysystem.u_lambda[0]:0,mysystem.u_lambda[1]:0,mysystem.u_lambda[2]:0,mysystem.u_lambda[3]:0,
     mysystem.u_lambda[4]:0,mysystem.u_lambda[5]:0,mysystem.u_lambda[6]:0,mysystem.u_lambda[7]:0,
     mysystem.u_m[0]:0,mysystem.u_m[1]:0,mysystem.u_m[2]:0,mysystem.u_m[3]:0,
     mysystem.h:0.1}
epsilon_max=np.array([20,20,0.22,10,10,10,10,50,50,50,1,1,1,1,500,500,500,500,70,70,70,70]).reshape(22,1)
epsilon_min=-np.array([20,20,0.3,10,10,10,10,50,50,50,1,1,1,1,500,500,500,500,70,70,70,70]).reshape(22,1)
sys.add_environment(Eta_0,epsilon_max,epsilon_min)


# Add a numerical environment
Eta_1=environment(1)
Eta_1.dict_of_values={x:0,y:a,theta:0.3,x_1:-b,x_2:b,y_1:a,y_2:a,
     mysystem.v_o[0]:0,mysystem.v_o[1]:0,mysystem.v_o[2]:0,
     mysystem.u_lambda[0]:0,mysystem.u_lambda[1]:0,mysystem.u_lambda[2]:0,mysystem.u_lambda[3]:0,
     mysystem.u_lambda[4]:0,mysystem.u_lambda[5]:0,mysystem.u_lambda[6]:0,mysystem.u_lambda[7]:0,
     mysystem.u_m[0]:0,mysystem.u_m[1]:0,mysystem.u_m[2]:0,mysystem.u_m[3]:0,
     mysystem.h:0.1}
epsilon_max=np.array([20,20,0.5,10,10,10,10,50,50,50,1,1,1,1,500,500,500,500,70,70,70,70]).reshape(22,1)
epsilon_min=-np.array([20,20,-0.2,10,10,10,10,50,50,50,1,1,1,1,500,500,500,500,70,70,70,70]).reshape(22,1)
sys.add_environment(Eta_1,epsilon_max,epsilon_min)
    
if True:
    up_shift=0
    right_shift=0
    theta_shift=0.4
    x_goal=np.array([right_shift,a+up_shift,theta_shift,-b+right_shift,a+up_shift,b+right_shift,a+up_shift,0,0,0]).reshape(10,1)
    x0=np.array([0,a,0,-b,a,b,a,0,0,0]).reshape(10,1)
    #x0=x_goal
    sys.goal=zonotope(x_goal.reshape(10,1),100*np.diag([1,1,0,1,1,1,1,1,1,1]))
    T=10
    x_traj,u,u_lambda,x_tishcom,x_time=point_trajectory_tishcom(sys,x0,[sys.goal],T,optimize_controls_indices=[0,1,2,3],cost=1)





"""
Visualization
"""
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

def animate(X,ax1,fig):
    x,y,theta,x_1,y_1,x_2,y_2=X[0:7]
    a=1
    b=1
    R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    # Good now plot
    ax1.set_xlabel("x",fontsize=20)
    ax1.set_ylabel("y",fontsize=20)
    ax1.set_xlim([-3*a,3*a])
    ax1.set_ylim([-0.7,3])
    fig.gca().set_aspect('equal')
    # Bar
    up_right=np.array([x,y])+np.dot(R,[b,a])
    down_right=np.array([x,y])+np.dot(R,[b,-a])
    down_left=np.array([x,y])+np.dot(R,[-b,-a])
    up_left=np.array([x,y])+np.dot(R,[-b,a])
    bar=[patches.Polygon(np.array([[up_right[0],down_right[0],down_left[0],up_left[0]],[up_right[1],down_right[1],down_left[1],up_left[1]]]).reshape(2,4).T, True)]
    ax1.add_collection(PatchCollection(bar,color=(0.8,0.3,0.4),alpha=0.8,edgecolor=(0,0,0)))
    # Finger left
    L,h=0.1,0.3
    base_up,base_down,corner=np.array([x_1,y_1])+np.dot(R,[-h,-L]),np.array([x_1,y_1])+np.dot(R,[-h,L]),np.array([x_1,y_1])
    finger_left=[patches.Polygon(np.array([[base_up[0],base_down[0],corner[0]],[base_up[1],base_down[1],corner[1]]]).reshape(2,3).T, True)]
    ax1.add_collection(PatchCollection(finger_left,color=(0.2,0.2,0.8),alpha=0.8,edgecolor=(0,0,0)))
    # Finger right
    L,h=0.1,0.3
    base_up,base_down,corner=np.array([x_2,y_2])+np.dot(R,[h,-L]),np.array([x_2,y_2])+np.dot(R,[h,L]),np.array([x_2,y_2])
    finger_right=[patches.Polygon(np.array([[base_up[0],base_down[0],corner[0]],[base_up[1],base_down[1],corner[1]]]).reshape(2,3).T, True)]
    ax1.add_collection(PatchCollection(finger_right,color=(0.2,0.2,0.8),alpha=0.8,edgecolor=(0,0,0)))    
    # Star 
    ax1.plot([x_1,x_2],[y_1,y_2],'*',linewidth=3,markersize=10)
    ax1.plot([-3,3],[0,0],linewidth=3,color=(0.1,0.1,0.1))
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.5)
    # The contact forces
#    ax1.quiver([x-b*np.cos(theta)+a*np.sin(theta),x-b*np.cos(theta)+a*np.sin(theta)],[y-a*np.cos(theta)-b*np.sin(theta),y-a*np.cos(theta)-b*np.sin(theta)+U[0]],scale=10)
    
    
def generate_figures(): 
    for t in range(T+1):
        fig,ax = plt.subplots()
        fig.set_size_inches(10, 8)
        if t!=T:
            ax.set_title(r'%d/%d Box'%(t,T)+'\n'+
                         r'$\lambda^{lc}_n: %0.1f * \lambda^{rc}_n: %0.1f * \lambda^{lf}_n: %0.1f * \lambda^{rf}_n: %0.1f$'
                         %(u_lambda[t][0],u_lambda[t][2],u_lambda[t][4],u_lambda[t][6])
                         +'\n'+r'$\lambda^{lc}_t: %0.1f * \lambda^{rc}_t: %0.1f * \lambda^{lf}_t: %0.1f * \lambda^{rf}_t: %0.1f $'
                         %(u_lambda[t][1],u_lambda[t][3],u_lambda[t][5],u_lambda[t][7]))
        animate(x_traj[t],ax,fig)
        if t==-1:
            ax.arrow(-1.8, 0.1, 0.3, 0.0, head_width=0.3, head_length=0.3, linewidth=10, fc='k', ec='k')
        fig.savefig('figures/box_%d.png'%t, dpi=100)

generate_figures()