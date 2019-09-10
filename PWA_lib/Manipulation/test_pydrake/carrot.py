#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 16:37:42 2019

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
from PWA_lib.Manipulation.system_numeric import environment,merge_timed_vectors_glue,merge_timed_vectors_skip,\
    trajectory_to_list_of_linear_cells,trajectory_to_list_of_linear_cells_full_linearization,\
    PWA_cells_from_state,hybrid_reachable_sets_from_state,environment_from_state
from PWA_lib.trajectory.poly_trajectory import point_trajectory_tishcom,\
    polytopic_trajectory_given_modes,point_trajectory_given_modes_and_central_traj,\
    point_trajectory_tishcom_time_varying


mysystem=system_symbolic("Carrot", 2)
x,y,theta=sym.Variable("x"),sym.Variable("y"),sym.Variable("theta")
x_m,y_m,theta_m,s_m=sym.Variable("x_m"),sym.Variable("y_m"),sym.Variable("theta_m"),sym.Variable("s_m")
mysystem.q_o=np.array([x,y,theta])
mysystem.q_m=np.array([x_m,y_m,theta_m,s_m])


def Rotate(t):
    return np.array([[sym.cos(t),-sym.sin(t)],[sym.sin(t),sym.cos(t)]])

def ReLU(cond,K):
    return sym.log(1+sym.exp(K*cond))
    
    
R=1 # The radius of the carrot
p=4*R/(3*np.pi)
# Dynamics:
mysystem.C=np.zeros((3,3))
g=9.8
mysystem.tau_g=np.array([0,-g,0])
mysystem.B=np.zeros((3,0))

# Mass Matrix
M=1
I=M*R**2*(1/2.0-16/(9*np.pi**2))
h=0.05
mysystem.M=blk(*[M,M,I])
mysystem.M_inv=np.linalg.inv(mysystem.M)

x_c=x-p*sym.sin(theta)
y_c=y+p*sym.cos(theta)

x_R=x_c+R*sym.sin(theta_m)
y_R=y_c-R*sym.cos(theta_m)
"""
Ground
"""
phi_1= y_c-R
psi_1= - R * theta - x_c 
J_1n=np.array([0,1,-p*sym.sin(theta)]).reshape(3,1)
J_1t=np.array([1,0,R-p*sym.cos(theta)]).reshape(3,1)
J_1=np.hstack((J_1n,J_1t))
C1=contact_point_symbolic_2D(mysystem,phi=phi_1,psi=psi_1,J=J_1,friction=0.75,name="ground contact")
C1.sliding=False

"""
Finger Right: surface contact with curvy side
"""
_zeta=np.dot(Rotate(np.pi/2-theta_m),np.array([x_R-x_m,y_R-y_m]))
phi_2=s_m-_zeta[0]
psi_m=_zeta[1]
psi_o=R * (theta_m-theta)
psi_2=psi_o-psi_m
J_2n=np.array([-sym.sin(theta_m),sym.cos(theta_m),-(psi_m-p*sym.sin(theta_m-theta))]).reshape(3,1)
J_2t=np.array([sym.cos(theta_m),sym.sin(theta_m),R+phi_2-p*sym.cos(theta_m-theta)]).reshape(3,1)
J_2=np.hstack((J_2n,J_2t))
C2=contact_point_symbolic_2D(mysystem,phi=phi_2,psi=psi_2,J=J_2,friction=0.5,name="right finger")
#C2.sliding=False
C2.no_contact=False

#raise 1


"""
Finger left: point contact with flat side
"""
_xi_1=np.dot(Rotate(np.pi/2+theta_m),np.array([s_m,0]))+np.array([x_m-x_c,y_m-y_c])
_xi_2=np.dot(Rotate(-np.pi-theta),_xi_1)
phi_3=-_xi_2[1]
psi_3=_xi_2[0]
#phi_3+=ReLU(psi_3-R,100) # RELU ADDED TO FIGURE psi>1 
J_3n=np.array([sym.sin(theta),-sym.cos(theta),-psi_3]).reshape(3,1)
J_3t=np.array([-sym.cos(theta),-sym.sin(theta),phi_3+p]).reshape(3,1)
J_3=np.hstack((J_3n,J_3t))
C3=contact_point_symbolic_2D(mysystem,phi=phi_3,psi=psi_3,J=J_3,friction=0.5,name="left finger")
#C3=contact_point_symbolic_2D(mysystem,phi=phi_2,psi=psi_2,J=J_2,friction=0.5,name="right finger")
C3.no_contact=False
C3.sliding=False

# Build the Symbolic MLD system
mysystem.build_and_linearize()
sys=system_numeric(mysystem)

x_m_0=0.11
y_m_0=0.85
theta_m_0=1.39
s_m_0=0.85

Eta_1=environment(0)
Eta_1.dict_of_values={x:0,y:R-p,theta:0.125,x_m:x_m_0,y_m:y_m_0,theta_m:theta_m_0,s_m:s_m_0,
     mysystem.v_o[0]:0,mysystem.v_o[1]:0,mysystem.v_o[2]:0,
     mysystem.u_lambda[0]:0,mysystem.u_lambda[1]:0,mysystem.u_lambda[2]:0,mysystem.u_lambda[3]:0,
     mysystem.u_lambda[4]:0,mysystem.u_lambda[5]:0,
     mysystem.u_m[0]:0,mysystem.u_m[1]:0,mysystem.u_m[2]:0,mysystem.u_m[3]:0,
     mysystem.h:h}
    
epsilon_max=np.array([2,2,1,10,10,10,10,50,50,50,1,1,1,1,500,500,500,500,500,500]).reshape(20,1)
epsilon_min=-np.array([2,2,0,10,10,10,10,50,50,50,1,1,1,1,500,500,500,500,500,500]).reshape(20,1)
sys.add_environment(Eta_1,epsilon_max,epsilon_min)
#sys.add_environment(Eta_1)


#raise 1

up_shift=0
right_shift=0
theta_shift=0.25
x_goal=np.array([0,R-p,theta_shift,x_m_0,y_m_0,theta_m_0,s_m_0,0,0,0]).reshape(10,1)
x0=np.array([0,R-p,0,x_m_0,y_m_0,theta_m_0,s_m_0,0,0,0]).reshape(10,1).reshape(10,1)
sys.goal=zonotope(x_goal.reshape(10,1),100*np.diag([1,1,0,1,1,1,1,0,0,0]))
T=10
sys.scale=np.array([0,0,0,1,1,1,1,0,0,0])
x_traj,u_traj,u_lambda,mode=point_trajectory_tishcom(sys,x0,[sys.goal],T,optimize_controls_indices=[0,1,2,3],cost=1,eps=0.0)

list_of_linear_cells=trajectory_to_list_of_linear_cells(sys,Eta_1,x_traj,u_traj,u_lambda,mode)
list_of_linear_cells_full=trajectory_to_list_of_linear_cells_full_linearization(mysystem,x_traj,u_traj,u_lambda,mode,h,epsilon_min,epsilon_max)
list_of_linear_cells_PWA=PWA_cells_from_state(mysystem,x_traj[2],h,epsilon_min,epsilon_max)
list_of_reachable_sets=hybrid_reachable_sets_from_state(mysystem,x_traj[1],h,epsilon_min,epsilon_max)
# The time varying PWA system:
list_of_env=[environment_from_state(mysystem,x_traj[t],h) for t in range(T)]
sys=system_numeric(mysystem)
for env in list_of_env:
    sys.add_environment(env,epsilon_max,epsilon_min)
sys.scale=np.array([0,0,0,1,1,1,1,0,0,0])
sys.goal=zonotope(x_goal.reshape(10,1),100*np.diag([1,1,0,1,1,1,1,0,0,0]))
x_TV,u_TV,lambda_TV,mode_TV=point_trajectory_tishcom_time_varying(sys,list_of_env,x0,[sys.goal],T,optimize_controls_indices=[0,1,2,3],cost=1)
# Check the linear one
x,u=point_trajectory_given_modes_and_central_traj(x_traj,list_of_linear_cells_full,sys.goal)
#raise 1
#x,u,G,theta=polytopic_trajectory_given_modes(x0,list_of_linear_cells_full,sys.goal,eps=0,order=1,scale=[])
#raise 1

"""
Visualization
"""
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

def carrot_vertices(x,y,theta,N=50):
    x_c,y_c=x-p*np.sin(theta),y+p*np.cos(theta)
    v=np.empty((N,2))
    for k in range(N):
        phi=-np.pi/2+np.pi/(N-1)*k
        v[k,0]=x_c+R*np.sin(phi+theta)
        v[k,1]=y_c-R*np.cos(phi+theta)
    return v

    
def animate(X,ax,fig):
    fig.gca().set_aspect('equal')
    ax.set_xlabel("x",fontsize=20)
    ax.set_ylabel("y",fontsize=20)
    ax.set_xlim([-4*R,4*R])
    ax.set_ylim([-2,4*R])
    x,y,theta,x_m,y_m,theta_m,s_m=X[0:7]
    a,b=0.4,3
    ax.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    # Carrot First
    v=carrot_vertices(x,y,theta)
    ax.add_collection(PatchCollection([patches.Polygon(v, True)],color=(1,0,0),alpha=0.31,edgecolor=(1,0,0)))
    # Center of Mass
    ax.plot([x],[y],'+',color=(1,0,0))
    # Ground
    ax.plot([-12,12],[0,0],'black',linewidth=2) 
    # Fingers
    R_m=np.array([[np.cos(theta_m),-np.sin(theta_m)],[np.sin(theta_m),np.cos(theta_m)]]).reshape(2,2)
    RF_nominal=np.array([[0,0,b,b],[-s_m,-s_m-a,-s_m-a,-s_m]])
    LF_nominal=np.array([[0,0,b,b],[s_m,s_m+a,s_m+a,s_m]])
    RF=np.dot( R_m, RF_nominal ) + np.array([x_m,y_m]).reshape(2,1)
    LF=np.dot( R_m, LF_nominal ) + np.array([x_m,y_m]).reshape(2,1)
    ax.add_collection(PatchCollection([patches.Polygon(RF.T, True)],color=(0.2,0.2,0.2),alpha=0.91,edgecolor=(0,0,0)))
    ax.add_collection(PatchCollection([patches.Polygon(LF.T, True)],color=(0.2,0.2,0.2),alpha=0.91,edgecolor=(0,0,0)))
        
def generate_figures(trajectory_of_x,trajectory_of_lambda):    
    for t in range(T+1):
        fig,ax = plt.subplots()
        fig.set_size_inches(10, 10)
        if t!=T:
            ax.set_title(r'%d/%d %s'%(t,T,mysystem.name)+'\n'+
                         r'$\lambda^{ground}_n: %0.1f * \lambda^{ground}_t: %0.1f $'
                         %(trajectory_of_lambda[t][0],trajectory_of_lambda[t][1])
                         +'\n'+r'$\lambda^{rF}_n: %0.1f * \lambda^{rF}_t: %0.1f $'
                         %(trajectory_of_lambda[t][2],trajectory_of_lambda[t][3])
                         +'\n'+r'$\lambda^{lF}_n: %0.1f * \lambda^{lF}_t: %0.1f $'
                         %(trajectory_of_lambda[t][4],trajectory_of_lambda[t][5]))
        animate(trajectory_of_x[t],ax,fig)
        if t==-1:
            ax.arrow(-1.8, 0.1, 0.3, 0.0, head_width=0.3, head_length=0.3, linewidth=10, fc='k', ec='k')
        fig.savefig('figures/carrot_%d.png'%t, dpi=100)

generate_figures(x_TV,lambda_TV)
    
"""
Numericals
"""
#config={x:0,y:R-p,theta:0,x_m:0.18,y_m:0.82,theta_m:np.pi/2-0.3,s_m:0.8}
#X=np.array([x,y,theta,x_m,y_m,theta_m,s_m])
#X_n=sym.Evaluate(X,config)
#
#phi_1_n=phi_1.Evaluate(config)
#phi_2_n=phi_2.Evaluate(config)
#phi_3_n=phi_3.Evaluate(config)
#
#x_R_n=x_R.Evaluate(config)
#y_R_n=y_R.Evaluate(config)
#
#x_c_n=x_c.Evaluate(config)
#y_c_n=y_c.Evaluate(config)
#
#psi_1_n=psi_1.Evaluate(config)
#psi_2_n=psi_2.Evaluate(config)
#psi_3_n=psi_3.Evaluate(config)
#psi_m_n=psi_m.Evaluate(config)    