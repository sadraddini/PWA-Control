#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 16:56:05 2019

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
C1=contact_point_symbolic_2D(mysystem,phi=phi_1,psi=psi_1,J=J_1,friction=0.8,name="ground contact")
C1.sliding=False


"""
Finger Right: surface contact with curvy side
"""
_zeta=np.dot(Rotate(np.pi/2-theta_m),np.array([x_R-x_m,y_R-y_m]))
phi_2=s_m-_zeta[0]
psi_m=_zeta[1]
psi_o=R * (theta_m-theta)
psi_2=psi_o-psi_m
phi_2+=ReLU(-psi_m,100) # ReLU added to punish negative -psi_m
J_2n=np.array([-sym.sin(theta_m),sym.cos(theta_m),-(psi_m-p*sym.sin(theta_m-theta))]).reshape(3,1)
J_2t=np.array([sym.cos(theta_m),sym.sin(theta_m),R+phi_2-p*sym.cos(theta_m-theta)]).reshape(3,1)
J_2=np.hstack((J_2n,J_2t))
C2=contact_point_symbolic_2D(mysystem,phi=phi_2,psi=psi_2,J=J_2,friction=0.3,name="right finger")
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
C3=contact_point_symbolic_2D(mysystem,phi=phi_3,psi=psi_3,J=J_3,friction=0.3,name="left finger")
#C3=contact_point_symbolic_2D(mysystem,phi=phi_2,psi=psi_2,J=J_2,friction=0.5,name="right finger")
C3.no_contact=False


# Build the Symbolic MLD system
mysystem.build_and_linearize()
"""
Now the sequences
"""

N=30
delta_theta=0.05


x_traj={}
u_traj={}
lambda_traj={}
# Zero Env
h=0.03
sys=system_numeric(mysystem)
Eta_0=environment(0)
Eta_0.dict_of_values={x:0,y:R-p,theta:delta_theta/2,x_m:0.18,y_m:0.82,theta_m:np.pi/2.0-0.4,s_m:0.8,
     mysystem.v_o[0]:0,mysystem.v_o[1]:0,mysystem.v_o[2]:0,
     mysystem.u_lambda[0]:0,mysystem.u_lambda[1]:0,mysystem.u_lambda[2]:0,mysystem.u_lambda[3]:0,
     mysystem.u_lambda[4]:0,mysystem.u_lambda[5]:0,
     mysystem.u_m[0]:0,mysystem.u_m[1]:0,mysystem.u_m[2]:0,mysystem.u_m[3]:0,
     mysystem.h:0.03}
epsilon_max=np.array([20,20,0.7,10,10,10,10,50,50,50,2,2,2,2,500,500,500,500,70,70]).reshape(20,1)
epsilon_min=-np.array([20,20,0.7,10,10,10,10,50,50,50,2,2,2,2,500,500,500,500,70,70]).reshape(20,1)
sys.add_environment(Eta_0,epsilon_max,epsilon_min)

theta_shift=delta_theta
x_goal=np.array([0,2,theta_shift,0.12,0.82,np.pi/2-0.4,0.8,0,0,0]).reshape(10,1)
x0=np.array([0,R-p,0,0.18,0.82,np.pi/2-0.4,0.8,0,0,0]).reshape(10,1).reshape(10,1)
sys.goal=zonotope(x_goal.reshape(10,1),100*np.diag([1,1,0,1,1,1,1,0,0,0]))
T=5
x_traj[0],u_traj[0],lambda_traj[0],mode=point_trajectory_tishcom(sys,x0,[sys.goal],T,optimize_controls_indices=[0,1,2,3],cost=1)


for k in range(N):
    sys=system_numeric(mysystem)
    Eta_1=environment(0)
    Eta_1.dict_of_values={x:x_traj[k][T][0],y:x_traj[k][T][1],theta:x_traj[k][T][2]+delta_theta/2,
                          x_m:x_traj[k][T][3],y_m:x_traj[k][T][4],theta_m:x_traj[k][T][5],s_m:x_traj[k][T][5],
         mysystem.v_o[0]:x_traj[k][T][6],mysystem.v_o[1]:x_traj[k][T][8],mysystem.v_o[2]:x_traj[k][T][9],
#         mysystem.u_lambda[0]:lambda_traj[k][T-1][0],mysystem.u_lambda[1]:lambda_traj[k][T-1][1],
#         mysystem.u_lambda[2]:lambda_traj[k][T-1][2],mysystem.u_lambda[3]:lambda_traj[k][T-1][3],
#         mysystem.u_lambda[4]:lambda_traj[k][T-1][4],mysystem.u_lambda[5]:lambda_traj[k][T-1][5],
         mysystem.u_lambda[0]:0,mysystem.u_lambda[1]:0,
         mysystem.u_lambda[2]:0,mysystem.u_lambda[3]:0,
         mysystem.u_lambda[4]:0,mysystem.u_lambda[5]:0,
#         mysystem.u_m[0]:u_traj[k][T-1][0],mysystem.u_m[1]:u_traj[k][T-1][1],
#         mysystem.u_m[2]:u_traj[k][T-1][2],mysystem.u_m[3]:u_traj[k][T-1][3],
         mysystem.u_m[0]:0,mysystem.u_m[1]:0,
         mysystem.u_m[2]:0,mysystem.u_m[3]:0,
         mysystem.h:h}
    epsilon_max=np.array([20,20,3.7,10,10,10,10,50,50,50,1,1,1,1,500,500,500,500,70,70,70,70]).reshape(22,1)
    epsilon_min=-np.array([20,20,0.7,10,10,10,10,50,50,50,1,1,1,1,500,500,500,500,70,70,70,70]).reshape(22,1)
    sys.add_environment(Eta_1,epsilon_max,epsilon_min)
    
    theta_shift=delta_theta*(k+2)
    x_goal=np.array([0,2,theta_shift,0.18,0.82,np.pi/2-0.3,0.8,0,0,0]).reshape(10,1)
    x0=x_traj[k][T].reshape(mysystem.n,1)
    sys.goal=zonotope(x_goal.reshape(10,1),10*np.diag([1,1,0,1,1,1,1,0,0,0]))
    T=5
    x_traj[k+1],u_traj[k+1],lambda_traj[k+1],mode=point_trajectory_tishcom(sys,x0,[sys.goal],T,optimize_controls_indices=[0,1,2,3],cost=1)



x_traj=merge_timed_vectors_glue([x_traj[k] for k in range(k+1)])
u_lambda=merge_timed_vectors_skip([lambda_traj[k] for k in range(k+1)])
T=max(x_traj.keys())



"""
Visualization
"""
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

def carrot_vertices(x,y,theta,N=20):
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
        
def generate_figures():    
    for t in range(T+1):
        fig,ax = plt.subplots()
        fig.set_size_inches(10, 10)
        if t!=T:
            ax.set_title(r'%d/%d %s'%(t,T,mysystem.name)+'\n'+
                         r'$\lambda^{ground}_n: %0.1f * \lambda^{ground}_t: %0.1f $'
                         %(u_lambda[t][0],u_lambda[t][1])
                         +'\n'+r'$\lambda^{rF}_n: %0.1f * \lambda^{rF}_t: %0.1f $'
                         %(u_lambda[t][2],u_lambda[t][3])
                         +'\n'+r'$\lambda^{lF}_n: %0.1f * \lambda^{lF}_t: %0.1f $'
                         %(u_lambda[t][4],u_lambda[t][5]))
        animate(x_traj[t],ax,fig)
        if t==-1:
            ax.arrow(-1.8, 0.1, 0.3, 0.0, head_width=0.3, head_length=0.3, linewidth=10, fc='k', ec='k')
        fig.savefig('figures/carrot_%d.png'%t, dpi=100)

def single_figure(X):
    fig,ax = plt.subplots()
    fig.set_size_inches(10, 10)
    """
    Numericals
    """
    config={x:X[0],y:X[1],theta:X[2],x_m:X[3],y_m:X[4],theta_m:X[5],s_m:X[6]}
    
    phi_1_n=phi_1.Evaluate(config)
    phi_2_n=phi_2.Evaluate(config)
    phi_3_n=phi_3.Evaluate(config)

    psi_1_n=psi_1.Evaluate(config)
    psi_2_n=psi_2.Evaluate(config)
    psi_3_n=psi_3.Evaluate(config)
    psi_m_n=psi_m.Evaluate(config)
    
    print "phi_1:",phi_1_n
    print "phi_2:",phi_2_n
    print "phi_3:",phi_3_n
    
    print "psi_1:",psi_1_n
    print "psi_2:",psi_2_n
    print "psi_3:",psi_3_n
    print "psi_m(3):",psi_m_n
 
#    x_R_n=x_R.Evaluate(config)
#    y_R_n=y_R.Evaluate(config)
#    
#    x_c_n=x_c.Evaluate(config)
#    y_c_n=y_c.Evaluate(config)
    
    animate(X,ax,fig)
    
generate_figures()
    
    