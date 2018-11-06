#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 06:49:55 2018

@author: sadraddini
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection


R=3 # Carrot radius
g=9.8 # Gravity
mu_ground=0.4
mu_finger=0.6

K_ground=400
c_ground=150
K_finger=700
c_finger=10

K_torsion=5000
c_torsion=100
K_torsion_ground=1000
c_torsion_ground=100

t=0
global t


h=0.01 # Time step
p=4*R/(3*np.pi) # Center of mass distance from O
I=40 # Moment of Inerita

def evolve_carrot(X,U,gripper):
    x,y,theta,x_dot,y_dot,theta_dot=X[0:6]
    X_new=np.empty(10)
    X_new[0]=x+x_dot*h
    X_new[1]=y+y_dot*h
    X_new[2]=theta+theta_dot*h
    f=carrot_forces(X,U,gripper)
    f_x=f[0]
    f_y=f[1]
    f_theta=f[2] # momentum around CoM
    X_new[3]=x_dot+h*f_x
    X_new[4]=y_dot+h*f_y
    X_new[5]=theta_dot+h*f_theta/I
    # Fingers
    c,d,D,psi=X[6:10]
    v_c,v_d,v_D,v_psi=U[0:4]
    X_new[6]=c+h*v_c
    X_new[7]=d+h*v_d
    X_new[8]=D+h*v_D
    X_new[9]=psi+h*v_psi
    return X_new


def carrot_forces(X,U,gripper=1):
    f_g=np.array([0,-g,0])
    f_ground=carrot_forces_ground(X)
    f_left=left_finger_force(X,U)
    f_right=right_finger_force(X,U)
    print "ground force:",f_ground
    print "left finger force:",f_left
    print "right finger force:",f_right
    if gripper==1:
        return f_ground+f_g+f_left+f_right
    elif gripper==0:
        return f_ground+f_g
    

def carrot_forces_ground(X):
    x,y,theta,x_dot,y_dot,theta_dot=X[0:6]
    if np.cos(theta)>=0: # theta between -pi/2 and pi/2, cylinder down
        c="down"
        d=y-R
        v_d=y_dot
        v_contact=x_dot+R*theta_dot
    elif np.sin(theta)>0: # left half touching
        c="left"
        d=y-R*np.sin(theta)
        v_d=y_dot-R*theta_dot*np.cos(theta)
        v_contact=x_dot+R*theta_dot*np.sin(theta)
    else: # right
        c="right"
        d=y+R*np.sin(theta)
        v_d=y_dot+R*theta_dot*np.cos(theta)
        v_contact=x_dot-R*theta_dot*np.sin(theta)
    f=np.zeros(3)
    print "ground penetration",d, "mode:",c
    if d>=0:
        return f
    elif d<0:
        fy=-K_ground*d-c_ground*v_d
        fx=-fy*mu_ground*np.sign(v_contact)
        if c=="down":
            f_theta=-fy*p*np.sin(theta)+fx*(R-p*np.cos(theta))
        elif c=="left":
            f_theta=fy*(-R*np.cos(theta)-p*np.sin(theta))+fx*(R*np.sin(theta)-p*np.cos(theta))
        else:
            f_theta=-fy*(-R*np.cos(theta)+p*np.sin(theta))+fx*(-R*np.sin(theta)-p*np.cos(theta))
#        if y-R*np.sin(theta)<0 and y+R*np.sin(theta)<0:
        if c=="left" or c=="right": # Add Torsion Spring
            f_theta+=-c_torsion_ground*theta_dot-K_torsion_ground*(theta-np.pi)
            fx+=-c_ground*x_dot
#            f_theta=-I*theta_dot/h
        return np.array([fx,fy,f_theta])
        
def left_finger_force(X,U):
    x,y,theta,x_dot,y_dot,theta_dot=X[0:6]
    psi=X[9]
    c,d=X[6:8]
    v_c,v_d=U[0:2]
    # d is the penetration: d<0 penetration, d>0 is no penetration
    # Forces: normal and tan
    if d>=0 or c>R:
        f_n=0
        f_t=0
    elif d<0:
        f_n=-K_finger*d-c_finger*v_d
        f_t=f_n*mu_finger*np.sign(-v_c)
    fy=-f_n*np.cos(theta)+f_t*np.sin(theta)
    fx=f_t*np.cos(theta)+f_n*np.sin(theta)
    if psi>=0 or d>=0:
        f_theta_left=0
    else:
        f_theta_left=K_torsion*psi
    print "left finger contact: f_n:",f_n,"f_t:",f_t
    f_theta=-f_t*p+f_n*c+f_theta_left
    print "moment of finger",f_theta
    return np.array([fx,fy,f_theta])

    
def right_finger_force(X,U):
    x,y,theta,x_dot,y_dot,theta_dot=X[0:6]
    c,d,D,psi=X[6:10]
    v_c,v_d,v_D,v_psi=U[0:4]
    # What is the relative velocity of right finger?
    # It depends on theta dot and also on v_c,v_d, and also v_D!
#    v_penetration=x_dot*np.sin(psi+theta)-y_dot*np.cos(theta+psi)-v_d*np.cos(psi)+v_c*np.sin(psi)-v_D
#    v_slide=x_dot*np.cos(psi+theta)+y_dot*np.sin(theta+psi)+v_d*np.sin(psi)+v_c*np.cos(psi)
    v_penetration=-v_d*np.cos(psi)-v_c*np.sin(psi)+v_D # Positive: going out
    v_slide=v_d*np.sin(psi)-v_c*np.cos(psi)
    # Now let's find the amount of penetration
    d=D-R-c*np.sin(psi)-d*np.cos(psi)
    print "right finger penertation:",d
    # Now let's find the contact forces
    if d>=0:
        f_n=0
        f_t=0
    elif d<0:
        f_n=-K_finger*d-c_finger*v_penetration
        f_t=f_n*mu_finger*np.sign(v_slide)
    print "right finger contact: f_n:",f_n,"f_t:",f_t
    # Now transfer the coordinates:
    fy=f_n*np.cos(theta+psi)+f_t*np.sin(theta+psi)
    fx=-f_n*np.sin(theta+psi)+f_t*np.cos(theta+psi)
    f_theta=f_n*p*np.sin(psi)+f_t*(R-p)*np.cos(psi)
    return np.array([fx,fy,f_theta])    


"""
Visualization Tools
"""    
def visualize(X,U,gripper=1):
    fig,ax1 = plt.subplots()
    ax1.set_xlabel("x",fontsize=20)
    ax1.set_ylabel("y",fontsize=20)
    ax1.set_xlim([-12,12])
    ax1.set_ylim([-2,10])
    fig.gca().set_aspect('equal')
    p_list=[]
    v=vertices(X)
    p_list.append(patches.Polygon(v, True))
    p=PatchCollection(p_list,color=(1,0,0),alpha=0.31,edgecolor=(1,0,0))
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax1.set_title("carrot %d"%t)
    ax1.plot([-12,12],[0,0],'black')
    center_mass(ax1,X)
    if gripper==1:
        # left Finger position
        finger_left(ax1,X)
        # Right finger
        finger_right(ax1,X)
    global t
    t=t+1
    fig.savefig('carrot_fig/carrot_%d.png'%t, dpi=100)
    plt.close()
    return fig

def vertices(X,N=50):
    x=X[0]
    y=X[1]
    theta=X[2]
    v=np.empty((50,2))
    for k in range(50):
        phi=-np.pi/2+np.pi/(N-1)*k
        v[k,0]=x+R*np.sin(phi+theta)
        v[k,1]=y-R*np.cos(phi+theta)
    return v

def center_mass(ax,X):
    x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi=X[0:10]
    ax.plot([x+p*np.sin(theta)],[y-p*np.cos(theta)],'+',color=(1,0,0))
        
def finger_left(ax,X,L=5.0,M=1.5):
    x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi=X[0:10]
    x_left=x-c*np.cos(theta)-d*np.sin(theta)
    y_left=y-c*np.sin(theta)+d*np.cos(theta)
    ax.plot([x_left,x_left+L*np.cos(theta+psi)],[y_left,y_left+L*np.sin(theta+psi)],'black')
#    ax.plot(x_left,y_left,'o')
    x_right=x_left-M*np.sin(psi+theta)
    y_right=y_left+M*np.cos(psi+theta)
    v1=[x_left,x_left+L*np.cos(theta+psi),x_right+L*np.cos(theta+psi),x_right]
    v2=[y_left,y_left+L*np.sin(theta+psi),y_right+L*np.sin(theta+psi),y_right]
    left_finger=[]
    left_finger.append(patches.Polygon(np.array([v1+v2]).reshape(2,4).T, True))
    ax.add_collection(PatchCollection(left_finger,color=(0,0,0),alpha=0.8,edgecolor=(0,0,0)))

def finger_right(ax,X,L=5.0,M=1.5):
    x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi=X[0:10]
    x_left=x-c*np.cos(theta)-d*np.sin(theta)
    y_left=y-c*np.sin(theta)+d*np.cos(theta)
    x_right=x_left+D*np.sin(psi+theta)
    y_right=y_left-D*np.cos(psi+theta)
    ax.plot([x_left,x_right],[y_left,y_right],'--',linewidth=0.1)
    ax.plot([x_right,x_right+L*np.cos(theta+psi)],[y_right,y_right+L*np.sin(theta+psi)],'black')
    x_left=x_right+M*np.sin(psi+theta)
    y_left=y_right-M*np.cos(psi+theta)
    v1=[x_left,x_left+L*np.cos(theta+psi),x_right+L*np.cos(theta+psi),x_right]
    v2=[y_left,y_left+L*np.sin(theta+psi),y_right+L*np.sin(theta+psi),y_right]
    left_finger=[]
    left_finger.append(patches.Polygon(np.array([v1+v2]).reshape(2,4).T, True))
    ax.add_collection(PatchCollection(left_finger,color=(0,0,0),alpha=0.8,edgecolor=(0,0,0)))
    