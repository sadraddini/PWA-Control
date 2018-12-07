# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 11:33:29 2018

@author: sadra
"""

from sympy import Symbol,pi,sin,cos
import numpy as np

x=Symbol('x')
y=Symbol('y')
theta=Symbol('theta')
x_dot=Symbol('x_dot')
y_dot=Symbol('y_dot')
theta_dot=Symbol('theta_dot')
c=Symbol('c')
d=Symbol('d')
D=Symbol('D')
psi=Symbol('psi')

v_c=Symbol('v_c')
v_d=Symbol('v_d')
v_D=Symbol('v_D')
v_psi=Symbol('v_psi')

R=Symbol('R')
K_ground=Symbol('K_ground')
c_ground=Symbol('c_ground')
mu_ground=Symbol('mu_ground')
K_finger=Symbol('K_finger')
c_finger=Symbol('c_finger')
mu_finger=Symbol('mu_finger')
p=4*R/(3*pi)

# Mode 1: Carrot cylinder side down, no finger
X=np.array([x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi])
U=np.array([v_c,v_d,v_D,v_psi])

def carrot_ground_down_slide_left(X):
    x,y,theta,x_dot,y_dot,theta_dot=X[0:6]
    """
    v_contact<0
    """
    d=y-R
    v_d=y_dot
    v_contact=x_dot+R*theta_dot
    fy=-K_ground*d-c_ground*v_d
    fx=fy*mu_ground # v_contact<0
    f_theta=-fy*p*sin(theta)+fx*(R-p*cos(theta))
    return np.array([fx,fy,f_theta])

def carrot_ground_down_slide_right(X):
    x,y,theta,x_dot,y_dot,theta_dot=X[0:6]
    """
    v_contact<0
    """
    d=y-R
    v_d=y_dot
    v_contact=x_dot+R*theta_dot
    fy=-K_ground*d-c_ground*v_d
    fx=fy*mu_ground # v_contact<0
    f_theta=-fy*p*sin(theta)+fx*(R-p*cos(theta))
    return np.array([fx,fy,f_theta])
   
def left_finger_force_slide_positive(X,U):
    x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi=X[0:10]
    v_c,v_d,v_D,v_psi=U[0:4]
    f_n=-K_finger*d-c_finger*v_d
    f_t=f_n*mu_finger*(-1)
    fy=-f_n*cos(theta)+f_t*sin(theta)
    fx=f_t*cos(theta)+f_n*sin(theta)
    f_theta=-f_t*p+f_n*c
    return np.array([fx,fy,f_theta]) 
    
def left_finger_force_slide_right(X,U):
    x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi=X[0:10]
    v_c,v_d,v_D,v_psi=U[0:4]
    f_n=-K_finger*d-c_finger*v_d
    f_t=f_n*mu_finger*(1)
    fy=-f_n*cos(theta)+f_t*sin(theta)
    fx=f_t*cos(theta)+f_n*sin(theta)
    f_theta=-f_t*p+f_n*c
    return np.array([fx,fy,f_theta]) 

def right_finger_force_right(X,U):
    # v_slide>0
    x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi=X[0:10]
    v_c,v_d,v_D,v_psi=U[0:4]
    v_penetration=-v_d*cos(psi)-v_c*sin(psi)+v_D # Positive: going out
    v_slide=sin(psi)-cos(psi)
    # Now let's find the amount of penetration
    d=D-R-c*sin(psi)-d*cos(psi)
    # Now let's find the contact forces
    f_n=-K_finger*d-c_finger*v_penetration
    f_t=f_n*mu_finger*v_slide
    fy=f_n*cos(theta+psi)+f_t*sin(theta+psi)
    fx=-f_n*sin(theta+psi)+f_t*cos(theta+psi)
    f_theta=f_n*p*sin(psi)+f_t*(R-p)*cos(psi)
    return np.array([fx,fy,f_theta])
    
def right_finger_force_left(X,U):
    # v_slide<=0
    x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi=X[0:10]
    v_c,v_d,v_D,v_psi=U[0:4]
    v_penetration=-v_d*cos(psi)-v_c*sin(psi)+v_D # Positive: going out
    v_slide=sin(psi)-cos(psi)
    # Now let's find the amount of penetration
    d=D-R-c*sin(psi)-d*cos(psi)
    # Now let's find the contact forces
    f_n=-K_finger*d-c_finger*v_penetration
    f_t=-f_n*mu_finger*v_slide
    fy=f_n*cos(theta+psi)+f_t*sin(theta+psi)
    fx=-f_n*sin(theta+psi)+f_t*cos(theta+psi)
    f_theta=f_n*p*sin(psi)+f_t*(R-p)*cos(psi)
    return np.array([fx,fy,f_theta])

def determine_mode(X,U):
    """
    The mode is (g,f_l,f_r)
    1: slide to right
    0: no contact
    -1: slide to left
    """
    x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi=X[0:10]
    v_c,v_d,v_D,v_psi=U[0:4]
    
    