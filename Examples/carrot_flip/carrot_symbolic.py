# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 11:33:29 2018

@author: sadra
"""

from sympy import *
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

R=Symbol('R')
K_ground=Symbol('K_ground')
c_ground=Symbol('c_ground')
mu_ground=Symbol('mu_ground')
p=4*R/(3*np.pi)

# Mode 1: Carrot cylinder side down, no finger
def carrot_ground_down_slide_right(x,y,theta,x_dot,y_dot,theta_dot):
    """
    v_contact>0
    """
    d=y-R
    v_d=y_dot
    v_contact=x_dot+R*theta_dot
    fy=-K_ground*d-c_ground*v_d
    fx=-fy*mu_ground # v_contact>0
    f_theta=-fy*p*sin(theta)+fx*(R-p*cos(theta))
    return np.array([fx,fy,f_theta])

def carrot_ground_down_slide_left(x,y,theta,x_dot,y_dot,theta_dot):
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

def carrot_ground_down_slide_left(x,y,theta,x_dot,y_dot,theta_dot):
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
    
def 