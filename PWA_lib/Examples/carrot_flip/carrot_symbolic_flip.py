# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 14:12:54 2018

@author: sadra
"""
from sympy import Symbol,pi,sin,cos,Function,diff,Matrix
import numpy as np


#X=np.array([x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi])
#U=np.array([v_c,v_d,v_D,v_psi])
#
#
#x=Symbol('x')
#y=Symbol('y')
#theta=Symbol('theta')
#x_dot=Symbol('x_dot')
#y_dot=Symbol('y_dot')
#theta_dot=Symbol('theta_dot')
#c=Symbol('c')
#d=Symbol('d')
#D=Symbol('D')
#psi=Symbol('psi')
#
#v_c=Symbol('v_c')
#v_d=Symbol('v_d')
#v_D=Symbol('v_D')
#v_psi=Symbol('v_psi')
#
#R=Symbol('R')
#K_ground=Symbol('K_ground')
#c_ground=Symbol('c_ground')
#mu_ground=Symbol('mu_ground')
#K_finger=Symbol('K_finger')
#c_finger=Symbol('c_finger')
#mu_finger=Symbol('mu_finger')
#p=4*R/(3*pi)
#
## Mode 1: Carrot cylinder side down, no finger
#X=np.array([x,y,theta,x_dot,y_dot,theta_dot,c,d,D,psi])
#U=np.array([v_c,v_d,v_D,v_psi])

class carrot_symbolic:
    def __init__(self,name):
        self.name=name
#        self.x=Symbol('x')
#        self.y=Symbol('y')
#        self.theta=Symbol('theta')
#        self.x_dot=Symbol('x_dot')
#        self.y_dot=Symbol('y_dot')
#        self.theta_dot=Symbol('theta_dot')
#        self.c=Symbol('c')
#        self.d=Symbol('d')
#        self.D=Symbol('D')
#        self.psi=Symbol('psi')
#        
#        self.v_c=Symbol('v_c')
#        self.v_d=Symbol('v_d')
#        self.v_D=Symbol('v_D')
#        self.v_psi=Symbol('v_psi')
#        
#        self.X=np.array([self.x,self.y,self.theta,self.x_dot,self.y_dot,self.theta_dot,self.c,self.d,self.D,self.psi])
#        self.U=np.array([self.v_c,self.v_d,self.v_D,self.v_psi])
        
        self.R=Symbol('R')
        self.g=Symbol('g')
        self.K_ground=Symbol('K_ground')
        self.c_ground=Symbol('c_ground')
        self.mu_ground=Symbol('mu_ground')
        self.K_finger=Symbol('K_finger')
        self.c_finger=Symbol('c_finger')
        self.mu_finger=Symbol('mu_finger')
        self.p=4*self.R/(3*pi)
        self.I=1/6.0*self.R**2
        
        self.t=Symbol("t")
        t=self.t
        self.x=Function("x")(t)
        self.y=Function("y")(t)
        self.theta=Function("theta")(t)
        self.c=Function("c")(t)
        self.d=Function("d")(t)
        self.D=Function("D")(t)
        self.psi=Function("psi")(t)
        
        self.x_dot,self.y_dot,self.theta_dot=diff(self.x,t),diff(self.y,t),diff(self.theta,t)
        self.c_dot,self.d_dot,self.D_dot,self.psi_dot=diff(self.c,t),diff(self.d,t),diff(self.D,t),diff(self.psi,t)
        self.X=[self.x,self.y,self.theta,self.x_dot,self.y_dot,self.theta_dot,self.c,self.d,self.D,self.psi]
        self.U=[self.c_dot,self.d_dot,self.D_dot,self.psi_dot]
        

    def f_g(self):
        self.f_g=np.zeros(10)
        self.f_g[4]=-self.g
        return self.f_g
        
    def lambda_ground_curve(self):
        phi=self.y-self.R
        t=self.t
        v_x=diff(self.x,t)
        v_theta=diff(self.theta,t)
        slide=v_x+self.R*v_theta
#        lambda_phi=-self.K_ground*phi # Positive when penetrations is negative
#        lambda_slide=lambda_phi*v_slide
        return (phi,slide)
        
    def lambda_left_finger_point(self):
        phi=self.d
        slide=diff(self.c,self.t)
        return (phi,slide)
        
        
    def lambda_right_finger(self):
        phi=self.D-self.R-self.c*sin(self.psi)-self.d*cos(self.psi)
        v_d=diff(self.d,self.t)
        v_c=diff(self.c,self.t)
        v_theta=diff(self.theta,self.t)
        slide=-self.R*v_theta-v_d*sin(self.psi)+v_c*cos(self.psi)
        #+v_d*sin(self.psi)-v_c*cos(self.psi)
        return (phi,slide)
        
    
    def Jacobian_ground(self):
        J=Matrix(np.zeros((10,2)))
        J[4,0]=1
        J[3,1]=1
        J[5,0]=-self.p*sin(self.theta)/self.I
        J[5,1]=self.R-self.p*cos(self.theta)/self.I 
        return J
        
    def Jacobian_left_finger(self):
        J=Matrix(np.zeros((10,2)))
        J[3,0],J[3,1]=-cos(self.theta),-sin(self.theta)
        J[4,0],J[4,1]=sin(self.theta),-cos(self.theta)
        J[5,0]=self.c/self.I
        J[5,1]=self.p/self.I
        return J
        
    def Jacobian_right_finger(self):
        J=Matrix(np.zeros((10,2)))
        J[3,0],J[3,1]=-cos(self.theta+self.psi),-sin(self.theta+self.psi)
        J[4,0],J[4,1]=-sin(self.theta+self.psi),cos(self.theta+self.psi)
        J[5,0]=-self.p*sin(self.psi)/self.I
        J[5,1]=-1/self.I*(self.R-self.p*cos(self.psi))
        return J
        
    
        
        
        
        
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
    
    