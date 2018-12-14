# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 17:53:53 2018

@author: sadra
"""

from carrot_symbolic_flip import carrot_symbolic
from sympy import Symbol,pi,sin,cos,Function,diff,Matrix
import numpy as np

import time as time
from pypolycontain.lib.polytope import polytope
 
O=carrot_symbolic("my carrot")


class PWA_cell:
    def __init__(self,A,B,c,p,name="PWA cell"):
        self.A=A
        self.B=B
        self.c=c
        self.p=p
        self.name=name

def GRADIENT(phi,x):
#    print len(phi),len(x)
    D=np.empty((len(phi),(len(x))),dtype='object')
    for i in range(len(phi)):
        for j in range(len(x)):
            D[i,j]=diff(phi[i],x[j])
    return D
    
def SUB(symbolic_carrot,phi,x,u):
    # phi is a vector function of O.x and O.u
    q=np.zeros(phi.shape)
    mu_ground=0.3
    mu_finger=0.5
    
    K_ground=800/2
    c_ground=150/2
    K_finger=700
    c_finger=10*0
    
    K_torsion=5000
    c_torsion=100
    R=3
    for i in range(phi.shape[0]):
        for j in range(phi.shape[1]):
#            print phi[i,j]
            O=symbolic_carrot
            z=phi[i,j].subs({diff(O.c,O.t):O.c_dot,diff(O.d,O.t):O.d_dot})
#            print "1:",z
            z=z.subs({symbolic_carrot.U[i]:u[i] for i in range(len(u))})
#            print "3:",z
            z=z.subs({symbolic_carrot.X[i]:x[i] for i in range(len(x))})
#            print "3:",z
            z=z.subs({O.R:R,  O.K_ground:K_ground,   O.K_finger:K_finger,  O.c_finger:c_finger,\
                    O.c_ground:c_ground,   O.mu_ground:mu_ground,   O.mu_finger:mu_finger})
#            print "4:",z
            try:
                q[i,j]=float(z)
            except:
                q[i,j]=float(z[0])
    return q

def ground_polyhedrons(symbolic_carrot,X,U):
    """
    Returns polyhedron for ground force being zero
    """
    t=time.time()
    x,y,theta,x_dot,y_dot,theta_dot=X[0:6]
    c,d,D,psi=X[6:10]
    v_c,v_d,v_D,v_psi=U[0:4]
    assert theta<np.pi/2 # Curve side is below!
    phi,v_slide=symbolic_carrot.lambda_ground_curve()
    # Numbers first
    phi_num=SUB(symbolic_carrot,np.array([phi]).reshape(1,1),X,U)
    v_slide_num=SUB(symbolic_carrot,np.array([v_slide]).reshape(1,1),X,U)
    """
    Returns polyhedron for ground force being zero. No contact
    """
    H_x_symbolic=GRADIENT([phi],symbolic_carrot.X)
    H_u_symbolic=GRADIENT([phi],symbolic_carrot.U)
    H_x=SUB(symbolic_carrot,H_x_symbolic,X,U)
    H_u=SUB(symbolic_carrot,H_u_symbolic,X,U)
#    print H_x,H_x.shape
#    print H_u,H_u.shape
    H=np.hstack((H_x,H_u))
    h=np.dot(H_x,X)+np.dot(H_u,U)-phi_num
    poly_no_contact=polytope(-H,-h)
    """
    Returns polyhedron for ground force with sliding toward positive
    v_slide=positive
    H(x,u)<=h
    """
    G_x_symbolic=GRADIENT([v_slide],symbolic_carrot.X)
    G_u_symbolic=GRADIENT([v_slide],symbolic_carrot.U)
    G_x=SUB(symbolic_carrot,G_x_symbolic,X,U)
    G_u=SUB(symbolic_carrot,G_u_symbolic,X,U)
    G=np.hstack((G_x,G_u))
    g=np.dot(G_x,X)+np.dot(G_u,U)-v_slide_num
    poly_positive_slide=polytope(np.vstack((H,-G)),np.vstack((h,-g)))           
    poly_negative_slide=polytope(np.vstack((H,G)),np.vstack((h,g)))     
    print "time taken",time.time()-t, "seconds"
    return (poly_no_contact,poly_positive_slide,poly_negative_slide)
    
def left_polyhedrons(symbolic_carrot,X,U):
    """
    Returns polyhedron for ground force being zero
    """
    t=time.time()
    x,y,theta,x_dot,y_dot,theta_dot=X[0:6]
    c,d,D,psi=X[6:10]
    v_c,v_d,v_D,v_psi=U[0:4]
    assert theta<np.pi/2 # Curve side is below!
    phi,v_slide=symbolic_carrot.lambda_left_finger_point()
    # Numbers first
    phi_num=SUB(symbolic_carrot,np.array([phi]).reshape(1,1),X,U)
    v_slide_num=SUB(symbolic_carrot,np.array([v_slide]).reshape(1,1),X,U)
    """
    Returns polyhedron for ground force being zero. No contact
    """
    H_x_symbolic=GRADIENT([phi],symbolic_carrot.X)
    H_u_symbolic=GRADIENT([phi],symbolic_carrot.U)
    H_x=SUB(symbolic_carrot,H_x_symbolic,X,U)
    H_u=SUB(symbolic_carrot,H_u_symbolic,X,U)
#    print H_x,H_x.shape
#    print H_u,H_u.shape
    H=np.hstack((H_x,H_u))
    h=np.dot(H_x,X)+np.dot(H_u,U)-phi_num
    poly_no_contact=polytope(-H,-h)
    """
    Returns polyhedron for ground force with sliding toward positive
    v_slide=positive
    H(x,u)<=h
    """
    G_x_symbolic=GRADIENT([v_slide],symbolic_carrot.X)
    G_u_symbolic=GRADIENT([v_slide],symbolic_carrot.U)
    G_x=SUB(symbolic_carrot,G_x_symbolic,X,U)
    G_u=SUB(symbolic_carrot,G_u_symbolic,X,U)
    G=np.hstack((G_x,G_u))
    g=np.dot(G_x,X)+np.dot(G_u,U)-v_slide_num
    poly_positive_slide=polytope(np.vstack((H,-G)),np.vstack((h,-g)))           
    poly_negative_slide=polytope(np.vstack((H,G)),np.vstack((h,g)))     
    print "time taken",time.time()-t, "seconds"
    return (poly_no_contact,poly_positive_slide,poly_negative_slide)

def right_polyhedrons(symbolic_carrot,X,U):
    """
    Returns polyhedron for ground force being zero
    """
    t=time.time()
    x,y,theta,x_dot,y_dot,theta_dot=X[0:6]
    c,d,D,psi=X[6:10]
    v_c,v_d,v_D,v_psi=U[0:4]
    assert theta<np.pi/2 # Curve side is below!
    phi,v_slide=symbolic_carrot.lambda_right_finger()
    # Numbers first
    phi_num=SUB(symbolic_carrot,np.array([phi]).reshape(1,1),X,U)
    v_slide_num=SUB(symbolic_carrot,np.array([v_slide]).reshape(1,1),X,U)
    """
    Returns polyhedron for ground force being zero. No contact
    """
    H_x_symbolic=GRADIENT([phi],symbolic_carrot.X)
    H_u_symbolic=GRADIENT([phi],symbolic_carrot.U)
    H_x=SUB(symbolic_carrot,H_x_symbolic,X,U)
    H_u=SUB(symbolic_carrot,H_u_symbolic,X,U)
#    print H_x,H_x.shape
#    print H_u,H_u.shape
    H=np.hstack((H_x,H_u))
    h=np.dot(H_x,X)+np.dot(H_u,U)-phi_num
    poly_no_contact=polytope(-H,-h)
    """
    Returns polyhedron for ground force with sliding toward positive
    v_slide=positive
    H(x,u)<=h
    """
    G_x_symbolic=GRADIENT([v_slide],symbolic_carrot.X)
    G_u_symbolic=GRADIENT([v_slide],symbolic_carrot.U)
    G_x=SUB(symbolic_carrot,G_x_symbolic,X,U)
    G_u=SUB(symbolic_carrot,G_u_symbolic,X,U)
    G=np.hstack((G_x,G_u))
    g=np.dot(G_x,X)+np.dot(G_u,U)-v_slide_num
    poly_positive_slide=polytope(np.vstack((H,-G)),np.vstack((h,-g)))           
    poly_negative_slide=polytope(np.vstack((H,G)),np.vstack((h,g)))     
    print "time taken",time.time()-t, "seconds"
    return (poly_no_contact,poly_positive_slide,poly_negative_slide)
    
def force_ground(symbolic_carrot):
    J=symbolic_carrot.Jacobian_ground()
    phi,v_slide=symbolic_carrot.lambda_ground_curve()
    lambda_phi=-symbolic_carrot.K_ground*phi-symbolic_carrot.c_ground*diff(phi,symbolic_carrot.t) # Positive when penetrations is negative
    lambda_slide=-lambda_phi*symbolic_carrot.mu_ground
    return linearization(symbolic_carrot,J,lambda_phi,lambda_slide)


def force_left(symbolic_carrot):
    J=symbolic_carrot.Jacobian_left_finger()
    phi,v_slide=symbolic_carrot.lambda_left_finger_point()
    lambda_phi=-symbolic_carrot.K_finger*phi-symbolic_carrot.c_finger*diff(phi,symbolic_carrot.t) # Positive when penetrations is negative
    lambda_slide=-lambda_phi*symbolic_carrot.mu_finger
    return linearization(symbolic_carrot,J,lambda_phi,lambda_slide)


def force_right(symbolic_carrot):
    J=symbolic_carrot.Jacobian_right_finger()
    phi,v_slide=symbolic_carrot.lambda_right_finger()
    lambda_phi=-symbolic_carrot.K_finger*phi-symbolic_carrot.c_finger*diff(phi,symbolic_carrot.t) # Positive when penetrations is negative
    lambda_slide=-lambda_phi*symbolic_carrot.mu_finger
    # Three Forces
    return linearization(symbolic_carrot,J,lambda_phi,lambda_slide)

def linearization(symbolic_carrot,J,lambda_phi,lambda_slide):
    f0=np.dot(J,Matrix([np.zeros((2,1))]).reshape(2,1))
    f1=np.dot(J,Matrix([[lambda_phi,lambda_slide]]).reshape(2,1)) 
    f2=np.dot(J,Matrix([[lambda_phi,-lambda_slide]]).reshape(2,1))
    A0=GRADIENT(f0,symbolic_carrot.X)
    B0=GRADIENT(f0,symbolic_carrot.U)
    A1=GRADIENT(f1,symbolic_carrot.X)
    B1=GRADIENT(f1,symbolic_carrot.U)
    A2=GRADIENT(f2,symbolic_carrot.X)
    B2=GRADIENT(f2,symbolic_carrot.U)
    return [(A0,B0,f0),(A1,B1,f1),(A2,B2,f2)]

    
D0,D1,D2={},{},{}   
O=carrot_symbolic("my carrot")
for i in ["ground","left","right"]:
    if i=="ground":
        (D0[i],D1[i],D2[i])=force_ground(O)
    elif i=="right":
        (D0[i],D1[i],D2[i])=force_right(O)
    elif i=="left":
        (D0[i],D1[i],D2[i])=force_left(O)
 
def PWA(x,u):
    t=time.time()
    CELL={}
    (p0,p1,p2)=ground_polyhedrons(O,x,u)
    (q0,q1,q2)=left_polyhedrons(O,x,u)
    (r0,r1,r2)=right_polyhedrons(O,x,u)
    A0_num,B0_num,A1_num,B1_num,A2_num,B2_num,c0_num,c1_num,c2_num={},{},{},{},{},{},{},{},{}
    for i in ["ground","left","right"]: 
        A0_num[i],B0_num[i],c0_num[i]=[SUB(O,D0[i][j],x,u) for j in [0,1,2]]
        A1_num[i],B1_num[i],c1_num[i]=[SUB(O,D1[i][j],x,u) for j in [0,1,2]]
        A2_num[i],B2_num[i],c2_num[i]=[SUB(O,D2[i][j],x,u) for j in [0,1,2]]
    for i in ["ground","left","right"]: 
        if i=="ground":
            CELL[0,i]=PWA_cell(A0_num[i],B0_num[i],c0_num[i],p0,name=i+" no contact")
            CELL[1,i]=PWA_cell(A1_num[i],B1_num[i],c1_num[i],p1,name=i+" positive slide")
            CELL[2,i]=PWA_cell(A2_num[i],B2_num[i],c2_num[i],p2,name=i+" negative slide")
        elif i=="left":
            CELL[0,i]=PWA_cell(A0_num[i],B0_num[i],c0_num[i],q0,name=i+" no contact")
            CELL[1,i]=PWA_cell(A1_num[i],B1_num[i],c1_num[i],q1,name=i+" positive slide")
            CELL[2,i]=PWA_cell(A2_num[i],B2_num[i],c2_num[i],q2,name=i+" negative slide")
        elif i=="right":
            CELL[0,i]=PWA_cell(A0_num[i],B0_num[i],c0_num[i],r0,name=i+" no contact")
            CELL[1,i]=PWA_cell(A1_num[i],B1_num[i],c1_num[i],r1,name=i+" positive slide")
            CELL[2,i]=PWA_cell(A2_num[i],B2_num[i],c2_num[i],r2,name=i+" negative slide")
    print "time taken",time.time()-t, "seconds"
    return CELL