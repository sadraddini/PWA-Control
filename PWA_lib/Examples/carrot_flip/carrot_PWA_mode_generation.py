# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 09:50:37 2018

@author: sadra
"""
from carrot_symbolic_flip import carrot_symbolic
from sympy import Symbol,pi,sin,cos,Function,diff,Matrix
import numpy as np

import time as time
from pypolycontain.lib.polytope import polytope
 
O=carrot_symbolic("my carrot")


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
            z=z.subs({symbolic_carrot.R:3,symbolic_carrot.t:0})
#            print "4:",z
            q[i,j]=float(z)
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

#def polyhedrons_of_contact(symbolic_carrot,phi,v_slide,X,U):
#    # Numbers first
#    phi_num=SUB(symbolic_carrot,np.array([phi]).reshape(1,1),X,U)
#    v_slide_num=SUB(symbolic_carrot,np.array([v_slide]).reshape(1,1),X,U)
#    """
#    Returns polyhedron for ground force being zero. No contact
#    """
#    H_x_symbolic=GRADIENT([phi],symbolic_carrot.X)
#    H_u_symbolic=GRADIENT([phi],symbolic_carrot.U)
#    H_x=SUB(symbolic_carrot,H_x_symbolic,X,U)
#    H_u=SUB(symbolic_carrot,H_u_symbolic,X,U)
##    print H_x,H_x.shape
##    print H_u,H_u.shape
#    H=np.hstack((H_x,H_u))
#    h=np.dot(H_x,X)+np.dot(H_u,U)-phi_num
#    print "time taken",time.time()-t, "seconds"
#    poly_no_contact=polytope(-H,-h)
#    """
#    Returns polyhedron for ground force with sliding toward positive
#    v_slide=positive
#    H(x,u)<=h
#    """
#    G_x_symbolic=GRADIENT([v_slide],symbolic_carrot.X)
#    G_u_symbolic=GRADIENT([v_slide],symbolic_carrot.U)
#    G_x=SUB(symbolic_carrot,G_x_symbolic,X,U)
#    G_u=SUB(symbolic_carrot,G_u_symbolic,X,U)
#    G=np.hstack((G_x,G_u))
#    g=np.dot(G_x,X)+np.dot(G_u,U)-v_slide_num
#    poly_positive_slide=polytope(np.vstack((H,-G)),np.vstack((h,-g)))           
#    poly_negative_slide=polytope(np.vstack((H,G)),np.vstack((h,g)))     
#    return (poly_no_contact,poly_positive_slide,poly_negative_slide)    
    
#x=X_traj[:,150].reshape(10,1)
#u=U_traj[:,150].reshape(4,1)
#(p0,p1,p2)=ground_polyhedrons(O,x,u)
#(q0,q1,q2)=left_polyhedrons(O,x,u)
#(r0,r1,r2)=right_polyhedrons(O,x,u)