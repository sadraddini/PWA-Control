# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:06:20 2018

@author: sadra
"""

# External imports
import numpy as np

# My modules
from pypolycontain.lib.zonotope import zonotope
from pypolycontain.lib.polytope import polytope


# Internal imports
from PWA_lib.trajectory.system import system,linear_cell
from PWA_lib.trajectory.poly_trajectory import point_trajectory,polytopic_trajectory_given_modes

sys=system()
sys.name="inverted pendulum with two walls"


sys.A[0,0]=np.array([[1,0.01],[0.1,1]])
sys.B[0,0]=np.array([[0,0.01]]).T
sys.c[0,0]=np.array([[0,0]]).T

sys.A[1,0]=np.zeros((2,2))
sys.B[1,0]=np.zeros((2,1))
sys.c[1,0]=np.zeros((2,1))
sys.A[1,1]=np.array([[0,0],[-10,0]])
sys.B[1,1]=np.array([[0,0]]).T
sys.c[1,1]=np.array([[0,1]]).T


sys.A[2,0]=np.zeros((2,2))
sys.B[2,0]=np.zeros((2,1))
sys.c[2,0]=np.zeros((2,1))
sys.A[2,1]=np.array([[0,0],[-10,0]])
sys.B[2,1]=np.array([[0,0]]).T
sys.c[2,1]=np.array([[0,-1]]).T

H=np.array([[1,0,0],[0,1,0],[-1,0,0],[0,-1,0],[0,0,1],[0,0,-1]])
h=np.array([[0.12,1,0.12,1,4,4]]).T   
sys.C[0,0]=polytope(H,h)

H=np.array([[1,0,0],[0,1,0],[-1,0,0],[0,-1,0],[0,0,1],[0,0,-1]])
h=np.array([[0.1,1,0.12,1,4,4]]).T 
sys.C[1,0]=polytope(H,h)

H=np.array([[1,0,0],[0,1,0],[-1,0,0],[0,-1,0],[0,0,1],[0,0,-1]])
h=np.array([[0.12,1,0.1,1,4,4]]).T 
sys.C[2,0]=polytope(H,h)

H=np.array([[1,0,0],[0,1,0],[-1,0,0],[0,-1,0],[0,0,1],[0,0,-1]])
h=np.array([[0.12,1,-0.1,1,4,4]]).T  
sys.C[1,1]=polytope(H,h)

H=np.array([[1,0,0],[0,1,0],[-1,0,0],[0,-1,0],[0,0,1],[0,0,-1]])
h=np.array([[-0.1,1,0.12,1,4,4]]).T  
sys.C[2,1]=polytope(H,h)

sys.goal=zonotope(np.array([0,0.0]).reshape(2,1),np.array([[0,0],[0,0]]))

sys.n=2
sys.m=1
sys.list_of_sum_indices=[0,1,2]
sys.list_of_modes={}
sys.list_of_modes[0]=[0]
sys.list_of_modes[1]=[0,1]
sys.list_of_modes[2]=[0,1]

sys.build()


import matplotlib.pyplot as plt

x0=np.array([-0.00,-0.95]).reshape(2,1)
T=50
(x_n,u,delta_PWA,mu)=point_trajectory(sys,x0,[sys.goal],T)

plt.plot([x_n[t][0,0] for t in range(T+1)],[x_n[t][1,0] for t in range(T+1)])
plt.plot([x_n[t][0,0] for t in range(T+1)],[x_n[t][1,0] for t in range(T+1)],'+')
plt.plot([0.1,0.1],[-1,1],'black')
plt.plot([-0.1,-0.1],[-1,1],'black')
plt.plot([0],[0],'o')
#plt.plot([-0.05],[0.5],'o')


# Now build a funnel
from pypolycontain.utils.redundancy_reduction import canonical_polytope
list_of_cells=[]
for t in range(T):
    A=sum([sys.A[n,i]*delta_PWA[t,n,i] for n in sys.list_of_sum_indices for i in sys.list_of_modes[n]])
    B=sum([sys.B[n,i]*delta_PWA[t,n,i] for n in sys.list_of_sum_indices for i in sys.list_of_modes[n]])
    c=sum([sys.c[n,i]*delta_PWA[t,n,i] for n in sys.list_of_sum_indices for i in sys.list_of_modes[n]])
    H=np.vstack([sys.C[n,i].H for n in sys.list_of_sum_indices for i in sys.list_of_modes[n] if delta_PWA[t,n,i]==1])
    h=np.vstack([sys.C[n,i].h for n in sys.list_of_sum_indices for i in sys.list_of_modes[n] if delta_PWA[t,n,i]==1])
    H,h=canonical_polytope(H,h)
    cell=linear_cell(A,B,c,polytope(H,h))    
    list_of_cells.append(cell)

#for t in range(T):
#    cell=list_of_cells[t]
#    A,B,c,p=cell.A,cell.B,cell.c,cell.p
#    print all(abs(x_n[t+1]-np.dot(A,x_n[t])-np.dot(B,u[t])-c)<10**-8)
#    print all(p.h-np.dot(p.H,np.vstack((x_n[t],u[t])))>=-10**-8)
    

scale=np.array([0.12,1])    
(x,u,G,theta)=polytopic_trajectory_given_modes(x0,list_of_cells,sys.goal,eps=0.5,order=1,scale=scale)


from visualization.visualize import add_tube
fig,ax=plt.subplots()
ax.set_xlim([-0.12,0.12])
ax.set_ylim([-1,1])
add_tube(ax,x,G,eps=0.0001,list_of_dimensions=[0,1])
ax.plot([x[t][0,0] for t in range(T+1)],[x[t][1,0] for t in range(T+1)])
#ax.plot([x[t][0,0] for t in range(T+1)],[x[t][1,0] for t in range(T+1)],'+')
ax.plot([0.1,0.1],[-1,1],'black')
ax.plot([-0.1,-0.1],[-1,1],'black')
ax.plot([0],[0],'o')
plt.plot([x_n[t][0,0] for t in range(T+1)],[x_n[t][1,0] for t in range(T+1)])
plt.plot([x_n[t][0,0] for t in range(T+1)],[x_n[t][1,0] for t in range(T+1)],'+')
#plt.plot([-0.05],[0.5],'o')