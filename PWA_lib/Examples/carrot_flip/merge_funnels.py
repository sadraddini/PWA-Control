# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 09:28:31 2018

@author: sadra
"""
import pickle
from PWA_lib.visualization.visualize import add_tube
import matplotlib.pyplot as plt

x_funnel_1=pickle.load(open("x_funnel_1.pkl","r"))
x_funnel_2=pickle.load(open("x_funnel_2.pkl","r"))
x_funnel_3=pickle.load(open("x_funnel_3.pkl","r"))
x_funnel_4=pickle.load(open("x_funnel_4.pkl","r"))

G_funnel_1=pickle.load(open("G_funnel_1.pkl","r"))
G_funnel_2=pickle.load(open("G_funnel_2.pkl","r"))
G_funnel_3=pickle.load(open("G_funnel_3.pkl","r"))
G_funnel_4=pickle.load(open("G_funnel_4.pkl","r"))

fig,ax1 = plt.subplots()
ax1.set_xlim([-0.025,0.1])
ax1.set_ylim([2.925,3])
#fig.gca().set_aspect('equal')
#    print "x_funnel is",x_funnel
#    print "G_funnel is",G_funnel
add_tube(ax1,x_funnel_1,G_funnel_1,eps=0.0002,list_of_dimensions=[0,1],axis=2)
add_tube(ax1,x_funnel_2,G_funnel_2,eps=0.0002,list_of_dimensions=[0,1],axis=2)
add_tube(ax1,x_funnel_3,G_funnel_3,eps=0.0001,list_of_dimensions=[0,1],axis=2)
add_tube(ax1,x_funnel_4,G_funnel_4,eps=0.0001,list_of_dimensions=[0,1],axis=2)
ax1.set_xlabel(r"x",fontsize=20)
ax1.set_ylabel(r"y",fontsize=20)
ax1.set_title(r"Central Polytopic Trajectory",fontsize=20)
