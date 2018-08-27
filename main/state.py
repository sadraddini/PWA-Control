#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 18:34:17 2018

@author: sadra
"""

### Internal imports
import sys
sys.path.append('../..')


class state:
    def __init__(self,x,G,mode,ID,t,character,epsilon=10**-8):
        self.name="s"+str(ID)+"-"+str(mode)+"-"+str(t)+"-"+str(character)
        self.x=x
        self.G=G
        self.Ginv=np.linalg.pinv(G)
        if np.linalg.det(G)<10**-6:
            (U,s,V)=np.linalg.svd(G)
            self.G_eps=np.dot(U,np.dot(np.diag(s+epsilon),V))
            self.G_eps_inv=np.linalg.pinv(self.G_eps)
        else: # Large volume:
            self.G_eps=self.G
            self.G_eps_inv=self.Ginv
        if G.shape[0]==G.shape[1]:
            self.volume=abs(np.linalg.det(G))
        self.volume_flag=self.volume>10**-8
        self.mode=mode
        self.ID=ID
        self.t=t
        self.character=character # 0 for goal, 1 for ordinary paralleltope, 2 for 
        v=vertices_cube(G.shape[1])
        self.vertices=(np.dot(G,v.T)).T
        self.backward_zero=-1 # 0 for not computed, -1 for computed but not useful, 1 some large regions leading it have been computed!
        self.cost_to_go=0
        self.time_to_go=0
        self.cost_to_child=0
        self.parent=[]
            
    def __repr__(self):
        return self.name