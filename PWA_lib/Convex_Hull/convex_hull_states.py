#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 15:28:10 2018

@author: sadra
"""

class hull_state:
    def __init__(self,list_of_pre,list_of_post):
        self.pres=list_of_pre
        self.posts=list_of_post
        self.pres_mode=list_of_pre[0].mode
        self.posts_mode=list_of_post[0].mode