#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 19:09:04 2019

@author: sadra
"""

import numpy as np
import pydrake.symbolic as sym

from PWA_lib.Manipulation.contact_point_pydrake import contact_point_symbolic
from PWA_lib.Manipulation.system_symbolic_pydrake import system_symbolic


mysystem=system_symbolic("my system")
q=np.array([sym.Variable("q%i"%i) for i in range(3)])

J=np.array([   q[0]**2*sym.sin(q[2])   , q[1]**3*sym.cos(q[2])  ])

z=sym.Jacobian(J,q)

E={q[i]:i+1 for i in range(3)}

z_v=sym.Evaluate(z,E)