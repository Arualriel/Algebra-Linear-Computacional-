# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 15:37:10 2016

@author: gastao
"""

import numpy as np

from scipy import linalg


A = np.matrix([[1,1,0],[1,1,1],[0,1,1]])

Q,R = linalg.qr(A)

print (Q,R)

for i in range(10):
    Q,R = linalg.qr(A)
    
    A1 = np.dot(R,Q)
    
    print (A1)    

    print( linalg.norm(A - A1) )   
    
    A = A1