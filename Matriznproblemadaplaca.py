#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 01:14:02 2018

@author: laura
"""

import numpy as np

n=3
N=n*n

M=np.zeros((N,N))

for i in range(N):
    for j in range(N):
        if i==j:
            M[i,j]=-4
        else:
            if i==(j-1):
                if i%2==0:
                    M[i,j]=1
                else:
                    M[i,j]=0
            else:
                if i==(j+1):
                    if i%2==1:
                        M[i,j]=1
                    else:
                        M[i,j]=0
                else:
                    if i==(j+2):
                        M[i,j]=1
                    else:
                        if i==(j-2):
                            M[i,j]=1


print(M[:,:])