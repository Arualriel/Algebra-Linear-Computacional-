#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 20:37:09 2018

@author: laura
"""

import numpy as np
from scipy import linalg


def Maior(y,ordem):
    maior=0
    for i in range(ordem):
        if y[i]>=maior:
            maior=y[i]
    return maior

def Normalizar(y,ordem,sigma):
    l=np.zeros((ordem,1))
    for i in range(ordem):
        l[i]=(y[i]/sigma)
    return l

def Erro(n,A):
    epsilon=0.000001
    k=2
    A1=np.identity
    A2=A
    erro=10
    while erro>=epsilon:
        A1=A2
        A2=A2*A
        erro=np.linalg.norm(A2-A1,np.inf)
        k=k+1
    return k,erro,epsilon

def Potencias(A,ordem):
    z=np.ones((ordem,1))
    zk=np.zeros((ordem,1))
    sigmak=0
    norma=10
    s=0
    sigma=0
    modulo=10
    epsilon=0.000001
 
    while (norma>=epsilon)and(modulo>=epsilon):
        
        y=np.dot(A,z)
        zk=z
        
        sigmak=sigma
        s=Maior(y,ordem)
        sigma=s[0]
        z=Normalizar(y,ordem,sigma)
        norma=linalg.norm(z-zk)
        modulo=abs(sigma-sigmak)

    zk=z*(1/linalg.norm(z))
    
    return zk

n=20

#A=np.matrix([[0.8,0.3],[0.2,0.7]])
#A=np.matrix([[0.2,0.4,0.1,0.4],[0.3,0.2,0.2,0.1],[0.1,0.3,0.5,0.3],[0.4,0.1,0.2,0.2]])

A=np.matrix([[1/3,1/3,0.0,1/5,0.0,0.0,0.0,0.0],
             [1/3,1/3,0.0,0.0,1/4,0.0,0.0,0.0],
             [0.0,0.0,1/3,1/5,0.0,1/3,0.0,0.0],
             [1/3,0.0,1/3,1/5,1/4,0.0,0.0,1/3],
             [0.0,1/3,0.0,1/5,1/4,0.0,0.0,1/3],
             [0.0,0.0,1/3,0.0,0.0,1/3,1/4,0.0],
             [0.0,0.0,0.0,1/5,0.0,1/3,1/4,1/3],
             [0.0,0.0,0.0,0.0,1/4,0.0,1/4,1/3]])
    
num=int(np.sqrt(np.size(A)))    

autovalores, P = np.linalg.eig(A)
k,erro,epsilon=Erro(n,A)

D=np.zeros((num,num))
for i in range(num):
    D[i,i]=autovalores[i]

#x0=np.matrix([[0.25],[0.25],[0.25],[0.25]])
#x0=np.matrix([[0.5],[0.5]])
x0=np.matrix([[0.0],[0.25],[0.0],[0.25],[0.25],[0.0],[0.0],[0.25]])
l1=Potencias(A,num)
An=P*(D**n)*np.linalg.inv(P)
xn=P*(D**n)*np.linalg.inv(P)*x0
print("************Cadeia de Markov************")
print("\nA matriz convergiu em",k,"iteracoes, com erro=",erro,"para epsilon=",epsilon)
print("\nn=",n,"\n \nA^",n,"= \n",An)
print("\nAutovetor dominante l1 (autovalor=1) = \n",l1)
print("\nO",n,"-esimo estado da cadeia é o vetor x",n,"=\n",xn)