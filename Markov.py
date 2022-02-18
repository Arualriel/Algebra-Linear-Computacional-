#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 23 18:24:54 2018

@author: laura
"""

import numpy as np

from scipy import linalg

def Autovalores(A,ordem):
    epsilon=0.000000001
    Q,R = linalg.qr(A)
    l=np.zeros(ordem)
    norm=1
    i=0
    while (i<=100)and(norm>=epsilon):
        Q,R = linalg.qr(A)
        
        A1 = np.dot(R,Q)
    
        norm=(linalg.norm(A - A1))    
        
        A = A1
        i=i+1
        
    for i in range(ordem):
        l[i]=A[i,i]
    return l

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

def Potencias(A,ordem):
    z=np.ones((ordem,1))
    zk=np.zeros((ordem,1))
    sigmak=0
    norma=10
    s=0
    sigma=0
    modulo=10
    epsilon=0.00000001
    
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




def SVD(A,AtA,ordemAtA,lAtA,r,s):
    pare=0
    U=np.zeros((ordemAtA,ordemAtA))
    Vt=np.zeros((ordemAtA,ordemAtA))
    S=np.zeros((r,s))#*************ordem de sigma
    v=np.zeros((ordemAtA,1))
    w=np.zeros((ordemAtA,1))
    I=np.identity(ordemAtA,float)
    mi=np.zeros(ordemAtA)
    for i in range(ordemAtA):
        if lAtA[i]==0:
            pare=1
    if pare!=1:
        for i in range(ordemAtA):
            
            mi[i]=np.sqrt(lAtA[i])
            v=Gradiente(AtA-lAtA[i]*I,ordemAtA)
            print(v)
            w=(1/mi[i])*np.dot(A,v)
            print(w)
            for j in range(ordemAtA):
                U[j,i]=w[j]
                Vt[j,i]=v[j]
            S[i,i]=mi[i]
            
    else:
        U,S,Vt=linalg.svd(AtA)

    return U,S,Vt




def Gradiente(A,ordem):
    norma=0
    epsilon=0.00000000001
    tmin=0
    t=0
    pare=0
    r=np.zeros(ordem)
    errorel=0
    k=0
    Ar=np.zeros(ordem)
    
    
    b=np.zeros((ordem,1))
    v=np.ones((ordem,1))
    vk=np.zeros((ordem,1))
    r=np.dot(A,v)-b
    norma=linalg.norm(r)
    k=0
    num=1000000
    pare=0
    if(norma!=0):
        while (k<=num)and(pare==0):
            Ar=np.dot(A,r)
            tmin=(((linalg.norm(r))**2)/(np.dot(Ar.T,r)))
            t=tmin[0,0]
            vk=v-t*r
            r=r-t*Ar
            norma=linalg.norm(r,np.inf)
            errorel=((linalg.norm(vk-v))/(linalg.norm(vk)))
            if(norma<epsilon) and (errorel<epsilon):
                pare=1
            v=vk
            k=k+1
    else:
        print("a norma do residuo é zero!")
    norma=linalg.norm(vk)
    vk=vk*(1/norma)
#        for i in range(ordem):
#            
#            S[i,j]=(vk[i]/norma)
#    
        
    return vk

#def Markov(n,S,l,ordem): 
#    Pn=np.zeros((ordem,ordem))
#    
#    for i in range(ordem):
#        Pn[i,i]=l[i]**n
#    An=np.dot(S,np.dot(Pn,linalg.inv(S)))
#    return An
    
def Markovs(U,S,Vt,n,r,s):
    Si=np.zeros((r,s))
    for i in range(r):
        for j in range(s):
            if(i==j):
                Si[i,j]=S[i,i]**n
    An=np.dot(U,np.dot(Si,Vt))
    return An




n=2
A = np.matrix([[3/10,7/10,0],[7/20,12/20,1/20],[0,3/7,4/7]])
AtA=np.dot(A.T,A)
ordem=int(np.sqrt(np.size(A)))
r,s=np.shape(A)
ordemAtA=int(np.sqrt(np.size(AtA)))
l=Autovalores(A,ordem)
lAtA=Autovalores(AtA,ordemAtA)
U,S,Vt=SVD(A,AtA,ordemAtA,lAtA,r,s)

Vt=Vt.T
An=Markovs(U,S,Vt,n,r,s)
An=An.T
print("A^",n,"=",An)
print("U=",U,"\n S=",S,"\n Vt=",Vt)
x=Potencias(A,ordem)
print ("vetor dominante=",x)
print("\n\n U*Vt=\n",np.dot(Vt,U))
U,S,Vt=linalg.svd(A)
print("U=",U,"\n S=",S,"\n Vt=",Vt)
