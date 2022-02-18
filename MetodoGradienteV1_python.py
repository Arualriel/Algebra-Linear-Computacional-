
import numpy as np

def Tmin(A,rk1):
    rr = np.dot(rk1.T,rk1)
    Ar = np.dot(A,rk1)
    Arr = np.dot(Ar.T,rk1)
    
    return np.float(rr/Arr)

def erroVk(vk,vk1):
    numerador = np.linalg.norm(vk - vk1)
    denominador = np.linalg.norm(vk)

    return numerador/denominador

    



##########################################
##########################################
# DADOS INICIAIS


A = np.matrix([[10,1],[1,10]])
#b = np.matrix([[1],[2],[3]])
b = np.matrix([[11],[1]])



v0 = 20*np.ones((2,1))
#v0 = np.zeros((2,1))

epsilon = 0.000001

kmax = 100

r0 = np.dot(A,v0) - b

sol = v0


for k in range(kmax):
    
    tmin = Tmin(A,r0)
   
    v1 = v0 - tmin*r0
        
    pk = -r0
    
    r1 = r0 + tmin*np.dot(A,pk)

    erroR = np.linalg.norm(r1)
    erroV = erroVk(v1,v0)        
    

    if (erroR < epsilon) or (erroV < epsilon):
        print "v[",k,"] = \n",v1
        sol = np.concatenate((sol,v1),axis=1)
        N = k
        break
    else:
        v0 = v1        
        r0 = r1
        
        sol = np.concatenate((sol,v1),axis=1)

print "Atingiu o numero maximo de iteracoes sem convergir!"              



X,Y = sol[0,:].tolist(),sol[1,:].tolist()

X,Y = X[0],Y[0]

numX = len(X)
Z = numX*[0]
        
        
        
        
        

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(X,Y,Z)
ax.plot(X,Y,Z,c="r")

for i in range(N):
    ax.text(X[i],Y[i],Z[i],r"$v_%d$"%(i))



def F(x,y):
        
    X=np.matrix([[x],[y]])
    
    return 0.5*np.dot(X.T,np.dot(A,X)) - np.dot(X.T,b)


Xf = np.arange(-3,3,0.1)
Yf = Xf

X,Y = np.meshgrid(Xf,Yf)

n,m = np.shape(X)
Z = np.zeros((n,m))


for i in range(n):
    for j in range(m):
        Z[i,j] = F(X[i,j],Y[i,j])


ax.plot_surface(X, Y, Z, rstride=10, cstride=10, 
                       linewidth=0.1, antialiased=False,alpha=0.2)

ax.contour(X, Y, Z, zdir='z', offset=0)




plt.show()

  
