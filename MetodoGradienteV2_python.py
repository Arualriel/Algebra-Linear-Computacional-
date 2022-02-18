
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


A = np.matrix([[10,1,0],[1,10,1],[0,1,10]])
#b = np.matrix([[1],[2],[3]])
b = np.matrix([[11],[11],[1]])



v0 = np.ones((3,1))
epsilon = 0.0001

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



X,Y,Z = sol[0,:].tolist(),sol[1,:].tolist(),sol[2,:].tolist()

X,Y,Z = X[0],Y[0],Z[0]

        
        
        
        
        

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(X,Y,Z)
ax.plot(X,Y,Z,c="r")

for i in range(N):
    ax.text(X[i],Y[i],Z[i],r"$v_%d$"%(i))



plt.show()

  
