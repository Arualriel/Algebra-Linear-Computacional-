

import numpy as np

def Tmin(A,rk1):
    rr = np.dot(rk1.T,rk1)
    Ar = np.dot(A,rk1)
    print "Ar(1) = \n",Ar
    Arr = np.dot(Ar.T,rk1)
    
    return np.float(rr/Arr)

def erroVk(vk,vk1):
    numerador = np.linalg.norm(vk - vk1,np.inf)
    denominador = np.linalg.norm(vk,np.inf)

    return numerador/denominador

    



##########################################
##########################################
# DADOS INICIAIS


A = np.matrix([[10,1,0],[1,10,1],[0,1,10]])
#b = np.matrix([[1],[2],[3]])
b = np.matrix([[11],[11],[1]])



v0 = np.zeros((3,1))
epsilon = 0.1

kmax = 100

r0 = np.dot(A,v0) - b

for k in range(kmax):
    
    tmin = Tmin(A,r0)
   
    v1 = v0 - tmin*r0
        
    pk = -r0
    
    r1 = r0 + tmin*np.dot(A,pk)

    erroR = np.linalg.norm(r1)
    erroV = erroVk(v1,v0)        
    

    if (erroR < epsilon) or (erroV < epsilon):
        print "v[",k,"] = \n",v1
        break
    else:
        v0 = v1
        r0 = r1


print "Atingiu o numero maximo de iteracoes sem convergir!"        
