
 
#biblioteca numerica
import numpy as np
 
#biblioteca para graficos
import matplotlib.pyplot as plt

from matplotlib import cm
 
from mpl_toolkits.mplot3d import Axes3D
 
######################################
### Definindo funcoes 
######################################
 
def montaMatriz(n):
  N = n*n
  M = np.zeros((N,N))
  D = np.zeros((N,N))
  
  for i in range(N):
    D[i,i] = -4.0
    
    if i < n*(n-1):
      M[i,i+n] = 1.0
      M[i+n,i] = 1.0
      
    if (i + 1)%n != 0 :
      M[i,i+1] = 1.0
      M[i+1,i] = 1.0
  
  return M,D

def montaB(n,e0,d0,c0,b0):
    N = n*n
    b    = np.zeros((N,1))
    bAux = np.zeros((N,1))
    
    paE = []
    paD = []
    aE1 = 1
    aD1 = n
    for i in range(n):
        aE = aE1 + n*i 
        aD = aD1 + n*i
 
        paE.append(aE)
        paD.append(aD)
        
    
    for i in range(N):    
        if (i + 1) <= n:
            b[i] = b0
            
        if (i+1)<= N and (i+1)>(N-n):
            b[i] = c0
 
        if ((i+1) in paE) and ((i+1) != paE[0]) and ((i+1) != paE[-1]):
            b[i] = e0
        if ((i+1) in paD) and ((i+1) != paD[0]) and ((i+1) != paD[-1]): 
            b[i] = d0
            
        if (i+1) == paE[0] or (i+1) == paE[-1]:
            bAux[i] = e0
        if (i+1) == paD[0] or (i+1) == paD[-1]: 
            bAux[i] = d0
    return -(b + bAux)
    
def resolveSistema(A,b):
	X = np.linalg.solve(A,b)
	return X

def Jacobi(X0,C,g,epsilon):
    for i in range(1000):
        X1 = g - np.dot(C,X0)
    
        erro = (np.linalg.norm(X1-X0,np.inf))/np.linalg.norm(X1,np.inf)
        #print i,erro
        if erro < epsilon:
            return X1,erro        
        else: 
            X0 = X1
    
    print ("nao convergiu")
    
 
 ######################################
### GRAFICOS
######################################

def geraGraficos(X1,n):

	x = range(n+2)
	y = range(n+2)

	X,Y = np.meshgrid(x,y)

	Zinter = np.reshape(X1,(n,n))

	Z = np.zeros((n+2,n+2))

	for i in range(n+2):
	    if (i != 0) and (i != n+1):
        		Z[i,0] = 20
        		Z[i,n+1] = 20
        		
        		Z[0,i]   = 10
        		Z[n+1,i] = 40


	for i in range(n):
	    for j in range(n):
            Z[i+1,j+1]=Zinter[i,j]

        	fig = plt.figure(figsize=plt.figaspect(0.5))

	ax0 = fig.add_subplot(1, 2, 1, projection='3d')
	ax1 = fig.add_subplot(1, 2, 2, projection='3d')

	for i in range(n+2):
	    for j in range(n+2):
		x = X[i,j]
		y = Y[i,j]
		z = Z[i,j]
		
		ax0.scatter(x,y,z)

	ax1.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
		linewidth=0.2, antialiased=False)

	plt.show()
    
    
    

######################################
### DADOS
######################################
n = 4
epsilon = 0.01

M,D = montaMatriz(n)     
b = montaB(n,20,20,40,10)

D1 = np.linalg.inv(D)

C = np.dot(D1,M)
g = np.dot(D1,b)


X0 = np.zeros((n*n,1))

X1 = resolveSistema(M+D,b)
#X1,erro = Jacobi(X0,C,g,epsilon)
    
#print X1, erro 
 
geraGraficos(X1,n)

