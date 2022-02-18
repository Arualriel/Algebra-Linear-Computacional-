#Biblioteca Numerica
import numpy as np
from scipy import linalg



#Funcao que calcula o ERRO RELATIVO entre X(k+1) e X(k).
def erroRelativo(X0,X1):
    n1 = np.linalg.norm(X1-X0) #Calcula || X(k+1) - X(k) ||
    #n2 = np.linalg.norm(X1,np.inf)    #Calcula || X(k+1)|| (ambos com a norma infinito)
    erro = n1#n1/n2                   #erro = ||X(k+1) - X(k)||/||X(k+1)||
    
    #erro = np.linalg.norm(X1-X0)
    
    return erro

##########################################
### Montando as Matrizes L, D e R 
### a partir da matriz A.
##########################################

#Funcao para montar a matriz R.
def montaMatrizR(A):
        
    n,m = np.shape(A) #ordem da matriz A
    R = np.zeros((n,m))    
    
    for i in range(n):
        for j in range(i,m):
            if i != j :
                R[i,j] = A[i,j]
    
    return R
    
#Funcao para montar a matriz L.
def montaMatrizL(A):
    
    n,m = np.shape(A) #ordem da matriz A
    L = np.zeros((n,m))
    
    for j in range(m):
        for i in range (j,n):
            if j != i:
                L[i,j] = A[i,j]
    return L

#Funcao para montar a matriz D.    
def montaMatrizD(A):

    n,m = np.shape(A) #ordem da matriz A
    D = np.zeros((n,m))    
    
    for i in range(n):
        D[i,i] = A[i,i]
    
    return D
    
    
    
    
    
    
#######################################    
####    Metodos Iterativos
#######################################    

def metodoJacobi(A,b,X0,epsilon,N):
    # A = L + D + R
    L = montaMatrizL(A)
    D = montaMatrizD(A)
    R = montaMatrizR(A)
    
    D1 = np.linalg.inv(D) #matriz inversa da D    
    
    #X0 = X(k) e X1 = X(k+1) 
    C = np.dot(D1,-(L+R))
    print ("C = \n",C)
    
    for i in range(1,N + 1):
        X1 = np.dot(D1,(b - L*X0 - R*X0))
        
        erro = erroRelativo(X1,X0)
        
        ##Para visualizar as iteracoes descomente as duas linhas abaixo:
        print ("X(",i,") = \n",X1)
        print ("Erro Relativo entre X(",i,") e X(",i-1,") :::", erro)
        if erro < epsilon:            
            return i,erro,X1
        else:
            X0 = X1 #Atualizando X0
            
        if i == N :
            #Numero maximo de iteracoes alcancado retorna -1 para 
            #que seja impressa uma mensagem de aviso.
        
            return -1

def metodoSeidel(A,b,X0,epsilon,N):
    # A = L + D + R
    L = montaMatrizL(A)
    D = montaMatrizD(A)
    R = montaMatrizR(A)
    
    LD1 = np.linalg.inv(L+D) #matriz inversa da L + D    
    
    C = np.dot(LD1,-R)
    print ("C = \n",C)
    
    #X0 = X(k) e X1 = X(k+1)
    

    for i in range(1,N + 1):
        X1 = np.dot(LD1,(b - R*X0))
        
        erro = erroRelativo(X1,X0)
        
        ##Para visualizar as iteracoes descomente as duas linhas abaixo:
        print ("X(",i,") = \n",X1)
        print ("Erro Relativo entre X(",i,") e X(",i-1,") :::", erro)
        
        if erro < epsilon:            
            return i,erro,X1
        else:
            X0 = X1 #Atualizando X0
            
        if i == N :
            #Numero maximo de iteracoes alcancado retorna -1 para 
            #que seja impressa uma mensagem de aviso.
        
            return -1






#######################################    
####    Aplicando os Metodos.
#######################################    

# Dado o sistema na forma Ax = b :

#matriz quadrada
#A = np.matrix([[5,1,1],[3,4,1],[3,3,6]])
#A = np.matrix([[2,5],[3,1]])
#A = np.matrix([[11,2],[2,11]])


B = np.matrix([[-0.00019069, -0.00196143],[-0.00196143, -0.00827477]])
I=np.identity(2,float)
K=0.5*I-B
Q,A = linalg.qr(K)
#matriz coluna
#b = np.matrix([[5],[6],[0]]) 
#b = np.matrix([[13],[13]]) 
b=np.zeros((2,1))
b[0,0],b[1,0]=0.08024518, 0.1611033



#precisao epsilon:
epsilon = 10.0**(-8.0)

#Aproximacao inicial X0
#X0 = np.matrix([[0],[0],[0]]) 
X0 = np.matrix([[1],[1]]) 
 
#Numero maximo de iteracoes
N = 100 


print ("***********************************")
print ("  Metodo Iterativo de Gauss-Jacobi ")
print ("***********************************")
 
resultadoJacobi = metodoJacobi(A,b,X0,epsilon,N)

if resultadoJacobi == -1:
    print ("Numero maximo de iteracoes alcancado!")
    print ("A sequencia nao convergiu para as ",N," iteracoes.")
else: 
    iteracaoJ   = resultadoJacobi[0]
    erroJacobi  = resultadoJacobi[1]
    Xjacobi     = resultadoJacobi[2]
    #Xjacobi=Xjacobi*(1/(linalg.norm(Xjacobi)))
    print ("Convergiu com ", iteracaoJ," iteracoes.")
    print ("solucao aproximada X (Gauss-Jacobi):: \n",Xjacobi)
    print ("erro Relativo (Gauss-Jacobi):::", erroJacobi) 
    
    
    
    
    
print ("***********************************")
print ("  Metodo Iterativo de Gauss-Seidel ")
print ("***********************************")
 
resultadoSeidel = metodoSeidel(A,b,X0,epsilon,N)

if resultadoSeidel == -1:
    print ("Numero maximo de iteracoes alcancado!")
    print ("A sequencia nao convergiu para as ",N," iteracoes.")
else: 
    iteracaoS  = resultadoSeidel[0]
    erroSeidel = resultadoSeidel[1]
    Xseidel    = resultadoSeidel[2]
    #Xseidel=Xseidel/(linalg.norm(Xseidel))
    print ("Convergiu com ", iteracaoS," iteracoes.")
    print ("solucao aproximada X (Gauss-Seidel):: \n",Xseidel)
    print ("erro Relativo (Gauss-Seidel):::", erroSeidel)