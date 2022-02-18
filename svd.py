

import numpy as np

from scipy import linalg


#Biblioteca para manipulacao de 
#imagens
from PIL import Image

im1 = Image.open("rosto1.jpg")

M1 = np.asarray(im1.convert('L'))

U, s, Vh = linalg.svd(M1, full_matrices=False)

r = len(s)
S = np.zeros((r,r))


c = 50

for i in range(c):
    S[i,i] = s[i]


M2 = np.dot(U,np.dot(S,Vh))


print ("Reconstrucao com ",100.0*c/r,"%")
print ("com erro de ",linalg.norm((M2-M1))/linalg.norm((M1)))

im1 = Image.fromarray(M1)
im2 = Image.fromarray(M2)

im1.show()
im2.show()
