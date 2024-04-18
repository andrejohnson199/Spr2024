#Harmonic Oscillator Perturbation with V(r) = V delta(s)
#Andre Johnson

import numpy as np
import sympy as sympy
from sympy.physics.sho import E_nl
from sympy.physics.sho import R_nl

#specify V
V = 5


class nlmbasis:
    def __init__(self,n,l,m):
        self.n=n
        self.l=l
        self.m=m

    def print_nlm(self):
        print(self.n,self.l,self.m)

class nlm2state:
    def __init__(self,nl,l1,m1,n2,l2,m2):
        self.n1 = n1
        self.l1 = l1
        self.m1 = m1
        self.n2 = n2
        self.l2 = l2
        self.m2 = m2

    def print_nlm2(self):
        print(self.n1,self.l1,self.m1,' ',self.n2,self.l2,self.m2)


#Getting occupied states
#Input is N=2n+l
N=2

#Making an array of all the occupied states 
nlmarr = np.zeros(3*N,dtype=nlmbasis)

#array index
i = 0


for n in range(int(N/2+1)):
    for l in range(N+1):
        if 2*n + l == N:
            for m in range(-l,l+1):
                nlmarr[i] = nlmbasis(n,l,m)
                nlmarr[i].print_nlm()
                i+=1

                
print('Done with single state')

#Making the two-state vector
#Brute force method of all states
i = 0
nlm2arr = np.zeros((3*N)**2,dtype=nlm2state)
for n1 in range(nlmarr.size):
    s1 = nlmarr[n1]
    n1 = s1.n
    l1 = s1.l
    m1 = s1.m
    for n2 in range(nlmarr.size):
        s2 = nlmarr[n2]
        nlm2arr[i] = nlm2state(n1,l1,m1,s2.n,s2.l,s2.m)
        nlm2arr[i].print_nlm2()
        i+=1
        
        





        

