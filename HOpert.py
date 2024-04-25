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

#This is a different way of doing the nlm 2 state basis
#This will hold integer index values for the appropiate single nlm state
class twoState:
    def __init__(self,i,j):
        self.i = i
        self.j = j


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
                #nlmarr[i].print_nlm()
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
       # nlm2arr[i].print_nlm2()
        i+=1
        
        
print('Done with double states')
#Surface delta potential
#Set R
R = 5
vmat = np.zeros((3*N)**2)
for s in range(nlm2arr.size):
    vmat[s] = (-V)*(1/R**2)*R_nl(nlm2arr[s].n1,nlm2arr[s].l1,1,R)*R_nl(nlm2arr[s].n2,nlm2arr[s].l2,1,R)*(4*np.pi)**2
    # print(vmat[s])
print('Done with V')

#Energy terms, can combine with previous loop later
#I have N as the input parameter so all states share the same N 
#So this loop is not really needed, but I will keep it in case I end up 
#modifying code to use multiple N values
# emat = np.zeros((3*N)**2)
# for s in range(nlm2arr.size):
#     #I have two print statements because E_nl has a 3/2 in its calculation
#     #While the H I was given does not, ASK
#     emat[s] = E_nl(nlm2arr[s].n1,nlm2arr[s].l1,1)+E_nl(nlm2arr[s].n2,nlm2arr[s].l2,1)
#     # print(emat[s])
#     # print(2*nlm2arr[s].n1+nlm2arr[s].l1+2*nlm2arr[s].n2+nlm2arr[s].l2)

#Printing out the energy of each state
#ASK
#Right now the energy from the H is way larger than the energy shifts
#Is it the scale of R_nl? The given R value?
elevel = np.zeros((3*N)**2)
print('Printing Energy Levels')
for s in range(nlm2arr.size):
    elevel[s] = vmat[s] + 2*nlm2arr[s].n1+nlm2arr[s].l1+2*nlm2arr[s].n2+nlm2arr[s].l2
    print(elevel[s],': ',end=''),nlm2arr[s].print_nlm2()

#Printing just the first order shifts
print('Printing First order shifts')
for s in range(nlm2arr.size):
    print(vmat[s],': ',end=''),nlm2arr[s].print_nlm2()






        

