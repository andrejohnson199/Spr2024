#I will store all my functions in this file then call them in the main file
#Andre Johnson
import numpy as np
from sympy.physics.sho import R_nl


#Class to store all the quantum numbers for each state
#I will use ms=0,1 for up/down
class nlmbasis:
    def __init__(self,n,l,m,ms):
        self.n=n
        self.l=l
        self.m=m
        self.ms=ms

    def print_nlm(self):
        print(self.n,self.l,self.m,self.ms)

#function to fill all single particle states
#returns an array of all single particle states
def fillSingleStates(N):
    nlmarr = np.array([],dtype=nlmbasis)
    for n in range(int(N/2+1)):
        for l in range(N+1):
            if 2*n + l <= N:
                for m in range(-l,l+1):
                   nlmarr = np.append(nlmarr,[nlmbasis(n,l,m,0)])
                   nlmarr = np.append(nlmarr,[nlmbasis(n,l,m,1)])
    return nlmarr

#function to get all two particle states
#takes single state array size as argument
#Array has two columns and K(K-1) rows
#K is number of single states
def fillTwoStates(numSingle):
    twoArr = np.ones((numSingle*(numSingle-1),2),dtype=int)
    row = 0
    for i in range(numSingle):
        for j in range(numSingle):
            if i!=j:
                twoArr[row][0] = i
                twoArr[row][1] = j
                row+=1
    return twoArr

#function to make the potential matrix
#using surface delta potential
def fillPotential(two,one,v,r):
    vMat = np.zeros((len(two),len(two)))
    #Printing variable
    i=0
    for row in range(len(vMat)):
        bra1 = one[two[row][0]]
        bra2 = one[two[row][1]]
        for col in range(len(vMat)):
            ket1 = one[two[col][0]]
            ket2 = one[two[col][1]]

            #Parity Rule
            if (bra1.l-ket1.l)%2==(bra2.l-ket2.l)%2:
                vMat[row][col] = -v * r**2 * R_nl(bra1.n,bra1.l,1,r) * R_nl(bra2.n,bra2.l,1,r) * R_nl(ket1.n,bra1.l,1,r) * R_nl(ket1.n,bra1.l,1,r)

                #Printing nonzero elements to check
                if i<10:
                    print(vMat[row][col],'Element:',row,col,' Ket1:',ket1.n,ket1.l,ket1.m,ket1.ms, ' Ket2:',ket2.n,ket2.l,ket2.m,ket2.ms, ' Bra1:',bra1.n,bra1.l,bra1.m,bra1.ms, ' Bra2:',bra2.n,bra2.l,bra2.m,bra2.ms,)
                    i+=1

            # #Rules
            # if bra1.l==ket1.l and bra2.l==ket2.l and bra1.m==ket1.m and bra2.m==ket2.m and bra1.ms==ket1.ms and bra2.ms==ket2.ms:
            #     vMat[row][col] = -v * r**2 * R_nl(bra1.n,bra1.l,1,r) * R_nl(bra2.n,bra2.l,1,r) * R_nl(ket1.n,bra1.l,1,r) * R_nl(ket1.n,bra1.l,1,r)

            #     #Printing nonzero elements to check
            #     print(vMat[row][col],'Element:',row,col,' Ket1:',ket1.n,ket1.l,ket1.m,ket1.ms, ' Ket2:',ket2.n,ket2.l,ket2.m,ket2.ms, ' Bra1:',bra1.n,bra1.l,bra1.m,bra1.ms, ' Bra2:',bra2.n,bra2.l,bra2.m,bra2.ms,)
    return vMat
