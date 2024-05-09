#I will store all my functions in this file then call them in the main file
#Andre Johnson
import numpy as np


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
def fillTwoStates(numSingle):
    twoArr = np.ones((numSingle**2,2),dtype=int)
    row = 0
    for i in range(numSingle):
        for j in range(numSingle):
            twoArr[row][0] = i
            twoArr[row][1] = j
            row+=1
    return twoArr
