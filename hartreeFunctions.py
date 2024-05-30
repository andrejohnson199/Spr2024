#Methods for Hartree-Fock
#Mostly reformatting coupledFunction files
#Andre Johnson
import numpy as np

class coupled:
    def __init__(self,n,l,j,m,t):
        self.n = n
        self.l = l
        self.j = j
        self.m = m
        self.t = t

#function to fill all single particle states
#returns an array of all single particle states
#input is number of protons and neutrons
def fillSingleStates(N):
    nljarr = np.array([],dtype=coupled)
    for n in range(int(N/2+1)):
        for l in range(N+1):
            if 2*n + l <= N:
                j = l+1/2
                if l == 0:
                    nljarr = np.append(nljarr,[coupled(n,l,j)])
                else:
                    for i in range(2):
                        nljarr= np.append(nljarr,[coupled(n,l,j)])
                        j-=1        
    return nljarr

#Another method
def fillSingleStates(N,pTot,nTot):
    pFill = 0
    nFill = 0
    nljarr = np.array([],dtype=coupled)
    for n in range(int(N/2+1)):
        for l in range(N+1):
            if 2*n + l <= N:
                j = l+1/2
                m=-j
                if l == 0:
                    for k in range(int(2*j)+1):
                        #Putting a proton in
                        if pFill < pTot:
                            t=1/2
                            nljarr = np.append(nljarr,[coupled(n,l,j,m,t)])
                            print(n,l,j,m,t)
                            pFill+=1
                        #Putting a neutron in
                        if nFill < nTot:
                            t=-1/2
                            nljarr = np.append(nljarr,[coupled(n,l,j,m,t)])
                            print(n,l,j,m,t)
                            nFill+=1
                        m+=1
                else:
                    for i in range(2):
                        for k in range(int(2*j)+1):
                            #Putting a proton in
                            if pFill < pTot:
                                t=1/2
                                nljarr = np.append(nljarr,[coupled(n,l,j,m,t)])
                                print(n,l,j,m,t)
                                pFill+=1
                            #Putting a neutron in
                            if nFill < nTot:
                                t=-1/2
                                nljarr = np.append(nljarr,[coupled(n,l,j,m,t)])
                                print(n,l,j,m,t)
                                nFill+=1
                            m+=1
                        j-=1
                        m=-j        
    return nljarr

#New single state fill method
def newSingleFill(P,N):
    numP = 0
    numN = 0
    #Starting in the Os state
    n = 0
    l = 0
    
    nljarr = np.array([],dtype=coupled)
    #Loop to fill states
    #Start with protons
    for i in range(P):
        j = l+1/2
        m = -j
        if l == 0:
            for k in range(int(2*j)+1):
                nljarr = np.append(nljarr,[coupled(n,l,j,m,1/2)])
                m+=1
        else:
            for i in range(2):
                for k in range(int(2*j)+1):
                    nljarr= np.append(nljarr,[coupled(n,l,j,m,1/2)])
                    j-=1