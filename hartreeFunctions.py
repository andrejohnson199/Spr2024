#Methods for Hartree-Fock
#Mostly reformatting coupledFunction files
#Andre Johnson
import numpy as np
#from coupledFunctions import *
from sympy.physics.sho import R_nl
from numpy import linalg as LA
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import clebsch_gordan as cg


class coupled:
    def __init__(self,n,l,j,m,t):
        self.n = n
        self.l = l
        self.j = j
        self.m = m
        self.t = t

class twoCoupled:
    def __init__(self,a,b,J,T):
        self.a = a #index for single particle state
        self.b = b #index for single partice state
        self.J = J
        self.T = T

#Filling all single particle states
#Takes in Energy level,Proton/Neutron number
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
                            t=-1/2
                            nljarr = np.append(nljarr,[coupled(n,l,j,m,t)])
                            print(n,l,j,m,t)
                            pFill+=1
                        #Putting a neutron in
                        if nFill < nTot:
                            t=1/2
                            nljarr = np.append(nljarr,[coupled(n,l,j,m,t)])
                            print(n,l,j,m,t)
                            nFill+=1
                        m+=1
                else:
                    for i in range(2):
                        for k in range(int(2*j)+1):
                            #Putting a proton in
                            if pFill < pTot:
                                t=-1/2
                                nljarr = np.append(nljarr,[coupled(n,l,j,m,t)])
                                print(n,l,j,m,t)
                                pFill+=1
                            #Putting a neutron in
                            if nFill < nTot:
                                t=1/2
                                nljarr = np.append(nljarr,[coupled(n,l,j,m,t)])
                                print(n,l,j,m,t)
                                nFill+=1
                            m+=1
                        j-=1
                        m=-j        
    return nljarr

def fillTwoStates(oneStates):
    numSize = oneStates.size
    twoArr = np.array([],dtype=twoCoupled)
    for i in range(numSize):
        a = oneStates[i]
        for k in range(numSize):
            # a and b are the single particle states
            # a and b cannot be the same state(i and k)
            if i!=k:
                b = oneStates[k]
                #need to get the J coupling
                T = abs(a.t+b.t)
                for J in range(int(abs(a.j-b.j)),int(a.j+b.j)+1):
                    #Included isospin in single states
                    twoArr = np.append(twoArr,[twoCoupled(i,k,J,T)])
                    #print(i,k,J,T)
    return twoArr

#New normalization factor
#includes mj
def N2(a,b,J,T):
    s = 0
    #Might have to include ms in if condition, not sure
    if a.n==b.n and a.l==b.l and a.j==b.j and a.m==b.m:
        s = 1
        #J+T is even is not allowed
        if (J+T)%2==0:
            return -1
    return ((1-s*(-1)**(J+T))**(1/2))/(1+s)


#Function for building the Hamiltonian
#Hartree Fock Basis
# def buildHF(twoStates,singleStates):
#     hMat = np.zeros((len(twoStates),len(twoStates)))
#     for i in range(len(hMat)):
#         for j in range(len(hMat)):
#             hMat[i][j] = hElement(twoStates,singleStates)
#     return 0

#Potential element in uncoupled basis
#a,b,c,d are all single particle states in coupled basis
def vUncouple(a,b,c,d):
    #TESTING EFFICIENCY
    #return 0
    #Return value
    vReturn = 0
    # a!=b and c!=d
    if a.n==b.n and a.l==b.l and a.j==b.j and a.m==b.m and a.t==b.t:
        return vReturn
    if c.n==d.n and c.l==d.l and c.j==d.j and c.m==d.m and c.t==d.t:
        return vReturn
    
    #To go from coupled to uncoupled, need to sum over J,Mj,T,Mt
    #Will include selection rules for J,T,parity

    #Parity check first
    if (a.l-c.l)%2 != (b.l-d.l)%2:
        #print("Parity Failed")
        return vReturn
    
    #Isospin check
    if a.t+b.t != c.t+d.t:
        return vReturn
    
    #Summation over J, J for both two particle states must be equal
    for J1 in range(int(abs(a.j-b.j)),int(a.j+b.j)+1):
        for J2 in range(int(abs(c.j-d.j)),int(c.j+d.j)+1):
            if J1 == J2:
                #print("J's are equal: ",J1)

                #Mj loop
                Mj = -J1
                p3 = 0
                for i in range(2*J1+1):
                    #print("Hitting Mj loop: ",Mj)
                    p3 += cg(a.j,b.j,J1,a.m,b.m,Mj)*cg(c.j,d.j,J1,c.m,d.m,Mj)
                    Mj+=1
                #LOOP OVER T
                #DO I NEED BOTH T LOOPS
                #WE KNOW T IS EITHER 0 OR 1
                for T1 in range(int(abs(a.t+b.t))+1):
                    for T2 in range(int(abs(c.t+d.t))+1):
                        if T1 == T2:
                            #print("T's are equal: ",T1)
                            #Loop for all Mj and Mt projections

                            #Normalization check
                            p1 = 0
                            if N2(a,b,J1,T1)!=-1 and N2(c,d,J1,T1)!=-1:
                                p1 = (N2(a,b,J1,T1)*N2(c,d,J1,T1))**(-1)

                            #Coupled Matrix Element
                            p4 = vElement(a,b,c,d,J1,T1)

                            #Loop for all Mt 
                            Mt = -T1
                            p2 = 0
                            for j in range(2*T1+1):
                                #print("Hitting Mt loop: ",Mt)
                                p2 += cg(1/2,1/2,T1,a.t,b.t,Mt)*cg(1/2,1/2,T1,c.t,d.t,Mt)
                                Mt+=1
                            vReturn += p1*p2*p3*p4
    return vReturn

#Hamiltonian
#Think I only need single states
#Returns H and array of uncoupled potential matrix elements
#This is to save time
# def H(singleStates):
#     #Matrix of the singleStates
#     hReturn = np.zeros((singleStates.size,singleStates.size))
#     #Array for potential values
#     vReturn = np.zeros(singleStates.size**(4))
#     i=0
#     for row in range(singleStates.size):
#         for col in range(singleStates.size):
#             #Getting kinetic part
#             if row==col:
#                 hReturn[row][col] = 2*singleStates[row].n + singleStates[row].l

#             #HF wavefunctions and potential

#             #Skipping a bit to the summation over j1 and j2
#             #I believe this summation is over only the states below the fermi surface
#             #But we are using 16 O, so all states are used
#             #Might throw errors if j1=row and j2=col
#             for j1 in range(singleStates.size):
#                 for j2 in range(singleStates.size):
#                     p1 = vUncouple(singleStates[row],singleStates[j1],singleStates[col],singleStates[j2])
#                     # print(p1)
#                     vReturn[i] = p1
#                     i+=1
                    
#                     #Initial run, use HO radial wavefunctions
#                     hReturn[row][col] += R_nl(singleStates[j1].n,singleStates[j1].l,1,1)*p1*R_nl(singleStates[j2].n,singleStates[j2].l,1,1)
#     return hReturn, vReturn

#New initial H building method
def newH(singleStates):
    #C matrix is for beta and j values
    #Initially it is the identity
    #beta is column, j is row
    printSize = 0
    numPrint = 16
    C = np.identity(singleStates.size)

    #Matrix of the singleStates
    hReturn = np.zeros((singleStates.size,singleStates.size))
    #Array for potential values
    

    #Row and column indeces for hReturn
    for j in range(singleStates.size):
        for k in range(singleStates.size):
            #Getting kinetic part
            if j==k:
                hReturn[j][k] = 2*singleStates[j].n + singleStates[j].l

            #Loop to cycle through beta terms(occupied orbits)
            #Hard coded 16, will fix later
            for b in range(16):
                #Loops to cycle through all orbits
                for j1 in range(singleStates.size):
                    for j2 in range(singleStates.size):
                        p1 = C[b][j1]
                        p3 = C[j2][b]
                        #Most of these will be zero
                        if p1!=0 and p3!=0:
                            p2 = vUncouple(singleStates[j],singleStates[j1],singleStates[k],singleStates[j2])
                            if printSize < numPrint:
                                #print('Uncoupled V: ',j,j1,k,j2,':',p2)
                                printSize+=1

                            hReturn[j][k] += p1*p2*p3
                            #Printing off-diagonal terms
                            #There shouldnt be any
                            if j!=k and hReturn[j][k]!=0:
                                print('Off-Diagonal V: ',j,j1,k,j2,':',p2)


    return hReturn

#NEW 7/15
#Adding in J and T check
#if a and b are equal
def W(a,b,J,T):
    w = 0
    if a.n==b.n and a.l==b.l and a.j==b.j and (-1)**(J+T+1)==1:
       w = 1
    return w

#This is the simplified potential
def K(a,b,c,d):
    #Set A_T equal to 1 to match book
    A_T = 1
    return -(1/4)*A_T*(-1)**(a.n+b.n+c.n+d.n)

def N(a,b,J,T):
    s = 0
    #Might have to include ms in if condition, not sure
    if a.n==b.n and a.l==b.l and a.j==b.j:
        s = 1
    return ((1-s*(-1)**(J+T))**(1/2))/(1+s)

def jHat(a):
    return (2*a.j+1)**(1/2)

def vElement(a,b,c,d,J,T):
    #One long expression broken up in pieces
    p1 = K(a,b,c,d)*N(a,b,J,T)*N(c,d,J,T)
    p2 = (1+(-1)**(a.l+b.l+c.l+d.l))*jHat(a)*jHat(b)*jHat(c)*jHat(d)
    #Next two go p3-p4
    p3 = (1+(-1)**T)*wigner_3j(a.j,b.j,J,1/2,1/2,-1)*wigner_3j(c.j,d.j,J,1/2,1/2,-1)
    p4 = (-1)**(a.l+c.l+b.j+d.j)*(1-(-1)**(c.l+d.l+J+T))*wigner_3j(a.j,b.j,J,1/2,-1/2,0)*wigner_3j(c.j,d.j,J,1/2,-1/2,0)
    p5 = W(a,b,J,T)*W(c,d,J,T)
    return p1*p2*(p3-p4)*p5


#Overloading H
#New parameters are eigenvectors
#Since the original H happens before this, should store vUncouple values somewhere
# def H2(singleStates,eVectors,pots):
#     #Matrix of the singleStates
#     hReturn = np.zeros((singleStates.size,singleStates.size))
#     i=0
#     for row in range(singleStates.size):
#         for col in range(singleStates.size):
#             #Getting kinetic part
#             if row==col:
#                 hReturn[row][col] = 2*singleStates[row].n + singleStates[row].l

#             #HF wavefunctions and potential

#             #Skipping a bit to the summation over j1 and j2
#             #I believe this summation is over only the states below the fermi surface
#             #But we are using 16 O, so all states are used
#             #Might throw errors if j1=row and j2=col
#             for j1 in range(singleStates.size):
#                 for j2 in range(singleStates.size):
#                     p1 = pots[i]
#                     i+=1
#                     # p1 = vUncouple(singleStates[row],singleStates[j1],singleStates[col],singleStates[j2])
                    
#                     #Initial run, use HO radial wavefunctions
#                     #hReturn[row][col] += R_nl(singleStates[j1].n,singleStates[j1].l,1,1)*p1*R_nl(singleStates[j2].n,singleStates[j2].l,1,1)

#                     #Need to get correct eigenvector and its transpose
#                     colVec = eVectors[:,j2]
#                     rowVec = eVectors[:,j1]
#                     #Build Storage for p1
#                     hReturn[row][col] += np.vdot(rowVec,p1*colVec)
#     return hReturn

#New H2 process
def newH2(ham,C,singleStates):
    #The C matrix is eVec from newHF
    #C matrix is for beta and j values
    #rows and cols are equal
    rows,cols = ham.shape

    #Matrix of the singleStates
    hReturn = np.zeros((rows,cols))

    #Row and column indeces for hReturn
    for j in range(rows):
        for k in range(cols):
            #Getting kinetic part
            if j==k:
                hReturn[j][k] = 2*singleStates[j].n + singleStates[j].l

            #Loop to cycle through beta terms(occupied orbits)
            #Hard coded 16, will fix later
            for b in range(16):
                #Loops to cycle through all orbits
                for j1 in range(singleStates.size):
                    for j2 in range(singleStates.size):
                        p1 = C[b][j1]
                        p3 = C[j2][b]
                        #Most of these will be zero
                        if p1!=0 and p3!=0:
                            p2 = vUncouple(singleStates[j],singleStates[j1],singleStates[k],singleStates[j2])
                            hReturn[j][k] += p1*p2*p3
    return hReturn



#HF process
#Diagonalize and record eigenvalues and eigenkets
# def HF(ham,singleStates,pots):
#     runs = 0
#     eVal, eVec = LA.eig(ham)

#     #Not sure how important eVal is for 16 O, as this is full
#     #Need to rerun Hamiltonian calculation with new eigenvectors
#     newH = H2(singleStates,eVec,pots)
#     eVal2, eVec2 = LA.eig(newH)
#     while(not np.allclose(eVal,eVal2)):
#         eVal, eVec = eVal2, eVec2
#         eVal2, eVec2 = LA.eig(H2(singleStates,eVec,pots))
#         runs+=1
    
#     print('Number of runs: ',runs)
#     return eVal, eVec

#new HF method process for new H process
def newHF(ham,singleStates):
    runs = 1
    eVal, eVec = LA.eig(ham)

    #Rebuild Hamiltonian
    nextH = newH2(ham,eVec,singleStates)
    eVal2, eVec2 = LA.eig(nextH)
    while(not np.allclose(eVal,eVal2)):
        eVal, eVec = eVal2, eVec2
        eVal2, eVec2 = LA.eig(newH2(ham,eVec,singleStates))
        runs+=1

    print('Number of runs: ', runs)
    return eVal, eVec
