#Functions for finding matrix elements of surface delta potential
#This time in the coupled basis!
#Andre Johnson

import numpy as np
from tabulate import tabulate
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import clebsch_gordan as cg


#Maybe add isospin later
class coupled:
    def __init__(self,n,l,j):
        self.n = n
        self.l = l
        self.j = j

class twoCoupled:
    def __init__(self,a,b,J,T):
        self.a = a #index for single particle state
        self.b = b #index for single partice state
        self.J = J
        self.T = T

#function to fill all single particle states
#returns an array of all single particle states
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

def fillTwoStates(oneStates):
    numSize = oneStates.size
    twoArr = np.array([],dtype=twoCoupled)
    for i in range(numSize):
        a = oneStates[i]
        for k in range(numSize):
            #WITHOUT MS
            # a and b are the single particle states
            b = oneStates[k]
            #need to get the J coupling
            for J in range(int(abs(a.j-b.j)),int(a.j+b.j)+1):
                #T=0,1
                twoArr = np.append(twoArr,[twoCoupled(i,k,J,0)])
                twoArr = np.append(twoArr,[twoCoupled(i,k,J,1)])
    return twoArr

def fillPotential(two,one):
    vMat = np.zeros((len(two),len(two)))
    #Table
    table = np.array(['abcd','J','T','V'])
    aname = ''
    bname = ''
    cname = ''
    dname = ''
    for row in range(len(vMat)):
        a = one[two[row].a]
        b = one[two[row].b]
        J1 = two[row].J
        T1 = two[row].T

        #Printing name
        if a.l==0:
            aname = 1
        elif a.j==1.5:
            aname = 2
        else:
            aname = 3
        if b.l==0:
            bname = 1
        elif b.j==1.5:
            bname = 2
        else:
            bname = 3

        for col in range(len(vMat)):
            c = one[two[col].a]
            d = one[two[col].b]
            J2 = two[col].J
            T2 = two[col].T

            #Printing name
            if c.l==0:
                cname = 1
            elif c.j==1.5:
                cname = 2
            else:
                cname = 3
            if d.l==0:
                dname = 1
            elif d.j==1.5:
                dname = 2
            else:
                dname = 3

            #Parity, J, and T need to match
            #Nontrivial elements
            #Not printing the zeros
            if (a.l-c.l)%2 == (b.l-d.l)%2 and J1 == J2 and T1 == T2:
                vMat[row][col] = vElement(a,b,c,d,J1,T1)

                #Do not need to print redundant states
                if abs(vMat[row][col])>0 and cname<=dname and aname<=bname: 
                    #Adding to the table
                    tNum = 1000*aname+100*bname+10*cname+dname
                    table = np.vstack([table,[tNum,J1,T1,vMat[row][col]]])
    #Print to terminal
    print(tabulate(table,headers='firstrow'))
    print('1=0s1/2, 2=0p3/2, 3=0p1/2')
    # #Print to output file
    # #Uncomment if you want
    # with open("output.txt", "a") as f:
    #     print(tabulate(table,headers='firstrow'),file=f)
    #     print('1=0s1/2, 2=0p3/2, 3=0p1/2',file=f)
    return vMat

#NEW 7/15
#Adding in J and T check
#if a and b are equal
#def W(a,b,J,T):
#    if a.n==b.n and a.l==b.l and a.j==b.j and a.t==b.t

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
    return p1*p2*(p3-p4)
