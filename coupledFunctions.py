#Functions for finding matrix elements of surface delta potential
#This time in the coupled basis!
#Andre Johnson

import numpy as np
from tabulate import tabulate
from sympy.physics.wigner import wigner_3j

#Ill right out all the different functions 
#Then make it look neat later

class coupled:
    def __init__(self,n,l,j,ms):
        self.n = n
        self.l = l
        self.j = j
        self.ms = ms

class twoCoupled:
    def __init__(self,a,b,J,T):
        self.a = a #index for single particle state
        self.b = b #index for single partice state
        self.J = J
        self.T = T

#GETTING RID OF MS
#I BELIEVE IT IS JUST DOUBLE THE AMOUNT I HAVE TO STORE

#function to fill all single particle states
#returns an array of all single particle states
def fillSingleStates(N):
    nljarr = np.array([],dtype=coupled)
    for n in range(int(N/2+1)):
        for l in range(N+1):
            if 2*n + l <= N:
                j = l+1/2
                if l == 0:
                    nljarr = np.append(nljarr,[coupled(n,l,j,0)])
                    #nljarr = np.append(nljarr,[coupled(n,l,j,1)])
                else:
                    for i in range(2):
                        nljarr= np.append(nljarr,[coupled(n,l,j,0)])
                        #nljarr = np.append(nljarr,[coupled(n,l,j,1)])
                        j-=1        
    return nljarr

def fillTwoStates(oneStates):
    numSize = oneStates.size
    #Setting T=1, like nucleons
    #T=1
    twoArr = np.array([],dtype=twoCoupled)
    for i in range(numSize):
        a = oneStates[i]
        for k in range(numSize):
            # if i!=k:
            #     # a and b are the single particle states
            #     b = oneStates[k]
            #     #need to get the J coupling
            #     for J in range(int(abs(a.j-b.j)),int(a.j+b.j)+1):
            #         #T=0,1
            #         twoArr = np.append(twoArr,[twoCoupled(i,k,J,0)])
            #         twoArr = np.append(twoArr,[twoCoupled(i,k,J,1)])

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
    #Printing index
    i = 0
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
            aname = '0s1/2'
            aname = 1
        elif a.j==1.5:
            aname = '0p3/2'
            aname = 2
        else:
            aname = '0p1/2'
            aname = 3
        if b.l==0:
            bname = '0s1/2'
            bname = 1
        elif b.j==1.5:
            bname = '0p3/2'
            bname = 2
        else:
            bname = '0p1/2'
            bname = 3
        for col in range(len(vMat)):
            c = one[two[col].a]
            d = one[two[col].b]
            J2 = two[col].J
            T2 = two[col].T

            #Printing name
            if c.l==0:
                cname = '0s1/2'
                cname = 1
            elif c.j==1.5:
                cname = '0p3/2'
                cname = 2
            else:
                cname = '0p1/2'
                cname = 3
            if d.l==0:
                dname = '0s1/2'
                dname = 1
            elif d.j==1.5:
                dname = '0p3/2'
                dname = 2
            else:
                dname = '0p1/2'
                dname = 3

            #Parity, J, and T need to match
            #Nontrivial elements
            #Not printing the zeros
            if (a.l-c.l)%2 == (b.l-d.l)%2 and J1 == J2 and T1 == T2:
                vMat[row][col] = vElement(a,b,c,d,J1,T1)

                #Do not need to print redundant states
                if abs(vMat[row][col])>0 and cname<=dname and aname<=bname: 
                    #print(vMat[row][col],'Element: ',row,col,' a: ',a.n,a.l,a.j,' b: ',b.n,b.l,b.j,' c: ',c.n,c.l,c.j,' d: ',d.n,d.l,d.j,' J: ',J1,' T: ',T1)

                    #New Printing
                    #print(vMat[row][col],'Element:',row,col,' abcd: ',orbits[0],orbits[1],orbits[2],orbits[3])
                    #print(vMat[row][col],'Element:',row,col,' abcd: ',aname,bname,cname,dname,' J: ',J1,' T: ',T1)
                    #print(aname,bname,cname,dname, ' ',J1,T1,' ',vMat[row][col])

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