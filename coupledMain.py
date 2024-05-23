#Main file for the surface delta potential
#This time using the coupled basis
#Andre Johnson

import numpy as np
from coupledFunctions import *

#The first step is the get the max energy(N) level
#The states will fill up to this level
Nmax = 1

#Get all single particle states
singleStates = fillSingleStates(Nmax)


#Getting all two particle states
twoStates = fillTwoStates(singleStates)

# #testing vElement
# #All states in 0s with J=0 and T=1, v=-1
a = singleStates[0]
b = singleStates[1]
c = singleStates[0]
d = singleStates[1]
vTest = vElement(a,b,c,d,0,1)
print(vTest)

#Potential Matrix
vMatrix = fillPotential(twoStates,singleStates)