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

#Potential Matrix
vMatrix = fillPotential(twoStates,singleStates)