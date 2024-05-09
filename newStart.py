#Here I am reformatting Harmonic Oscillator Perturbation
#The goal is to make a cleaner code
#Andre Johnson

import numpy as np
import sympy as sympy
from sympy.physics.sho import E_nl
from sympy.physics.sho import R_nl
from newFunctions import *

#The first step is the get the max energy(N) level
#The states will fill up to this level
Nmax = 1

#Get all single particle states
singleStates = fillSingleStates(Nmax)

#Getting two particle states
#array whose entries are the index for the single particle array
twoStates = fillTwoStates(singleStates.size)

