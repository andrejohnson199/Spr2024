#Hartree-Fock Program
#Will make use of surface delta potential
#Andre Johnson

#THIS IS A ROUGH DRAFT
#WILL SEPARATE INTO MAIN AND FUNCTION FILE LATER

from hartreeFunctions import *

#Need to know what element
#Starting with 16 O
numP = 8
numN = 8
#Including isospin in single particle would help this
#I know this is 1 because 16O fills to this level
#But how would computer know?
N=1
#Single Particle states
#Fills Os and Op states
singleStates = fillSingleStates(N,numP,numN)
print(singleStates.size)
#Need new fill SingleStates to include m projection and isospin
