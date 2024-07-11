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
print('Single Particle Orbits: n,l,j,m,tz')
singleStates = fillSingleStates(N,numP,numN)
#Need new fill SingleStates to include m projection and isospin

#Also need two states
#twoStates = fillTwoStates(singleStates)

# #Testing uncouple method
# vTest = vUncouple(singleStates[0],singleStates[2],singleStates[0],singleStates[2])
# print(vTest)
#print(vElement(singleStates[0],singleStates[2],singleStates[0],singleStates[2],0,1))
#print(vElement(singleStates[0],singleStates[2],singleStates[0],singleStates[2],1,1))

#Testing clebsch gordan
#print(cg(1/2,1/2,0,-1/2,1/2,1))

# #Testing Hamiltonian
# Ham, pots = H(singleStates)
# print(Ham)
# #print(pots)
# print('Finished initial H build, now starting HF ')
# # #Run HF
# Val,Vec = HF(Ham,singleStates,pots)
# print(Val)
# print(Vec)

#Testing new H method
Ham = newH(singleStates)
print(Ham)
#Val, Vec = newHF(Ham,singleStates)
#print(Val)
#print(Vec)



#Build the Hamiltonian(HF)
#H = buildHF()
