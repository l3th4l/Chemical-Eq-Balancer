#chemical equation balancer
import numpy as np 
from sympy import *
import re

#eq = 'Fe + H_2O -> Fe_3O_4 + H_2'

eq = input("Enter your equation : \n")

eql = [x for x in eq.split(' ') if x[0].isalpha()]

elems = list(set(re.findall('[A-Z][a-z]*', eq))) #Gets all the elements required in the eq

mvecs = []
for mol in eql:
    els = re.findall('[A-Z][a-z]*[_]?[0-9]?', mol)

    comp = []

    for el in els:
        if '_' in el:
            e, c = el.split('_')
            c = int(c)
        else:
            c = 1
            e = el

        comp.append(np.array([int(x == e)*c for x in elems]))
    mvecs.append(np.sum(comp, axis = 0))

mvecs.append(np.zeros(len(elems)).astype(int)) #creates augmented matrix for the linear system

mvecs = Matrix(np.array(mvecs).T)

coefs = mvecs.rref() #reduces augmented matrix to reduced echolon form
coefl = []
denom_list = []
for coef in coefs[0][:,-2]:
    denom_list.append(fraction(coef)[-1])

l = lcm(denom_list) # gets lcm of denominators of coefficients

for coef in coefs[0][:,-2]:
    coefl.append(np.abs(coef*l)) #multiplies the absolute value of the coefficient with the lcm
coefl.append(l) #adds the free variable * lcm to the coefficient list

print("\nBalenced Molecules : \n")

for mols in list(zip(coefl,eql)):
    print(('%s%s'%(mols[0], mols[1])).replace('1', ''))