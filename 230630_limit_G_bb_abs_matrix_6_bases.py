"""
21-11-22, Muriel Louman, AMOLF

analysis 221121_grow_mRNA_wr_energetic_Jennys_way_test_rates
"""

#import needed packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import datetime
import os.path
from random import randint

def end_pur1(E1, E2, E3, E4,E5,E6):
    return 1- 1/(6*(1 + E1)*(1 - E1/(6*(1 + E1)) - E2/(6*(1 + E2)) - E3/(
   6*(1 + E3)) - E4/(6*(1 + E4)) - E5/(6*(1 + E5)) - E6/(6*(1 + E6))))

def end_py1(E1, E2, E3, E4,E5,E6):
    return 1- 1/(6*(1 + E2)*(1 - E1/(6*(1 + E1)) - E2/(6*(1 + E2)) - E3/(
   6*(1 + E3)) - E4/(6*(1 + E4)) - E5/(6*(1 + E5)) - E6/(6*(1 + E6))))

def end_pur2(E1, E2, E3, E4,E5,E6):
    return 1- 1/(6*(1 + E3)*(1 - E1/(6*(1 + E1)) - E2/(6*(1 + E2)) - E3/(
   6*(1 + E3)) - E4/(6*(1 + E4)) - E5/(6*(1 + E5)) - E6/(6*(1 + E6))))

def end_py2(E1, E2, E3, E4,E5,E6):
    return 1 -1/(6*(1 + E4)*(1 - E1/(6*(1 + E1)) - E2/(6*(1 + E2)) - E3/(
   6*(1 + E3)) - E4/(6*(1 + E4)) - E5/(6*(1 + E5)) - E6/(6*(1 + E6))))

def end_pur3(E1, E2, E3, E4,E5,E6):
    return 1- 1/(6*(1 + E5)*(1 - E1/(6*(1 + E1)) - E2/(6*(1 + E2)) - E3/(
   6*(1 + E3)) - E4/(6*(1 + E4)) - E5/(6*(1 + E5)) - E6/(6*(1 + E6))))

def end_py3(E1, E2, E3, E4,E5,E6):
    return 1 -1/(6*(1 + E6)*(1 - E1/(6*(1 + E1)) - E2/(6*(1 + E2)) - E3/(
   6*(1 + E3)) - E4/(6*(1 + E4)) - E5/(6*(1 + E5)) - E6/(6*(1 + E6))))

L_rep = [0,1,2,3,4,5,6,20,100]
energies = [[+3,-3,-1,+1,1,-1],[+3,-3,-1,+1,-1,+1],[+3,-3,-1,+1,-3,3]]
L_E = []
L_nucleotides = ["aTGCLambdaSigma", "aATGCiGiC", "aATGCXK"]

for rep in L_rep:
    print("rep", rep)
    for energie in energies:
        print("energies", energie)
        L_Ea= []
        for i,E in enumerate(energie):
            if i ==0 or i==2 or i==4:
                E = E+rep
            print(E)
            E = np.exp(E)
            L_Ea.append(E)
        L_E.append(L_Ea)

errors = []
for config in L_E:
    error = end_py1(config[0], config[1], config[2], config[3],config[4],config[5])
    errors.append(error)

L_errors = []
for i in range(len(L_rep)):
    print("cost = %s, errors for aATGCLambdaSIgma, aATGCiGiC, aATGCXK:"%L_rep[i])
    print(errors[0 + 3*i], errors[1+3*i], errors[2+3*i])
    e_1 =  errors[0 + 3*i]
    e_2 = errors[1 + 3*i]
    e_3 = errors[2 + 3*i]
    L_errors.append([e_1,e_2,e_3])

print(L_errors)