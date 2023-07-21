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


def incorrect_prob(x, Ee):
    return x/(2*(Ee + x)*(1 - 1/(2*(1 + x)) - Ee/(2*(Ee + x))))


L_g_tt = [0,2,4,6,8,10]


errors = []
for g_tt in L_g_tt:
    exp_g_tt = np.exp(g_tt)
    x = np.exp(g_tt + 5)
    error = incorrect_prob(x, exp_g_tt)
    errors.append(error)



print(errors)

