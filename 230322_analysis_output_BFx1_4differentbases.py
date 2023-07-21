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


def date_now(now):
    year = '{:02d}'.format(now.year)
    month = f'{now.month}' 
    day = f'{now.day}' 
    string_date = year + month + day
    print(string_date)
    return string_date

def graph_equilibration_check(length_polymer, error_prob,delta_g_bb, delta_g_tt, model):
    plt.plot(length_polymer, error_prob, label = "$\Delta G_{bb}$ = %s, $\Delta G_{tt}$ = %s"%(delta_g_bb, delta_g_tt))
    plt.title("Energetic: model: %s, Check when error is equilibrated as function of the length of the mRNA grown for different $\Delta G_{bb}$ and $\Delta G_{tt}$" %model)
    plt.xlabel("length mRNA")
    plt.ylabel("$\epsilon$")
    plt.legend(bbox_to_anchor=(1.04, 1))
    plt.subplots_adjust(right=0.7)
    return

def graph_error_g_bb(L_delta_g_bb, L_epsilon, model_title, bases, L_std_high, L_std_low, color_def):
    plt.plot(L_delta_g_bb, L_epsilon, label = "bases = %s" %(bases), color= color_def)
    plt.fill_between(L_delta_g_bb, L_std_low, L_std_high, color= color_def, alpha=0.1)
    plt.title("%s: Effect $\Delta G_{pol}$ on error probability of the copy polymer"%(model_title))
    plt.ylim(0,0.6)
    plt.xlabel("Backbone strength $\Delta G_{pol}$")
    plt.ylabel("Error probability $\epsilon$")
    plt.legend(bbox_to_anchor=(1.04, 1))
    plt.subplots_adjust(right=0.7)
    return




def make_model_names(bases_combies, model, models_title):
    L_models= []
    L_model_titles = []
    for bases in bases_combies:
        L_models.append(model + bases)
        string = models_title + " " + bases
        L_models_titles.append(string)

    return L_models, L_models_titles


L_delta_g_pol_label = ["-ln2", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
L_delta_g_pol_number = [-np.log(2), 0,1,2,3,4,5,6,7,8,9,10]



all_steps = True
DIR = "230323output"
model = "muriel_big_abs_BF_x1"
model_title = "BF x=1"
base_combinations = [ "GCXK","iGiCXK","aATGC", "deltabetaXK","aATXK"]
ratio ='1111'
# number_steps = 0
# MX_models = ["x=$e^{\Delta G_{TT} +5}$, Model B","x=$e^{\Delta G_{TT} +5}$, Model F", "x=1, Model B", "x=1, Model F" ]

# L_models, L_models_titles = make_model_names(base_combinations, model, model_title)
L_models = [model]
L_models_titles = [model_title]
date = date_now(datetime.datetime.now())
date_options = ['2023323', '2023327']


#define colors
colormap = plt.cm.plasma #nipy_spectral, Set1,Paired   
color_list = [colormap(i) for i in np.linspace(0, 1,100)]

print(L_models)



# line_styles = ['dashed', 'solid']

# L_models = L_models[i:i+3] #for same M and X diffent models
line_styles = ['solid', 'dashdot', 'dashed', 'dotted', (5, (10, 3)), (0, (3, 5, 1, 5))]
print(L_models)


for k, model in enumerate(L_models):
    for j, base in enumerate(base_combinations):
        L_epsilon_mean = []
        L_epsilon_std_high = []
        L_epsilon_std_low = []
        L_delta_g_bb_graph = []
        for i,delta_g_pol in enumerate(L_delta_g_pol_number):     

            for date in date_options:
                string = '/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/' + DIR + '/' +date + 'model_' + model + base +ratio + '_delta_g_pol_' + L_delta_g_pol_label[i]+'_allsteps.csv'
                file_exists = os.path.exists(string)
                print(string)
                if file_exists == True:
                    data = pd.read_csv(string)
                    epsilon_mean = data['error probability'].mean()
                    epsilon_std = data['error probability'].std()
                    L_epsilon_mean.append(epsilon_mean)
                    L_epsilon_std_high.append(epsilon_mean + epsilon_std)
                    L_epsilon_std_low.append(epsilon_mean - epsilon_std)
                    L_delta_g_bb_graph.append(delta_g_pol)
                    
                    


                else:
                    print("file does not exist")


        #define error as all mismatches and plot the error as function of delta G_bb
        graph_error_g_bb(L_delta_g_bb_graph, L_epsilon_mean, L_models_titles[k], base, L_epsilon_std_high, L_epsilon_std_low, color_list[j*15])
plt.show()


