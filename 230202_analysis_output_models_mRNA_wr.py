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

def graph_error_g_bb(L_delta_g_bb, L_epsilon, model, delta_g_tt, model_title, L_std_high, L_std_low,gtt_theory, color_def):
    plt.plot(L_delta_g_bb, L_epsilon, label = "$\Delta G_{tt}$ = %s" %(delta_g_tt), color= color_def)
    if model_title[-1] =="1":
        plt.scatter([10], [gtt_theory], color= color_def)
    plt.fill_between(L_delta_g_bb, L_std_low, L_std_high, color= color_def, alpha=0.1)
    plt.title("%s: Effect $\Delta G_{pol}$ on error probability of the copy polymer"%(model_title))
    plt.ylim(0,0.6)
    plt.xlabel("Backbone strength $\Delta G_{pol}$")
    plt.ylabel("Error probability $\epsilon$")
    plt.legend(bbox_to_anchor=(1.04, 1))
    plt.subplots_adjust(right=0.7)
    return

def graph_error_g_bb_multiple_models_per_model(L_delta_g_bb, L_epsilon, model, delta_g_tt, model_title, L_std_high, L_std_low, color_def, line_style):
    plt.plot(L_delta_g_bb, L_epsilon, label = "$\Delta G_{tt}$ = %s, %s" %(delta_g_tt, model_title[4:12]), linestyle=line_style, color= color_def)
    plt.fill_between(L_delta_g_bb, L_std_low, L_std_high, color= color_def, alpha=0.1)
    plt.title("%s: Effect $\Delta G_{pol}$ on error probability of the polymer string" %(model_title[0:4]+model_title[12:]))
    plt.xlabel("Backbone strength $\Delta G_{pol}$")
    plt.ylabel("Error probability $\epsilon$")
    plt.legend(bbox_to_anchor=(1.04, 1))
    plt.subplots_adjust(right=0.7)
    return

def graph_error_g_bb_multiple_models_per_model2(L_delta_g_bb, L_epsilon, model, delta_g_tt, model_title, L_std_high, L_std_low, color_def, line_style):
    plt.plot(L_delta_g_bb, L_epsilon, label = "$\Delta G_{tt}$ = %s, %s" %(delta_g_tt, model_title[0:4]), linestyle=line_style, color= color_def)
    plt.fill_between(L_delta_g_bb, L_std_low, L_std_high, color= color_def, alpha=0.1)
    plt.title("%s: Effect $\Delta G_{pol}$ on error probability of the polymer string" %(model_title[4:]))
    plt.xlabel("Backbone strength $\Delta G_{pol}$")
    plt.ylabel("Error probability $\epsilon$")
    plt.legend(bbox_to_anchor=(1.04, 1))
    plt.subplots_adjust(right=0.7)
    return

def graph_error_g_bb_multiple_models_per_model3(L_delta_g_bb, L_epsilon, model, delta_g_tt, model_title, L_std_high, L_std_low, color_def, line_style):
    plt.plot(L_delta_g_bb, L_epsilon, label = "$\Delta G_{tt}$ = %s, %s" %(delta_g_tt, model_title), linestyle=line_style, color= color_def)
    plt.fill_between(L_delta_g_bb, L_std_low, L_std_high, color= color_def, alpha=0.1)
    plt.title("Comparison models: Effect $\Delta G_{pol}$ on error probability of the polymer string")
    plt.xlabel("Backbone strength $\Delta G_{pol}$")
    plt.ylim(0,0.6)
    plt.ylabel("Error probability $\epsilon$")
    plt.legend(bbox_to_anchor=(1.3, 1))
    plt.subplots_adjust(right=0.7)
    return


def make_model_names(x_def, M_model, models, models_lose, models_titles):
    L_models = []
    L_models_titles = []
    for x in x_def:
        if x == "x1":
            x_title = "1"
        else:
            x_title = "$e^{\Delta G_{TT} +5}$"
        for M in (M_model):
            if M =="M1":
                M_title = "Model B"
            else: 
                M_title = "Model F"
            for i,model in enumerate(models):
                model = model + M + "_" + x
                L_models.append(model)
                L_models_titles.append(models_titles[i] +" " + M_title+ " x=" + x_title)
    for mod in models_lose:
        L_models.append(mod)
        L_models_titles.append(models_titles[-1])
    return L_models, L_models_titles


L_delta_g_pol_label = ["-ln2", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
L_delta_g_pol_number = [-np.log(2), 0,1,2,3,4,5,6,7,8,9,10]
L_delta_g_tt = [0,2,4,6,8,10]

models = ["muriel_big_abs_", "no_abs_matrix_", "muriel_"]
models_title = ["FGA", "FGC", "Small absorption matrix", "CG"]
all_steps = True
DIR = "230221output" #"230201output"
x_def = ["xexp5", "x1"] #should be xdelta_g_tt_pow2
M_model = ["M1", "M2"]
# number_steps = 0
# MX_models = ["x=$e^{\Delta G_{TT} +5}$, Model B","x=$e^{\Delta G_{TT} +5}$, Model F", "x=1, Model B", "x=1, Model F" ]

L_models, L_models_titles = make_model_names(x_def, M_model, models, ["Jenny_full_steps"], models_title)
date = date_now(datetime.datetime.now())
date_options = ["2023210", "2023211", "2023213", "2023221", "2023223", "2023224"]#["202323", "202324", "202306"]


#define colors
colormap = plt.cm.plasma #nipy_spectral, Set1,Paired   
color_list = [colormap(i) for i in np.linspace(0, 1,100)]

i = 0
# L_models = [L_models[i], L_models[i+3], L_models[i+6], L_models[i+9]] #per model all different M and x
# L_models_titles = [L_models_titles[i], L_models_titles[i+3], L_models_titles[i+6], L_models_titles[i+9]]
# L_models = [L_models[i], L_models[i+3]] #per model all different M and x = exp
# L_models_titles = [L_models_titles[i], L_models_titles[i+3]] #per model all different M and x = exp
# L_models = [L_models[i+6], L_models[i+9]] #per model all different M and x = 1
# L_models_titles = [L_models_titles[i+6], L_models_titles[i+9]]#per model all different M and x = 1
# L_models = [L_models[i], L_models[i+1]] #FGA FGC for one x and for one M 
# L_models_titles = [L_models_titles[i], L_models_titles[i+1]] #(i=0,3,6 of 9)

L_models = [L_models[i], L_models[i+3], L_models[12]] #FGA for one x and for two M and CG
L_models_titles = [L_models_titles[i], L_models_titles[i+3],  L_models_titles[12]] ##FGA for one x and for two M and CG (i=0,3,6 of 9)
L_delta_g_tt = [4]



# line_styles = ['dashed', 'solid']

# L_models = L_models[i:i+3] #for same M and X diffent models
line_styles = ['solid', 'dashdot', 'dashed', 'dotted']
print(L_models)
L_gtt_theorie = [0.5, 0.19251027051493744, 0.034723337448322934, 0.004920911195322951, 0.0006702507235977487, 9.078749428737782e-05]

for k, model in enumerate(L_models):
    for j, delta_g_tt in enumerate(L_delta_g_tt):
        L_epsilon_mean = []
        L_epsilon_std_high = []
        L_epsilon_std_low = []
        L_delta_g_bb_graph = []
        for i,delta_g_pol in enumerate(L_delta_g_pol_number):          
            for date in date_options:
                # print('/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/%s/%smodel_%s_delta_g_pol_%s_delta_g_tt_%s_allsteps.csv'%(DIR, date, model, L_delta_g_pol_label[i], delta_g_tt))                                                
                file_exists = os.path.exists('/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/%s/%smodel_%s_delta_g_pol_%s_delta_g_tt_%s_allsteps.csv'%(DIR, date, model, L_delta_g_pol_label[i], delta_g_tt))
                
                if file_exists == True:
                    data = pd.read_csv('/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/%s/%smodel_%s_delta_g_pol_%s_delta_g_tt_%s_allsteps.csv'%(DIR, date, model, L_delta_g_pol_label[i], delta_g_tt))
                    epsilon_mean = data['error probability'].mean()
                    epsilon_std = data['error probability'].std()
                    L_epsilon_mean.append(epsilon_mean)
                    L_epsilon_std_high.append(epsilon_mean + epsilon_std)
                    L_epsilon_std_low.append(epsilon_mean - epsilon_std)
                    L_delta_g_bb_graph.append(delta_g_pol)


                else:
                    print("file does not exist")
                    # file_exists = os.path.exists('/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/%s/%smodel_%s_delta_g_pol_%s_delta_g_tt_%s_allsteps.csv'%(DIR, '202324', model, L_delta_g_pol_label[i], delta_g_tt))
                    # if file_exists == True:
                    #     data = pd.read_csv('/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/%s/%smodel_%s_delta_g_pol_%s_delta_g_tt_%s_allsteps.csv'%(DIR, date, model, L_delta_g_pol_label[i], delta_g_tt))
                    #     epsilon_mean = data['error probability'].mean()
                    #     epsilon_std = data['error probability'].std()
                    #     L_epsilon_mean.append(epsilon_mean)
                    #     L_epsilon_std_high.append(epsilon_mean + epsilon_std)
                    #     L_epsilon_std_low.append(epsilon_mean - epsilon_std)
                    #     L_delta_g_bb_graph.append(delta_g_pol)
                    # else:
                    #     print("file does not exist")


        #define error as all mismatches and plot the error as function of delta G_bb
        # graph_error_g_bb(L_delta_g_bb_graph, L_epsilon_mean, model, delta_g_tt, L_models_titles[k], L_epsilon_std_high, L_epsilon_std_low, L_gtt_theorie[j], color_list[j*10])
        # graph_error_g_bb_multiple_models_per_model(L_delta_g_bb_graph, L_epsilon_mean, model, delta_g_tt, L_models_titles[k], L_epsilon_std_high, L_epsilon_std_low, color_list[j*10], line_styles[k])
        # graph_error_g_bb_multiple_models_per_model2(L_delta_g_bb_graph, L_epsilon_mean, model, delta_g_tt, L_models_titles[k], L_epsilon_std_high, L_epsilon_std_low, color_list[j*10], line_styles[k])
        graph_error_g_bb_multiple_models_per_model3(L_delta_g_bb_graph, L_epsilon_mean, model, delta_g_tt, L_models_titles[k], L_epsilon_std_high, L_epsilon_std_low, color_list[20], line_styles[k])

        # print(model, delta_g_tt, L_epsilon_mean)
plt.show()


