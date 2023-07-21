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
from pylab import *


def date_now(now):
    year = '{:02d}'.format(now.year)
    month = f'{now.month}' 
    day = f'{now.day}' 
    string_date = year + month + day
    print(string_date)
    return string_date


def graph_error_g_bb(L_delta_g_bb, L_epsilon, model_title, bases, L_std_high, L_std_low, color_def, repulsion):

    plt.plot(L_delta_g_bb, L_epsilon, label = "bases = %s" %(bases), color= color_def)
    plt.fill_between(L_delta_g_bb, L_std_low, L_std_high, color= color_def, alpha=0.1)
    if repulsion =="-1":
        plt.title("%s: Effect $\Delta G_{pol}$ on error probability of the copy polymer, not allowing pur-pur or py-py binding"%(model_title))
    else:
        plt.title("%s: Effect $\Delta G_{pol}$ on error probability of the copy polymer, with repulsion %s"%(model_title,repulsion))
    plt.ylim(0,0.6)
    plt.xlabel("Backbone strength $\Delta G_{pol}$")
    plt.ylabel("Error probability $\epsilon$")
    plt.legend(bbox_to_anchor=(1.04, 1))
    plt.subplots_adjust(right=0.7)
    return


def graph_error_g_bb_3tables(L_delta_g_bb, L_epsilon, model_title, bases, L_std_high, L_std_low, color_def, repulsion, matrix_table1, colnames,iteration,graph,table1, table2, table3,error_limit,fontsizee, lineswidth,base_legend):
    # matrix_table1 = matrix_table1.round(decimals=2)
    if iteration ==0:
        table1.axis('tight')
        table1.axis('off')
        the_table1 = table1.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        the_table1.auto_set_font_size(False)
        the_table1.set_fontsize(fontsizee)
        table2.plot()
        table3.plot()
    if iteration ==1:
        table2.axis('tight')
        table2.axis('off')
        the_table2 = table2.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        the_table2.auto_set_font_size(False)
        the_table2.set_fontsize(fontsizee)
        table1.plot()
        table3.plot()
    if iteration ==2:
        table3.axis('tight')
        table3.axis('off')
        the_table3 =table3.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        the_table3.auto_set_font_size(False)
        the_table3.set_fontsize(fontsizee)
        table1.plot()
        table2.plot()

    if bases == "aATLambdaSigma" or bases == "aATGC" or bases == "aATXK":
        if repulsion != "-1":
            the_graph = graph.scatter(12, error_limit, color= color_def)
            labels = ["-2","0","2","4","6","8","10", "$\infty$"]
            graph.set_xticks([-2,0,2,4,6,8,10,12])
            graph.set_xticklabels(labels,fontsize =fontsizee)

    bases_str  = bases
    if bases == "aATLambdaSigma":
        bases_str = "aAT$\Lambda\Sigma$" 
    if bases == "ATLambdaSigma":
        bases_str = "AT$\Lambda\Sigma$"\
    
    the_graph = graph.plot(L_delta_g_bb, L_epsilon, label = "%s" %(base_legend), color= color_def,linewidth = lineswidth)
    graph.fill_between(L_delta_g_bb, L_std_low, L_std_high, color= color_def, alpha=0.1)
    graph.set_ylim([0,0.8])
    graph.set_xlabel("Backbone strength $\Delta G_{bb}$", fontsize =fontsizee)
    graph.set_ylabel("Error probability $\epsilon$",fontsize =fontsizee)
    graph.tick_params(axis='y', labelsize=fontsizee, width=lineswidth*0.5) 
    graph.tick_params(axis='x', labelsize=fontsizee, width=lineswidth*0.5)
    graph.spines['left'].set_linewidth(lineswidth*0.5)
    graph.spines['right'].set_linewidth(lineswidth*0.5)
    graph.spines['bottom'].set_linewidth(lineswidth*0.5)
    graph.spines['top'].set_linewidth(lineswidth*0.5)
    # graph.set_xticklabels(L_delta_g_bb,fontsize=16)
    box = graph.get_position()
    graph.set_position([box.x0, box.y0+box.height*0.1, box.width, box.height*0.9])
    graph.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3, fontsize =fontsizee)

    if repulsion =="-1":
        graph.set_title("Effect $\Delta G_{bb}$ on error probability of the copy polymer, not allowing same group matches", fontsize =fontsizee,y=1.0, pad=14)
    else:
        graph.set_title("Effect $\Delta G_{bb}$ on error probability of the copy polymer, with free energy cost %s"%(repulsion), fontsize =fontsizee, y=1.0, pad=14)
    
    return








def change_base_string_to_int(base):
    base_int = 0
    if base == "aA":
        base_int =1
    elif base =="A":
        base_int = 15
    elif base == "T":
        base_int = 2
    elif (base == "G"):
        base_int = 3
    elif (base == "C"):
        base_int = 4
    elif (base == "iG"):
        base_int = 5
    elif (base == "iC"):
        base_int = 6
    elif (base == "S"):
        base_int = 7
    elif (base == "L"):
        base_int = 8
    elif (base == "Lambda"):
        base_int = 9
    elif (base == "Sigma"):
        base_int = 10
    elif (base == "delta"):
        base_int = 11
    elif (base == "beta"):
        base_int = 12
    elif (base == "X"):
        base_int = 13
    elif (base == "K"):
        base_int = 14
    return base_int



def base_OH_groups(base):
    if base == 1:
        #base = 1 means, base is aA.
        OH_groups = [1, 0, 1]
    elif base == 15:
        OH_groups = [1,0,5]
    elif base == 2:
        #base = 2 means, base is T.
        OH_groups = [0, 1, 0]
    elif base == 3:
        #base = 3 means, base is G.
        OH_groups = [0, 1, 1]
    elif base ==  4:
        #base = 4 means, base is C.
        OH_groups = [1, 0, 0]
    elif base ==  5:
        #base = 5 means, base is iG.
        OH_groups = [1, 1, 0]
    elif base == 6:
        #base = 6 means, base is iC.
        OH_groups = [0, 0, 1]
    elif base == 7:
        #base = 7 means, base is S.
        OH_groups = [0, 0, 0]
    elif base == 8:
        #base = 8 means, base is L.
        OH_groups = [1, 1, 1]
    
    elif base == 9:
        #base = 9 means, base is Lambda.
        OH_groups = [1, 1, 1]
    
    elif base == 10:
        #base = 10 means, base is Sigma.
        OH_groups = [0, 0, 0]
    elif base == 11:
        #base = 11 means, base is delta.
        OH_groups = [1, 0, 0]
    
    elif base == 12:
        #base = 12 means, base is beta.
        OH_groups = [0, 1, 1]
    
    elif base == 13:
        #base = 12 means, base is X.
        OH_groups = [0, 1, 0]
    elif base == 14:
        #base = 14 means, base is L.
        OH_groups = [1, 0, 1]
    return OH_groups


def OH_matches_strength(base_i, base_j, rep):
    match = 0
    base = False
    strength = 0 
    for i in range(0,len(base_i)):
        if base_i[i]==5 or base_j[i]==5:
            base = True
    if base == False:
        for i in range(0,len(base_i)):
            if base_i[i] != base_j[i]:
                match += 1
        strength = -1*(2 * match - 3) +rep
    elif base == True:
        for i in range(0,len(base_i)):
            if base_i[i] ==5 or base_j[i] ==5:
                strength = strength
            else:
                if base_i[i] != base_j[i]:
                    match += 1
        strength = -1 * (2*match -2) +rep
    return strength

def make_E_OH_matrix(bases, rep):
    
    OH_matrix=np.zeros((len(bases), len(bases)))
    matrix = np.array([["        " for c in range(len(bases))] for r in range(len(bases))])
    strength = 0

    for i, basei in enumerate (bases):
        for j, basej in enumerate (bases):
            if basei =="L":
                basei = "Lambda"
            if basei =="S":
                basei = "Sigma"
            if basej =="L":
                basej = "Lambda"
            if basej =="S":
                basej = "Sigma"
            base_i_int = change_base_string_to_int(basei)
            base_i = base_OH_groups(base_i_int)
            base_j_int = change_base_string_to_int(basej)
            base_j = base_OH_groups(base_j_int)
            base_i_dev = base_i_int%2
            base_j_dev = base_j_int%2
            tot_dev = base_i_dev +base_j_dev 
            if tot_dev ==1:
                strength = OH_matches_strength(base_i,base_j,0)
                strength_str = str(round(strength,2))
            else: 
                strength = OH_matches_strength(base_i,base_j,rep)
                strength_str = str(round(strength,2))
                if rep ==-1:
                    strength = 100
                    strength_str = "$\infty$"
            # strength = OH_matches_strength(base_i,base_j,0)
            OH_matrix[i,j] = strength
            matrix[i,j] = strength_str
            print(matrix)
    
    # OH_matrix_str = [[str(ele) for ele in sub] for sub in OH_matrix] 
    # strengths = 'aa'
    # if rep == -1:
    #     for i, basei in enumerate (bases):
    #         for j, basej in enumerate (bases):
    #             if OH_matrix[i,j] ==100:
    #                 strengths = 'aa'
    #             else:
    #                 strengths = str(OH_matrix[i,j])
    #             OH_matrix_str[i,j] = strengths
    
        
    return OH_matrix, matrix


def matrix_purpy(bases):
    # pur1pur1, pur1pur2, pur2pur2, py1py1, py1py2, py2py2, py1pur1, py1pur2 ,py2pur1 ,py2pur2 = probs
    matrix = []
    for i, basei in enumerate(bases):
        row = []
        for j, basej in enumerate(bases):
            row.append(basei +basej)
        matrix.append(row)
    return matrix

def make_model_names(bases_combies, model, models_title):
    L_models= []
    L_model_titles = []
    for bases in bases_combies:
        L_models.append(model + bases)
        string = models_title + " " + bases
        L_models_titles.append(string)

    return L_models, L_models_titles

def graph_error_high_g_bb_all_bases(L_error_high_bb, L_repulsion, bases_str, colors_bases, fontsizee, lineswidth, base_legend):
    f, (ax, ax2) = plt.subplots(1, 2, sharey=True, facecolor='w', gridspec_kw={'width_ratios': (5,5)})
    L_errors_base1 = []
    L_errors_base2 = []
    L_errors_base3 = []
    L_errors_sem_base1 = []
    L_errors_sem_base2 = []
    L_errors_sem_base3 = []
    for i in range (0,len(L_repulsion)):
        if L_repulsion[i] ==-1:
            L_repulsion[i] = 10
        if L_repulsion[i] == 20:
            L_repulsion[i] = 8
        L_errors_base1.append(L_error_high_bb[3*i])
        L_errors_base2.append(L_error_high_bb[3*i+1])
        L_errors_base3.append(L_error_high_bb[3*i+2])
    if bases_str[0] =="aATGCLambdaSigma":
        bases_str[0] ="aATGC$\Lambda\Sigma$"
    if bases_str[0] =="ATLambdaSigma":
        bases_str[0] ="AT$\Lambda\Sigma$"

    graph = ax.errorbar(L_repulsion, L_errors_base1,label = "%s" %(base_legend[0]),fmt='-o', color =colors_bases[0], linewidth = lineswidth, markersize=lineswidth*3)
    graph = ax.errorbar(L_repulsion, L_errors_base2, label = "%s" %(base_legend[1]),fmt='-o',color =colors_bases[1], linewidth = lineswidth, markersize=lineswidth*3)
    graph = ax.errorbar(L_repulsion, L_errors_base3, label = "%s" %(base_legend[2]),fmt='-o',color =colors_bases[2], linewidth = lineswidth, markersize=lineswidth*3)
    graph = ax2.errorbar(L_repulsion, L_errors_base1,  label = "%s" %(base_legend[0]),fmt='-o',color =colors_bases[0], linewidth = lineswidth, markersize=lineswidth*3)
    graph = ax2.errorbar(L_repulsion, L_errors_base2, label = "%s" %(base_legend[1]),fmt='-o',color =colors_bases[1], linewidth = lineswidth, markersize=lineswidth*3)
    graph = ax2.errorbar(L_repulsion, L_errors_base3,label = "%s" %(base_legend[2]),fmt='-o',color =colors_bases[2], linewidth = lineswidth, markersize=lineswidth*3)

    ax.set_ylim([0,0.8])
    ax.set_xlim(-0.5, 6.5)
    ax2.set_xlim(7.5, 10.5)
    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.spines['left'].set_linewidth(lineswidth*0.5)
    ax2.spines['right'].set_linewidth(lineswidth*0.5)
    ax.spines['bottom'].set_linewidth(lineswidth*0.5)
    ax2.spines['bottom'].set_linewidth(lineswidth*0.5)
    ax.spines['top'].set_linewidth(lineswidth*0.5)
    ax2.spines['top'].set_linewidth(lineswidth*0.5)
    # ax2.spines['top'].set_visible(False)
    # ax.spines['top'].set_visible(False)
    ax.yaxis.tick_left()
    ax2.yaxis.tick_right()
    ax2.tick_params(labelright='off')
    ax2.set_xticks([8,10])
    labels = ["0","1","2","3","4","5","6","20", "100"]
    ax.set_xticks([0,1,2,3,4,5,6])
    ax.set_xticklabels(labels[0:7],fontsize =fontsizee)
    ax2.set_xticklabels(labels[7:], fontsize =fontsizee)
    ax.tick_params(axis='y', labelsize=fontsizee, width=lineswidth*0.5)
    ax2.tick_params(axis='y', labelsize=fontsizee,width=lineswidth*0.5)
    ax.tick_params(axis='x', labelsize=fontsizee, width=lineswidth*0.5)
    ax2.tick_params(axis='x', labelsize=fontsizee,width=lineswidth*0.5)
    ax.set_xlabel("Same group matches cost", fontsize =fontsizee)
    ax.xaxis.set_label_coords(1.0, -.1)
    ax.set_ylabel("Error probability $\epsilon$", fontsize =fontsizee)
    # ax2.set_ylabel("Error probability $\epsilon$", fontsize =fontsizee)
    ax2.yaxis.set_label_coords(1.15, 0.5)
    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1-d, 1+d), (-d, +d), **kwargs)
    ax.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
    ax2.plot((-d, +d), (-d, +d), **kwargs)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0+box.height*0.1, box.width, box.height*0.9])
    ax2.set_position([box.x0+box.width +2*d, box.y0+box.height*0.1, box.width, box.height*0.9])
    ax.legend(loc='lower center', bbox_to_anchor=(1.0, -0.28), ncol=3, fontsize =fontsizee)
    return

plot1 = True



if plot1 ==True:
    L_delta_g_pol_label = ["-ln4","-ln2", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    L_delta_g_pol_number = [-np.log(4),-np.log(2), 0,1,2,3,4,5,6,7,8,9,10]
    L_errors = [0.682475291059189, 0.682475291059189, 0.682475291059189,0.63217714067757, 0.6202528585865767, 0.6429263771115179,0.5814028034205805, 0.5319121191342184, 0.5814028034205805,0.5447099750872387, 0.44980382969388233, 0.49615699865355445,0.5255168981206141, 0.3994690101758146, 0.3994690101758146,0.5172851808956745, 0.3760832287354581, 0.32271793830029694,0.5140675227826148, 0.36662525281293135, 0.28008919549588873,0.5121444504545963, 0.36088452129225734, 0.24931594522843759,0.5121444488396317, 0.36088451644323927, 0.24931591736204894]

    all_steps = True
    # DIR = "230414output"
    L_DIR = ["230529output", "230529_equi_output"]#"230517_equi_output","230517output"]
    model = "muriel_big_abs_BF_x1"
    model_title = "BF x=1"
    bases_good=[["aA", "T", "G","C","$\Lambda$", "$\Sigma$"],["aA","T","G","C","iG","iC"],["aA","T","G","C","X","K"]]
    bases_pur1py1pur2py2=[["aA", "T","G","C", "L", "S"],["aA","T","G","C","iG","iC"],["aA","T","G","C","X","K"]]#[["A", "T", "L", "S"],["A","T","G","C"],["A","T","X","K"]]#
    base_combinations = ["aATGCLambdaSigma","aATGCiGiC", "aATGCXK"] #["ATLambdaSigma","ATGC", "ATXK"]# "GCXK","aATXK","aATGC", "deltabetaXK",
    Base_legend = ["Opposite group", "Equal", "Same group"]
    # Base_legend = base_combinations
    ratio ='1111'
    # number_steps = 0
    # MX_models = ["x=$e^{\Delta G_{TT} +5}$, Model B","x=$e^{\Delta G_{TT} +5}$, Model F", "x=1, Model B", "x=1, Model F" ]

    # L_models, L_models_titles = make_model_names(base_combinations, model, model_title)exit()\e
    L_models = [model]
    L_models_titles = [model_title]
    date = date_now(datetime.datetime.now())
    date_options = ["2023530", "2023531", "2023621"]#"2023517","2023518", "2023519", "2023524"]
    L_repulsion = ["0","1","2","3","4", "5", "6","20", "-1"]
    number_tables=3

    #define colors
    colormap = plt.cm.plasma #plt.cm.plasma #nipy_spectral, Set1,Paired   
    color_list = [colormap(i) for i in np.linspace(0, 1,100)]

    print(L_models)
    OH_matrix = np.zeros((4, 4))
    fontsizee = 23
    line_widths = 3

    # line_styles = ['dashed', 'solid']

    # L_models = L_models[i:i+3] #for same M and X diffent models
    line_styles = ['solid', 'dashdot', 'dashed', 'dotted', (5, (10, 3)), (0, (3, 5, 1, 5))]
    print(L_models)


    for k, model in enumerate(L_models):
        L_error_high_bb = []
        L_sem_error_high_g_bb = []
        L_rep_numbers= []
        for r, repulsion in enumerate(L_repulsion):
            L_colors = []
            for j, base in enumerate(base_combinations):
                L_epsilon_mean = []
                L_epsilon_std_high = []
                L_epsilon_std_low = []
                L_delta_g_bb_graph = []
                L_sem = []
                for i,delta_g_pol in enumerate(L_delta_g_pol_number): 
                    rep = float(repulsion)
                    label = L_delta_g_pol_label[i]
                    # if rep > -1:
                        # if L_delta_g_pol_label[i] == '-ln4':
                            # label = "-lnrep" 
                            # delta_g_pol = -np.log(2 + np.exp(-rep))  
                    for DIR in L_DIR:
                        for date in date_options:

                            L_delta_g_bb_graph.append(delta_g_pol)
                            repulsion_int = float(repulsion)
                            print(bases_pur1py1pur2py2[j], "repulsion:", repulsion_int, "G_pol:",delta_g_pol)
                            # print(matrix_purpy(bases_pur1py1pur2py2[j]))
                            OH_matrix, matrix = make_E_OH_matrix(bases_pur1py1pur2py2[j], repulsion_int)
                            print(OH_matrix)
                            
                if number_tables ==3:
                    L_colors.append(color_list[j*30])
                                    
            plt.show()
            L_rep_numbers.append(rep)

        graph_error_high_g_bb_all_bases(L_errors,L_rep_numbers, base_combinations, L_colors, fontsizee, line_widths, Base_legend)
        print("base:",Base_legend)
        print("error", L_errors[r+j])
        plt.show()
