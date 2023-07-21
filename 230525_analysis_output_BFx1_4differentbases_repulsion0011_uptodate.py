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

def graph_equilibration_check(length_polymer, error_prob,delta_g_bb, delta_g_tt, model):
    plt.plot(length_polymer, error_prob, label = "$\Delta G_{bb}$ = %s, $\Delta G_{tt}$ = %s"%(delta_g_bb, delta_g_tt))
    plt.title("Energetic: model: %s, Check when error is equilibrated as function of the length of the mRNA grown for different $\Delta G_{bb}$ and $\Delta G_{tt}$" %model)
    plt.xlabel("length mRNA")
    plt.ylabel("$\epsilon$")
    plt.legend(bbox_to_anchor=(1.04, 1))
    plt.subplots_adjust(right=0.7)
    return

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

def graph_error_g_bb_5tables(L_delta_g_bb, L_epsilon, model_title, bases, L_std_high, L_std_low, color_def, repulsion, matrix_table1, colnames,iteration,graph,table1, table2, table3,table4, table5, error_limit,fontsizee, lineswidth,base_legend):
    # matrix_table1 = matrix_table1.round(decimals=2)
    if iteration ==0:
        table1.axis('tight')
        table1.axis('off')
        the_table1 = table1.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        the_table1.auto_set_font_size(False)
        the_table1.set_fontsize(fontsizee)
        table2.plot()
        table3.plot()
        table4.plot()
        table5.plot()
    if iteration ==1:
        table2.axis('tight')
        table2.axis('off')
        the_table2 = table2.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        the_table2.auto_set_font_size(False)
        the_table2.set_fontsize(fontsizee)
        table1.plot()
        table3.plot()
        table4.plot()
        table5.plot()

    if iteration ==2:
        table3.axis('tight')
        table3.axis('off')
        the_table3 =table3.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        the_table3.auto_set_font_size(False)
        the_table3.set_fontsize(fontsizee)
        table1.plot()
        table2.plot()
        table4.plot()
        table5.plot()

    if iteration ==3:
        table4.axis('tight')
        table4.axis('off')
        the_table4 =table4.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        the_table4.auto_set_font_size(False)
        the_table4.set_fontsize(fontsizee)
        table1.plot()
        table2.plot()
        table3.plot()
        table5.plot()

    if iteration ==4:
        table5.axis('tight')
        table5.axis('off')
        the_table5 =table5.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        the_table5.auto_set_font_size(False)
        the_table5.set_fontsize(fontsizee)
        table1.plot()
        table2.plot()
        table4.plot()
        table3.plot()


    # if bases == "ATLambdaSigma" or bases == "ATGC" or bases == "ATXK" or bases == "ATiGiC" or bases == "ATalphabeta":
    #     if repulsion != "-1":
    #         the_graph = graph.scatter(12, error_limit, color= color_def)
    #         labels = ["-2","0","2","4","6","8","10", "$\infty$"]
    #         graph.set_xticks([-2,0,2,4,6,8,10,12])
    #         graph.set_xticklabels(labels,fontsize =fontsizee)

    bases_str  = bases
    if bases == "aATLambdaSigma":
        bases_str = "aAT$\Lambda\Sigma$" 
    if bases == "ATLambdaSigma":
        bases_str = "AT$\Lambda\Sigma$"
    if bases == "ATdeltabeta":
        bases_str = "AT$\delta\\beta$" 

    
    the_graph = graph.plot(L_delta_g_bb, L_epsilon, label = "%s" %(bases_str), color= color_def,linewidth = lineswidth)
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
    graph.legend(loc='lower center', bbox_to_anchor=(0.4, -0.5), ncol=3, fontsize =fontsizee)

    if repulsion =="-1":
        graph.set_title("Effect $\Delta G_{bb}$ on error probability of the copy polymer, not allowing same group matches", fontsize =fontsizee,y=1.05, x =0.8, pad=14)
    else:
        graph.set_title("Effect $\Delta G_{bb}$ on error probability of the copy polymer, with free energy cost %s"%(repulsion), fontsize =fontsizee, y=1.05, x = 0.8, pad=14)
    
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



def graph_error_g_bb_6tables(L_delta_g_bb, L_epsilon, model_title, bases, L_std_high, L_std_low, color_def, repulsion, matrix_table1, colnames,iteration,graph,table1, table2, table3, table4,table5,table6,OH_matrix):
    matrix_table1 = matrix_table1.round(decimals=2)
    if iteration ==0:
        table4.axis('tight')
        table4.axis('off')
        the_table4 = table4.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        the_table4.auto_set_font_size(False)
        the_table4.set_fontsize(16)
        table1.axis('tight')
        table1.axis('off')  
        the_table1 = table1.table(cellText=OH_matrix,colLabels=colnames,rowLabels = colnames, loc='center')
        the_table1.auto_set_font_size(False)
        the_table1.set_fontsize(16)
        table2.plot()
        table3.plot()
        table5.plot()
        table6.plot()
    if iteration ==1:
        table5.axis('tight')
        table5.axis('off')
        the_table5 = table5.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        the_table5.auto_set_font_size(False)
        the_table5.set_fontsize(16)
        table2.axis('tight')
        table2.axis('off')  
        the_table2 = table2.table(cellText=OH_matrix,colLabels=colnames,rowLabels = colnames, loc='center')
        the_table2.auto_set_font_size(False)
        the_table2.set_fontsize(16)
        table1.plot()
        table3.plot()
        table4.plot()
        table6.plot()
    if iteration ==2:
        table6.axis('tight')
        table6.axis('off')
        the_table6 = table6.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        the_table6.auto_set_font_size(False)
        the_table6.set_fontsize(16)
        table3.axis('tight')
        table3.axis('off')  
        the_table3 = table3.table(cellText=OH_matrix,colLabels=colnames,rowLabels = colnames, loc='center')
        the_table3.auto_set_font_size(False)
        the_table3.set_fontsize(16)
        table1.plot()
        table2.plot()
        table4.plot()
        table5.plot()
        

    bases_str  = bases
    if bases == "aATLambdaSigma":
        bases_str = "aAT$\Lambda\Sigma$"
    if bases == "ATLambdaSigma":
        bases_str = "AT$\Lambda\Sigma$"
    graph.plot(L_delta_g_bb, L_epsilon, label = "%s" %(bases_str), color= color_def)
    graph.fill_between(L_delta_g_bb, L_std_low, L_std_high, color= color_def, alpha=0.1)
    graph.set_ylim([0,0.8])
    graph.set_xlabel("Backbone strength $\Delta G_{bb}$", fontsize =16)
    graph.set_ylabel("Error probability $\epsilon$", fontsize =16)
    box = graph.get_position()
    graph.set_position([box.x0, box.y0+box.height*0.1, box.width, box.height*0.9])
    graph.legend(loc='lower center', bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize =16)

    if repulsion =="-1":
        graph.set_title("Effect $\Delta G_{bb}$ on error probability of the copy polymer, not allowing same group matches", fontsize =16)
    else:
        graph.set_title("Effect $\Delta G_{bb}$ on error probability of the copy polymer, with free energy cost %s"%(repulsion), fontsize =16)
    
    return

def graph_cummulative_3tables(L_delta_g_bb, L_pur_pur_incorrect, L_py_py_incorrect, L_pur_py_incorrect, L_pur_py_correct, model_title, bases, repulsion, matrix_table1, colnames,graph,table2):
    # matrix_table1 = matrix_table1.round(decimals=2)
    table2.axis('tight')
    table2.axis('off')
    the_table2 = table2.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
    the_table2.auto_set_font_size(False)
    the_table2.set_fontsize(16)
    # table1.plot()
    # table3.plot()

    bases_str  = bases
    if bases == "aATLambdaSigma":
        bases_str = "aAT$\Lambda\Sigma$" 
    if bases == "ATLambdaSigma":
        bases_str = "AT$\Lambda\Sigma$"
    the_graph = graph.plot(L_delta_g_bb, L_pur_pur_incorrect, label = "pur-pur incorrect", color = "r")
    the_graph = graph.plot(L_delta_g_bb, L_py_py_incorrect, label = "py-py incorrect", color = "orange")
    the_graph = graph.plot(L_delta_g_bb, L_pur_py_incorrect, label = "pur-py incorrect", color = "darkorchid")
    the_graph = graph.plot(L_delta_g_bb, L_pur_py_correct, label = "pur-py correct", color = "g")
    graph.set_ylim([0,1.0])
    graph.set_xlabel("Backbone strength $\Delta G_{pol}$", fontsize =16)
    graph.set_ylabel("cumulative error probability $\epsilon$",fontsize =16)
    graph.tick_params(axis='y', labelsize=16)
    graph.tick_params(axis='x', labelsize=16)
    # graph.set_xticklabels(L_delta_g_bb,fontsize=16)
    box = graph.get_position()
    graph.set_position([box.x0, box.y0+box.height*0.1, box.width, box.height*0.9])
    graph.legend(loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=3, fontsize =16)

    if repulsion =="-1":
        graph.set_title("%s, %s: Effect $\Delta G_{pol}$ on error probability of the copy polymer, not allowing pur-pur or py-py binding"%(model_title,bases_str), fontsize =16)
    else:
        graph.set_title("%s, %s: Effect $\Delta G_{pol}$ on error probability of the copy polymer, with free energy cost %s"%(model_title,bases_str, repulsion), fontsize =16)
    
    return

def graph_error_high_g_bb_all_bases(L_error_high_bb, L_SEM_error_high_g_bb, L_repulsion, bases_str, colors_bases, fontsizee, lineswidth, base_legend):
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
        L_errors_sem_base1.append(L_SEM_error_high_g_bb[3*i])
        L_errors_sem_base2.append(L_SEM_error_high_g_bb[3*i+1])
        L_errors_sem_base3.append(L_SEM_error_high_g_bb[3*i+2])
    if bases_str[0] =="aATLambdaSigma":
        bases_str[0] ="aAT$\Lambda\Sigma$"
    if bases_str[0] =="ATLambdaSigma":
        bases_str[0] ="AT$\Lambda\Sigma$"

    graph = ax.errorbar(L_repulsion, L_errors_base1,yerr=L_errors_sem_base1,label = "%s" %(base_legend[0]),fmt='-o', color =colors_bases[0], linewidth = lineswidth, markersize=lineswidth*3)
    graph = ax.errorbar(L_repulsion, L_errors_base2, yerr=L_errors_sem_base2,label = "%s" %(base_legend[1]),fmt='-o',color =colors_bases[1], linewidth = lineswidth, markersize=lineswidth*3)
    graph = ax.errorbar(L_repulsion, L_errors_base3, yerr=L_errors_sem_base3, label = "%s" %(base_legend[2]),fmt='-o',color =colors_bases[2], linewidth = lineswidth, markersize=lineswidth*3)
    graph = ax2.errorbar(L_repulsion, L_errors_base1, yerr=L_errors_sem_base1, label = "%s" %(base_legend[0]),fmt='-o',color =colors_bases[0], linewidth = lineswidth, markersize=lineswidth*3)
    graph = ax2.errorbar(L_repulsion, L_errors_base2,yerr=L_errors_sem_base2, label = "%s" %(base_legend[1]),fmt='-o',color =colors_bases[1], linewidth = lineswidth, markersize=lineswidth*3)
    graph = ax2.errorbar(L_repulsion, L_errors_base3, yerr=L_errors_sem_base3,label = "%s" %(base_legend[2]),fmt='-o',color =colors_bases[2], linewidth = lineswidth, markersize=lineswidth*3)

    ax.set_ylim([0,0.6])
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
    labels = ["0","1","2","3","4","5","6","20", "$\infty$"]
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

def graph_error_high_g_bb_all_bases_5(L_error_high_bb, L_SEM_error_high_g_bb, L_repulsion, bases_str, colors_bases, fontsizee, lineswidth, base_legend):
    f, (ax, ax2) = plt.subplots(1, 2, sharey=True, facecolor='w', gridspec_kw={'width_ratios': (5,5)})
    L_bases_aA = ["Opposite group", "Equal", "Same group"]#['aAT$\\Lambda\\Sigma$', 'aATGC', 'aATXK']
    L_error_high_bb_aA = [0.5218136237623762, 0.519714495049505, 0.5181401485148515, 0.4749326435643565, 0.45142219801980193, 0.49742857425742565, 0.44886822772277224, 0.36549699009900993, 0.4503376336633664, 0.4443550693069307, 0.2931462475247524, 0.36798096039603967, 0.43221495049504954, 0.24682503960396043, 0.24857436633663366, 0.43868737623762377, 0.2259735247524752, 0.1502291287128713, 0.4294860792079208, 0.22464403960396043, 0.09033335742574258, 0.4308855049504951, 0.21943115841584157, 0.04639121584158416, 0.4325648613861387, 0.22020082178217823, 0.04842039702970298]
    L_sem_aA = [0.027886576355325265, 0.032219735289080446, 0.0317402498792375, 0.031983744272522235, 0.025596617823257447, 0.03176628999248938, 0.027504092933868067, 0.02601922876735056, 0.026912508955285398, 0.033254640426640426, 0.024136681794068523, 0.028975486094600995, 0.028580157686190683, 0.024329819775707664, 0.02659083004861742, 0.031527622831685705, 0.027065309790798284, 0.02154382198806116, 0.027359710117500578, 0.02643509072347617, 0.017748301931798133, 0.028188180769118026, 0.02391630924525466, 0.011347798645646938, 0.028188180769118026, 0.02590029726752835, 0.011840416534822372]
    
    L_error_base1_aA= []
    L_error_base2_aA= []
    L_error_base3_aA =[]
    L_error_sem_base1_aA= []
    L_error_sem_base2_aA= []
    L_error_sem_base3_aA =[]
    L_errors_base1 = []
    L_errors_base2 = []
    L_errors_base3 = []
    L_errors_base4 = []
    L_errors_base5 = []
    L_errors_sem_base1 = []
    L_errors_sem_base2 = []
    L_errors_sem_base3 = []
    L_errors_sem_base4 = []
    L_errors_sem_base5 = []
    for i in range (0,len(L_repulsion)):
        if L_repulsion[i] ==-1:
            L_repulsion[i] = 10
        if L_repulsion[i] == 20:
            L_repulsion[i] = 8
        L_error_base1_aA.append(L_error_high_bb_aA[3*i])
        L_error_base2_aA.append(L_error_high_bb_aA[3*i+1])
        L_error_base3_aA.append(L_error_high_bb_aA[3*i+2])
        L_error_sem_base1_aA.append(L_sem_aA[3*i])
        L_error_sem_base2_aA.append(L_sem_aA[3*i+1])
        L_error_sem_base3_aA.append(L_sem_aA[3*i+2])
        L_errors_base1.append(L_error_high_bb[5*i])
        L_errors_base2.append(L_error_high_bb[5*i+1])
        L_errors_base3.append(L_error_high_bb[5*i+2])
        L_errors_base4.append(L_error_high_bb[5*i+3])
        L_errors_base5.append(L_error_high_bb[5*i+4])
        L_errors_sem_base1.append(L_SEM_error_high_g_bb[5*i])
        L_errors_sem_base2.append(L_SEM_error_high_g_bb[5*i+1])
        L_errors_sem_base3.append(L_SEM_error_high_g_bb[5*i+2])
        L_errors_sem_base4.append(L_SEM_error_high_g_bb[5*i+3])
        L_errors_sem_base5.append(L_SEM_error_high_g_bb[5*i+4])
    if bases_str[0] =="aATLambdaSigma":
        bases_str[0] ="aAT$\Lambda\Sigma$"
    if bases_str[0] =="ATLambdaSigma":
        bases_str[0] ="AT$\Lambda\Sigma$"
    if bases_str[1] =="ATdeltabeta":
        bases_str[1] ="AT$\delta\\beta$"

    graph1 = ax.errorbar(L_repulsion, L_errors_base1,yerr=L_errors_sem_base1,label = "%s" %(base_legend[0]),fmt='-o', color =colors_bases[0], linewidth = lineswidth, markersize=lineswidth*3)
    graph2 = ax.errorbar(L_repulsion, L_errors_base2, yerr=L_errors_sem_base2,label = "%s" %(base_legend[1]),fmt='-o',color =colors_bases[1], linewidth = lineswidth, markersize=lineswidth*3)
    graph3 = ax.errorbar(L_repulsion, L_errors_base3, yerr=L_errors_sem_base3, label = "%s" %(base_legend[2]),fmt='-o',color =colors_bases[2], linewidth = lineswidth, markersize=lineswidth*3)
    graph4 = ax.errorbar(L_repulsion, L_errors_base4, yerr=L_errors_sem_base4,label = "%s" %(base_legend[3]),fmt='-o',color =colors_bases[3], linewidth = lineswidth, markersize=lineswidth*3)
    graph5 = ax.errorbar(L_repulsion, L_errors_base5, yerr=L_errors_sem_base5, label = "%s" %(base_legend[4]),fmt='-o',color =colors_bases[4], linewidth = lineswidth, markersize=lineswidth*3)
    
    # graph1aA = ax.errorbar(L_repulsion, L_error_base1_aA, label = "%s" %(L_bases_aA[0]),fmt='--',color =colors_bases[0], linewidth = 0.5*lineswidth,alpha=1)
    # graph2aA = ax.errorbar(L_repulsion, L_error_base2_aA, label = "%s" %(L_bases_aA[1]),fmt='--',color =colors_bases[2], linewidth = 0.5*lineswidth, alpha=1)
    # graph3aA= ax.errorbar(L_repulsion, L_error_base3_aA, label = "%s" %(L_bases_aA[2]),fmt='--',color =colors_bases[4], linewidth = 0.5*lineswidth, alpha=1)
    # graph = ax2.errorbar(L_repulsion, L_error_base1_aA, label = "%s" %(L_bases_aA[0]),fmt='--',color =colors_bases[0], linewidth = 0.5*lineswidth,alpha=1)
    # graph = ax2.errorbar(L_repulsion, L_error_base2_aA, label = "%s" %(L_bases_aA[1]),fmt='--',color =colors_bases[2], linewidth = 0.5*lineswidth, alpha=1)
    # graph = ax2.errorbar(L_repulsion, L_error_base3_aA, label = "%s" %(L_bases_aA[2]),fmt='--',color =colors_bases[4], linewidth = 0.5*lineswidth, alpha=1)
    
    
    graph = ax2.errorbar(L_repulsion, L_errors_base1, yerr=L_errors_sem_base1, label = "%s" %(base_legend[0]),fmt='-o',color =colors_bases[0], linewidth = lineswidth, markersize=lineswidth*3)
    graph = ax2.errorbar(L_repulsion, L_errors_base2,yerr=L_errors_sem_base2, label = "%s" %(base_legend[1]),fmt='-o',color =colors_bases[1], linewidth = lineswidth, markersize=lineswidth*3)
    graph = ax2.errorbar(L_repulsion, L_errors_base3, yerr=L_errors_sem_base3,label = "%s" %(base_legend[2]),fmt='-o',color =colors_bases[2], linewidth = lineswidth, markersize=lineswidth*3)
    graph = ax2.errorbar(L_repulsion, L_errors_base4, yerr=L_errors_sem_base4, label = "%s" %(base_legend[3]),fmt='-o',color =colors_bases[3], linewidth = lineswidth, markersize=lineswidth*3)
    graph = ax2.errorbar(L_repulsion, L_errors_base5,yerr=L_errors_sem_base5, label = "%s" %(base_legend[4]),fmt='-o',color =colors_bases[4], linewidth = lineswidth, markersize=lineswidth*3)
    
    ax.set_ylim([0,0.6])
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
    labels = ["0","1","2","3","4","5","6","20", "$\infty$"]
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
    ax.set_position([box.x0, box.y0+box.height*0.2, box.width, box.height*0.8])
    ax2.set_position([box.x0+box.width +2*d, box.y0+box.height*0.2, box.width, box.height*0.8])
    # ax.legend()
    # l1 =ax.legend([graph1,graph2,graph3,graph4, graph5], base_legend, loc='lower center', bbox_to_anchor=(0.5, -0.39), ncol=3, fontsize =fontsizee)
    # ax.add_artist(l1)
    # ax2.legend([graph1aA,graph2aA,graph3aA],L_bases_aA, loc='lower center', bbox_to_anchor=(0.78, -0.39), ncol=2, fontsize =fontsizee)
    # # gca().add_artist(l1)
    l1 =ax.legend([graph1, graph2,graph3, graph4,graph5], base_legend, loc='lower center', bbox_to_anchor=(1, -0.39), ncol=3, fontsize =fontsizee)
    ax.add_artist(l1)
    # ax2.legend([graph2aA], [L_bases_aA[1]], loc='lower center', bbox_to_anchor=(0.5, -0.39), ncol=2, fontsize =fontsizee)
    # gca().add_artist(l1)
    return

def graph_error_high_g_bb_all_bases_5_small(L_error_high_bb, L_SEM_error_high_g_bb, L_repulsion, bases_str, colors_bases, fontsizee, lineswidth, base_legend):
    f, ax = plt.subplots()
    # f.set_figheight(5)
    # f.set_figwidth(5)
    L_bases_aA = ["Opposite group", "Equal", "Same group"]#['aAT$\\Lambda\\Sigma$', 'aATGC', 'aATXK']
    L_error_high_bb_aA = [0.5218136237623762, 0.519714495049505, 0.5181401485148515, 0.4749326435643565, 0.45142219801980193, 0.49742857425742565, 0.44886822772277224, 0.36549699009900993, 0.4503376336633664, 0.4443550693069307, 0.2931462475247524, 0.36798096039603967, 0.43221495049504954, 0.24682503960396043, 0.24857436633663366, 0.43868737623762377, 0.2259735247524752, 0.1502291287128713, 0.4294860792079208, 0.22464403960396043, 0.09033335742574258, 0.4308855049504951, 0.21943115841584157, 0.04639121584158416, 0.4325648613861387, 0.22020082178217823, 0.04842039702970298]
    L_sem_aA = [0.027886576355325265, 0.032219735289080446, 0.0317402498792375, 0.031983744272522235, 0.025596617823257447, 0.03176628999248938, 0.027504092933868067, 0.02601922876735056, 0.026912508955285398, 0.033254640426640426, 0.024136681794068523, 0.028975486094600995, 0.028580157686190683, 0.024329819775707664, 0.02659083004861742, 0.031527622831685705, 0.027065309790798284, 0.02154382198806116, 0.027359710117500578, 0.02643509072347617, 0.017748301931798133, 0.028188180769118026, 0.02391630924525466, 0.011347798645646938, 0.028188180769118026, 0.02590029726752835, 0.011840416534822372]
    
    L_error_base1_aA= []
    L_error_base2_aA= []
    L_error_base3_aA =[]
    L_error_sem_base1_aA= []
    L_error_sem_base2_aA= []
    L_error_sem_base3_aA =[]
    L_errors_base1 = []
    L_errors_base2 = []
    L_errors_base3 = []
    L_errors_base4 = []
    L_errors_base5 = []
    L_errors_sem_base1 = []
    L_errors_sem_base2 = []
    L_errors_sem_base3 = []
    L_errors_sem_base4 = []
    L_errors_sem_base5 = []
    for i in range (0,len(L_repulsion)-3):
        if L_repulsion[i] ==-1:
            L_repulsion[i] = 10
        if L_repulsion[i] == 20:
            L_repulsion[i] = 8
        L_error_base1_aA.append(L_error_high_bb_aA[3*i])
        L_error_base2_aA.append(L_error_high_bb_aA[3*i+1])
        L_error_base3_aA.append(L_error_high_bb_aA[3*i+2])
        L_error_sem_base1_aA.append(L_sem_aA[3*i])
        L_error_sem_base2_aA.append(L_sem_aA[3*i+1])
        L_error_sem_base3_aA.append(L_sem_aA[3*i+2])
        L_errors_base1.append(L_error_high_bb[5*i])
        L_errors_base2.append(L_error_high_bb[5*i+1])
        L_errors_base3.append(L_error_high_bb[5*i+2])
        L_errors_base4.append(L_error_high_bb[5*i+3])
        L_errors_base5.append(L_error_high_bb[5*i+4])
        L_errors_sem_base1.append(L_SEM_error_high_g_bb[5*i])
        L_errors_sem_base2.append(L_SEM_error_high_g_bb[5*i+1])
        L_errors_sem_base3.append(L_SEM_error_high_g_bb[5*i+2])
        L_errors_sem_base4.append(L_SEM_error_high_g_bb[5*i+3])
        L_errors_sem_base5.append(L_SEM_error_high_g_bb[5*i+4])
    if bases_str[0] =="aATLambdaSigma":
        bases_str[0] ="aAT$\Lambda\Sigma$"
    if bases_str[0] =="ATLambdaSigma":
        bases_str[0] ="AT$\Lambda\Sigma$"
    if bases_str[1] =="ATdeltabeta":
        bases_str[1] ="AT$\delta\\beta$"
    L_repulsion =L_repulsion[0:6]

    graph1 = ax.errorbar(L_repulsion, L_errors_base1,yerr=L_errors_sem_base1,label = "%s" %(base_legend[0]),fmt='-o', color =colors_bases[0], linewidth = lineswidth, markersize=lineswidth*3)
    graph2 = ax.errorbar(L_repulsion, L_errors_base2, yerr=L_errors_sem_base2,label = "%s" %(base_legend[1]),fmt='-o',color =colors_bases[1], linewidth = lineswidth, markersize=lineswidth*3)
    graph3 = ax.errorbar(L_repulsion, L_errors_base3, yerr=L_errors_sem_base3, label = "%s" %(base_legend[2]),fmt='-o',color =colors_bases[2], linewidth = lineswidth, markersize=lineswidth*3)
    graph4 = ax.errorbar(L_repulsion, L_errors_base4, yerr=L_errors_sem_base4,label = "%s" %(base_legend[3]),fmt='-o',color =colors_bases[3], linewidth = lineswidth, markersize=lineswidth*3)
    graph5 = ax.errorbar(L_repulsion, L_errors_base5, yerr=L_errors_sem_base5, label = "%s" %(base_legend[4]),fmt='-o',color =colors_bases[4], linewidth = lineswidth, markersize=lineswidth*3)
    
    graph1aA = ax.errorbar(L_repulsion, L_error_base1_aA, label = "%s" %(L_bases_aA[0]),fmt='--',color =colors_bases[0], linewidth = 0.5*lineswidth,alpha=1)
    graph2aA = ax.errorbar(L_repulsion, L_error_base2_aA, label = "%s" %(L_bases_aA[1]),fmt='--',color =colors_bases[2], linewidth = 0.5*lineswidth, alpha=1)
    graph3aA= ax.errorbar(L_repulsion, L_error_base3_aA, label = "%s" %(L_bases_aA[2]),fmt='--',color =colors_bases[4], linewidth = 0.5*lineswidth, alpha=1)
    ax.set_ylim(0,0.6)
    ax.set_xlim(-0.5, 5.5)

    ax.spines['right'].set_linewidth(lineswidth*0.5)
    ax.spines['left'].set_linewidth(lineswidth*0.5)
    
    ax.spines['bottom'].set_linewidth(lineswidth*0.5)
    
    ax.spines['top'].set_linewidth(lineswidth*0.5)
    
    ax.yaxis.tick_left()
    # labels = ["0","1","2","3","4","5","6","20", "$\infty$"]
    # ax.set_xticks([0,1,2,3,4,5])
    # ax.set_xticklabels(labels[0:6],fontsize =fontsizee)
    
    ax.tick_params(axis='y', labelsize=fontsizee, width=lineswidth*0.5)
    
    ax.tick_params(axis='x', labelsize=fontsizee, width=lineswidth*0.5)
    
    ax.set_xlabel("Same group matches cost", fontsize =fontsizee)
    # ax.xaxis.set_label_coords(1.0, -.1)
    ax.set_ylabel("Error probability $\epsilon$", fontsize =fontsizee)
    # ax2.set_ylabel("Error probability $\epsilon$", fontsize =fontsizee)
    # d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    # kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    # ax.plot((1-d, 1+d), (-d, +d), **kwargs)
    # ax.plot((1-d, 1+d), (1-d, 1+d), **kwargs)


    box = ax.get_position()
    ax.set_position([box.x0, box.y0+box.height*0.2, box.width*0.8, box.height*0.8])
    # ax.legend()
    l1 =ax.legend([graph1,graph2,graph3,graph4, graph5], base_legend, loc='lower center', bbox_to_anchor=(0.15, -0.449), ncol=2, fontsize =fontsizee)
    ax.add_artist(l1)
    ax.legend([graph1aA,graph2aA,graph3aA],L_bases_aA, loc='lower center', bbox_to_anchor=(0.8, -0.37), ncol=2, fontsize =fontsizee)
    
 
    return

def probability_base_base(data):
    matrix = np.array([[0,0,0,0,0,0,0,0,0,0]])
    for i in range(0,len(data)):
        if i % 4 ==0:
            row = data.iloc[i+2]

            pur1pur1= float(row["pur1pur1"])
            pur1pur2= float(row["pur1pur2"])
            pur2pur2= float(row["pur2pur2"])
            py1py1= float(row["py1py1"])
            py1py2= float(row["py1py2"])
            py2py2= float(row["py2py2"])
            py1pur1= float(row["py1pur1"])
            py1pur2= float(row["py1pur2"])
            py2pur1= float(row["py2pur1"])
            py2pur2= float(row["py2pur2"])
            total = pur1pur1 +pur1pur2 +pur2pur2 +py1py1 +py1py2 +py2py2 +py1pur1 + py1pur2 +py2pur1 +py2pur2
            list_values = [[pur1pur1, pur1pur2, pur2pur2, py1py1, py1py2, py2py2, py1pur1, py1pur2 ,py2pur1 ,py2pur2]]
            array_prob = np.array(list_values)/total
            matrix = np.concatenate((matrix, array_prob), axis=0)
    matrix = np.delete(matrix, 0, 0)
    prob_mean = np.mean(matrix, axis=0)
    prob_sem = np.std(matrix, axis=0)/(np.sqrt(len(data)))
    return prob_mean, prob_sem

def cumulative_correct_incorrect_probs(prob_mean):
    pur_pur_incorrect = prob_mean[0] + prob_mean[1] + prob_mean[2]
    py_py_incorrect = prob_mean[3] + prob_mean[4] + prob_mean[5]
    pur_py_incorrect = prob_mean[7] + prob_mean[8] 
    pur_py_correct = prob_mean[6] + prob_mean[9] 
    # [[pur1pur1, pur1pur2, pur2pur2, py1py1, py1py2, py2py2, py1pur1, py1pur2 ,py2pur1 ,py2pur2]]
    return [pur_pur_incorrect,py_py_incorrect, pur_py_incorrect, pur_py_correct]

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

def matrix_purpy_values(probs,bases):
    # pur1pur1, pur1pur2, pur2pur2, py1py1, py1py2, py2py2, py1pur1, py1pur2 ,py2pur1 ,py2pur2 = probs
    matrix = np.zeros((len(bases), len(bases)))
    matrix[0,0] = probs[0] #pur1pur1
    matrix[0,1] = probs[6] #py1pur1
    matrix[0,2] = probs[1] #pur1pur2
    matrix[0,3] = probs[8] #py2pur1

    matrix[1,0] = probs[6] #py1pur1
    matrix[1,1] = probs[3] #py1py1
    matrix[1,2] = probs[7] #py1pur2
    matrix[1,3] = probs[4] #py1py2

    matrix[2,0] = probs[1] #pur1pur2
    matrix[2,1] = probs[7] #py1pur2
    matrix[2,2] = probs[2] #pur2pur2
    matrix[2,3] = probs[9] #py2pur2

    matrix[3,0] = probs[8] #py2pur1
    matrix[3,1] = probs[4] #py1py2
    matrix[3,2] = probs[9] #py2pur2
    matrix[3,3] = probs[5] #py2py2
    return matrix

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

def graph_prob_vs_OHstrength(base_combination,L_prob, L_prob_sem, OH_matrix, delta_g_pol, repulsion,color):
    pur1pur1_strength = OH_matrix[0,0]
    pur1pur2_strength = OH_matrix[0,2]
    pur2pur2_strength = OH_matrix[2,2] 
    py1py1_strength = OH_matrix[1,1]
    py1py2_strength = OH_matrix[1,3]
    py2py2_strength = OH_matrix[3,3]
    py1pur1_strength = OH_matrix[1,0]
    py1pur2_strength = OH_matrix[1,2]
    py2pur1_strength = OH_matrix[3,0]
    py2pur2_strength = OH_matrix[3,2]
    L_OH_strength = [pur1pur1_strength, pur1pur2_strength, pur2pur2_strength, py1py1_strength, py1py2_strength, py2py2_strength, py1pur1_strength, py1pur2_strength ,py2pur1_strength ,py2pur2_strength]
    pur1,py1,pur2,py2 = base_combination[0], base_combination[1], base_combination[2],base_combination[3]
    L_labels = [pur1+pur1, pur1+pur2, pur2+pur2, py1+py1, py1+py2, py2+py2, py1+pur1, py1+pur2 ,py2+pur1 ,py2+pur2]
    markers = ["1", "+", "2", "3", "x", "4", "v", "<" ,">", "^"]
    begin =0
    if repulsion == "-1":
        begin = 6
    plt.title("For $\Delta G_{pol}=$%s, repulsion =%s: probability of base-base bond as function of the bond strength" %(delta_g_pol, repulsion))
    for i in range(0, len(L_labels)):
        plt.annotate(L_labels[i], (L_OH_strength[i], L_prob[i]))
    for i in range(begin, len(L_OH_strength)):
        plt.errorbar(L_OH_strength[i], L_prob[i], yerr= L_prob_sem[i] ,label = L_labels[i], marker=markers[i],c=color)
        # plt.scatter(L_OH_strength[i], L_prob[i], label = L_labels[i], marker=markers[i], c=color)
    plt.xlabel("$\Delta G_{OH}$ (base-base)")
    # plt.yscale('log')
    plt.ylabel("Probability (base-base)")
    plt.legend()
    return




plot1 = False
plot2 = False
plot3 = False
plot4 = True

if plot4 ==True:
    L_delta_g_pol_label = ["-ln4","-ln2", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    L_delta_g_pol_number = [-np.log(4),-np.log(2), 0,1,2,3,4,5,6,7,8,9,10]
    L_errors = [[0.5237129365887834, 0.5237129365887834, 0.5237129365887834, 0,0], [0.476843858505417, 0.45238679848588237, 0.4983239357157505,0,0], [0.4518352218724654, 0.36374161504240043, 0.4518352218724654,0,0], [0.4410076540717165, 0.29081265389302524, 0.3659956846194553,0,0], [0.43675880562301894, 0.2498544930193637, 0.2498544930193637,0,0], [0.4351570506859479, 0.23169413106203252, 0.14913661237245746,0,0], [0.43456241422467223, 0.22449770354088439, 0.09066412803169821,0,0], [0.43421500232773125, 0.2201702845430179, 0.04742591271129115,0,0],[0,0,0,0,0]]

    all_steps = True
    # DIR = "230414output"
    L_DIR = ["230529output", "230529_equi_output"]#"230517_equi_output","230517output"]
    model = "muriel_big_abs_BF_x1"
    model_title = "BF x=1"
    bases_good=[["A","T","$\delta$","$\\beta$"],["A", "T", "$\Lambda$", "$\Sigma$"],["A","T","iG","iC"],["A","T","G","C"],["A","T","X","K"]]
    bases_pur1py1pur2py2=[["A","T","delta","beta"],["A", "T", "L", "S"],["A","T","iG","iC"],["A","T","G","C"],["A","T","X","K"]]#[["A", "T", "L", "S"],["A","T","G","C"],["A","T","X","K"]]#
    base_combinations =["ATdeltabeta","ATLambdaSigma","ATiGiC","ATGC", "ATXK"] #["ATLambdaSigma","ATGC", "ATXK"]# "GCXK","aATXK","aATGC", "deltabetaXK",
    Base_legend = ["Opposite group S", "Opposite group O", "Equal S", "Equal O", "Same group S"]
    # Base_legend = base_combinations
    ratio ='1111'
    # number_steps = 0
    # MX_models = ["x=$e^{\Delta G_{TT} +5}$, Model B","x=$e^{\Delta G_{TT} +5}$, Model F", "x=1, Model B", "x=1, Model F" ]

    # L_models, L_models_titles = make_model_names(base_combinations, model, model_title)exit()\e
    L_models = [model]
    L_models_titles = [model_title]
    date = date_now(datetime.datetime.now())
    date_options = ["2023530", "2023531", "2023621", "2023625"]#"2023517","2023518", "2023519", "2023524"]
    L_repulsion = ["0","1","2","3","4", "5", "6","20", "-1"]
    number_tables=5

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
            if number_tables == 5:
                # fig = plt.figure()
                # fig.set_figheight(6)
                # fig.set_figwidth(6)       
                # graph = plt.subplot2grid(shape=(9, 9),  loc=(0, 0), colspan=3, rowspan=1)
                # table1 = plt.subplot2grid(shape=(3, 3), loc=(1, 0), colspan=1, rowspan=1)
                # table1 = plt.subplot2grid(shape=(3, 3), loc=(1, 0), colspan=1, rowspan=1)
                # table2 = plt.subplot2grid(shape=(3, 3), loc=(1, 1), colspan=1, rowspan=1)
                # table3 = plt.subplot2grid(shape=(3, 3), loc=(1, 2), colspan=1, rowspan=1)
                # table4 = plt.subplot2grid(shape=(3, 3), loc=(2, 0), colspan=1, rowspan=1)
                # table5 = plt.subplot2grid(shape=(3, 3), loc=(2, 1), colspan=1, rowspan=1) #9old


                # Create Blank Figure
                graph = plt.subplot2grid(shape=(3, 4), loc=(0, 0), colspan=2, rowspan=4)
                table1 = plt.subplot2grid(shape=(3, 4), loc=(0, 2), colspan=1)
                table2 = plt.subplot2grid(shape=(3, 4), loc=(1, 2), colspan=1)
                table3 = plt.subplot2grid(shape=(3, 4), loc=(2, 2), colspan=1)
                table4 = plt.subplot2grid(shape=(3, 4), loc=(0, 3), colspan=1)
                table5 = plt.subplot2grid(shape=(3, 4), loc=(1, 3), colspan=1)
                
                #eronder werkt niet
                # graph = plt.subplot2grid(shape=(5, 3), loc=(0, 0), colspan=3, rowspan=3 )
                # table1 = plt.subplot2grid(shape=(5, 3), loc=(3, 0), colspan=1)
                # table2 = plt.subplot2grid(shape=(5, 3), loc=(3, 1), colspan=1)
                # table3 = plt.subplot2grid(shape=(5, 3), loc=(3, 2), colspan=1)
                # table4 = plt.subplot2grid(shape=(5, 3), loc=(4, 0), colspan=1)
                # table5 = plt.subplot2grid(shape=(5, 3), loc=(4, 1), colspan=1)
            if number_tables ==1:
                plt.figure()
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
                            
                            string = '/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/' + DIR + '/' +date + 'model_' + model + base +ratio + '_delta_g_pol_' + label+ "_repulsion" + repulsion +'_allsteps.csv'
                            file_exists = os.path.exists(string)
                            # print(string)
                            if file_exists == True:
                                print(string)
                                data = pd.read_csv(string)
                                epsilon_mean = data['error probability'].mean()
                                epsilon_sem = data['error probability'].std()#/(np.sqrt(len(data['error probability'])))
                                print(epsilon_mean, epsilon_sem)
                                L_sem.append(epsilon_sem)
                                L_epsilon_mean.append(epsilon_mean)
                                L_epsilon_std_high.append(epsilon_mean + epsilon_sem)
                                L_epsilon_std_low.append(epsilon_mean - epsilon_sem)
                                L_delta_g_bb_graph.append(delta_g_pol)
                                
                            else:
                                print("file does not exist")
                                print(string)
                            
                            string2 = '/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/' + DIR + '/' +date + 'polymer_strings_model_' + model + base +ratio + '_delta_g_pol_' + label+ "_repulsion" + repulsion +'_allsteps.csv'
                            file_exists2 = os.path.exists(string2)
                            # print(string2)
                            if file_exists2 == True:
                                data2 = pd.read_csv(string2,on_bad_lines='skip')
                                prob_mean, prob_std = probability_base_base(data2)
                                matrix_prob= matrix_purpy_values(prob_mean,bases_pur1py1pur2py2[j])
                                matrix_sem = matrix_purpy_values(prob_std,bases_pur1py1pur2py2[j])
                                print(matrix_prob)
                                print(matrix_sem)
                                # print(matrix_purpy_values(prob_std,bases_pur1py1pur2py2[j]))
                                
                                repulsion_int = float(repulsion)
                                print(bases_pur1py1pur2py2[j], "repulsion:", repulsion_int, "G_pol:",delta_g_pol)
                                # print(matrix_purpy(bases_pur1py1pur2py2[j]))
                                OH_matrix, matrix = make_E_OH_matrix(bases_pur1py1pur2py2[j], repulsion_int)
                                print(OH_matrix)
                                
                            else:
                                print("polymer string file does not exist")
                            
                

                if number_tables ==1:
                    graph_error_g_bb(L_delta_g_bb_graph, L_epsilon_mean, L_models_titles[k], base, L_epsilon_std_high, L_epsilon_std_low, color_list[j*23], repulsion)
                # define error as all mismatches and plot the error as function of delta G_bb
                if number_tables ==5:
                    L_colors.append(color_list[j*23])
                    print(L_delta_g_bb_graph, L_epsilon_mean, L_epsilon_mean[10])    
                    L_error_high_bb.append(L_epsilon_mean[10])
                    L_sem_error_high_g_bb.append(epsilon_sem)
                    graph_error_g_bb_5tables(L_delta_g_bb_graph, L_epsilon_mean, L_models_titles[k], base, L_epsilon_std_high, L_epsilon_std_low, color_list[j*23],repulsion, matrix, bases_good[j],j,graph,table1,table2,table3, table4, table5, L_errors[r][j],fontsizee, line_widths, Base_legend[j])
                    print(L_errors[r][j])
            # plt.show()
            L_rep_numbers.append(rep)

        graph_error_high_g_bb_all_bases_5(L_error_high_bb, L_sem_error_high_g_bb, L_rep_numbers, base_combinations, L_colors, fontsizee, line_widths, Base_legend)   
        plt.show()
        # graph_error_high_g_bb_all_bases_5_small(L_error_high_bb, L_sem_error_high_g_bb, L_rep_numbers, base_combinations, L_colors, fontsizee, line_widths, Base_legend) 
        # plt.show()

if plot1 ==True:
    L_delta_g_pol_label = ["-ln4","-ln2", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    L_delta_g_pol_number = [-np.log(4),-np.log(2), 0,1,2,3,4,5,6,7,8,9,10]
    L_errors = [[0.5237129365887834, 0.5237129365887834, 0.5237129365887834], [0.476843858505417, 0.45238679848588237, 0.4983239357157505], [0.4518352218724654, 0.36374161504240043, 0.4518352218724654], [0.4410076540717165, 0.29081265389302524, 0.3659956846194553], [0.43675880562301894, 0.2498544930193637, 0.2498544930193637], [0.4351570506859479, 0.23169413106203252, 0.14913661237245746], [0.43456241422467223, 0.22449770354088439, 0.09066412803169821], [0.43421500232773125, 0.2201702845430179, 0.04742591271129115],[0,0,0]]

    all_steps = True
    # DIR = "230414output"
    L_DIR = ["230529output", "230529_equi_output"]#"230517_equi_output","230517output"]
    model = "muriel_big_abs_BF_x1"
    model_title = "BF x=1"
    bases_good=[["aA", "T", "$\Lambda$", "$\Sigma$"],["aA","T","G","C"],["aA","T","X","K"]]
    bases_pur1py1pur2py2=[["aA", "T", "L", "S"],["aA","T","G","C"],["aA","T","X","K"]]#[["A", "T", "L", "S"],["A","T","G","C"],["A","T","X","K"]]#
    base_combinations = ["aATLambdaSigma","aATGC", "aATXK"] #["ATLambdaSigma","ATGC", "ATXK"]# "GCXK","aATXK","aATGC", "deltabetaXK",
    Base_legend = ["Opposite group", "Equal", "Same group"]
    Base_legend = base_combinations
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
            if number_tables == 3:
                fig = plt.figure()
                fig.set_figheight(4)
                fig.set_figwidth(6)       
                graph = plt.subplot2grid(shape=(3, 3), loc=(0, 0), colspan=2, rowspan=3)
                table1 = plt.subplot2grid(shape=(3, 3), loc=(0, 2), colspan=1)
                table2 = plt.subplot2grid(shape=(3, 3), loc=(1, 2), colspan=1)
                table3 = plt.subplot2grid(shape=(3, 3), loc=(2, 2), colspan=1)
            if number_tables == 6:
                fig = plt.figure()
                # fig.set_figheight(4)
                # fig.set_figwidth(6)       
                graph = plt.subplot2grid(shape=(3, 4), loc=(0, 0), colspan=2, rowspan=3)
                table1 = plt.subplot2grid(shape=(3, 4), loc=(0, 2), colspan=1)
                table2 = plt.subplot2grid(shape=(3, 4), loc=(1, 2), colspan=1)
                table3 = plt.subplot2grid(shape=(3, 4), loc=(2, 2), colspan=1)
                table4 = plt.subplot2grid(shape=(3, 4), loc=(0, 3), colspan=1)
                table5 = plt.subplot2grid(shape=(3, 4), loc=(1, 3), colspan=1)
                table6 = plt.subplot2grid(shape=(3, 4), loc=(2, 3), colspan=1)
            if number_tables ==1:
                plt.figure()
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
                            
                            string = '/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/' + DIR + '/' +date + 'model_' + model + base +ratio + '_delta_g_pol_' + label+ "_repulsion" + repulsion +'_allsteps.csv'
                            file_exists = os.path.exists(string)
                            # print(string)
                            if file_exists == True:
                                print(string)
                                data = pd.read_csv(string)
                                epsilon_mean = data['error probability'].mean()
                                epsilon_sem = data['error probability'].std()#/(np.sqrt(len(data['error probability'])))
                                print(epsilon_mean, epsilon_sem)
                                L_sem.append(epsilon_sem)
                                L_epsilon_mean.append(epsilon_mean)
                                L_epsilon_std_high.append(epsilon_mean + epsilon_sem)
                                L_epsilon_std_low.append(epsilon_mean - epsilon_sem)
                                L_delta_g_bb_graph.append(delta_g_pol)
                                
                            else:
                                print("file does not exist")
                                print(string)
                            
                            string2 = '/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/' + DIR + '/' +date + 'polymer_strings_model_' + model + base +ratio + '_delta_g_pol_' + label+ "_repulsion" + repulsion +'_allsteps.csv'
                            file_exists2 = os.path.exists(string2)
                            # print(string2)
                            if file_exists2 == True:
                                data2 = pd.read_csv(string2,on_bad_lines='skip')
                                prob_mean, prob_std = probability_base_base(data2)
                                matrix_prob= matrix_purpy_values(prob_mean,bases_pur1py1pur2py2[j])
                                matrix_sem = matrix_purpy_values(prob_std,bases_pur1py1pur2py2[j])
                                print(matrix_prob)
                                print(matrix_sem)
                                # print(matrix_purpy_values(prob_std,bases_pur1py1pur2py2[j]))
                                
                                repulsion_int = float(repulsion)
                                print(bases_pur1py1pur2py2[j], "repulsion:", repulsion_int, "G_pol:",delta_g_pol)
                                # print(matrix_purpy(bases_pur1py1pur2py2[j]))
                                OH_matrix, matrix = make_E_OH_matrix(bases_pur1py1pur2py2[j], repulsion_int)
                                print(OH_matrix)
                                
                            else:
                                print("polymer string file does not exist")
                            
                

                if number_tables ==1:
                    graph_error_g_bb(L_delta_g_bb_graph, L_epsilon_mean, L_models_titles[k], base, L_epsilon_std_high, L_epsilon_std_low, color_list[j*23], repulsion)
                # define error as all mismatches and plot the error as function of delta G_bb
                if number_tables ==3:
                    L_colors.append(color_list[j*30])
                    print(L_delta_g_bb_graph, L_epsilon_mean, L_epsilon_mean[10])    
                    L_error_high_bb.append(L_epsilon_mean[10])
                    L_sem_error_high_g_bb.append(epsilon_sem)
                    graph_error_g_bb_3tables(L_delta_g_bb_graph, L_epsilon_mean, L_models_titles[k], base, L_epsilon_std_high, L_epsilon_std_low, color_list[j*30],repulsion, matrix, bases_good[j],j,graph,table1,table2,table3, L_errors[r][j],fontsizee, line_widths, Base_legend[j])
                    # print(L_errors[r][j])
                elif number_tables ==6:
                    graph_error_g_bb_6tables(L_delta_g_bb_graph, L_epsilon_mean, L_models_titles[k], base, L_epsilon_std_high, L_epsilon_std_low, color_list[j*30],repulsion, matrix_prob, bases_pur1py1pur2py2[j],j,graph,table1,table2,table3,table4,table5,table6,OH_matrix)
            plt.show()
            L_rep_numbers.append(rep)

        graph_error_high_g_bb_all_bases(L_error_high_bb, L_sem_error_high_g_bb, L_rep_numbers, base_combinations, L_colors, fontsizee, line_widths, Base_legend)
        print("base:",Base_legend)
        print("error", L_error_high_bb)
        print("sem error:", L_sem_error_high_g_bb)
        plt.show()


if plot2==True:
    L_delta_g_pol_label = ["-ln4","-ln2", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    L_delta_g_pol_number = [-np.log(4),-np.log(2), 0,1,2,3,4,5,6,7,8,9,10]

    all_steps = True
    # DIR = "230414output"
    L_DIR = ["230527output"]
    model = "muriel_big_abs_BF_x1"
    model_title = "BF x=1"
    bases_pur1py1pur2py2=[["aA", "T", "L", "S"],["aA","T","G","C"],["aA","T","X","K"]]#[["A", "T", "L", "S"],["A","T","G","C"],["A","T","X","K"]]#
    base_combinations = ["aATLambdaSigma","aATGC", "aATXK"] #["ATLambdaSigma","ATGC", "ATXK"]# "GCXK","aATXK","aATGC", "deltabetaXK",
    ratio ='1111'
    # number_steps = 0
    # MX_models = ["x=$e^{\Delta G_{TT} +5}$, Model B","x=$e^{\Delta G_{TT} +5}$, Model F", "x=1, Model B", "x=1, Model F" ]

    # L_models, L_models_titles = make_model_names(base_combinations, model, model_title)exit()\e
    L_models = [model]
    L_models_titles = [model_title]
    date = date_now(datetime.datetime.now())
    date_options = ["2023529"]
    L_repulsion = ["-1","0","1","2","3","4", "5", "6","20"]
    number_tables=3

    #define colors
    colormap = plt.cm.plasma #plt.cm.plasma #nipy_spectral, Set1,Paired   
    color_list = [colormap(i) for i in np.linspace(0, 1,100)]

    print(L_models)
    OH_matrix = np.zeros((4, 4))


    # line_styles = ['dashed', 'solid']

    # L_models = L_models[i:i+3] #for same M and X diffent models
    line_styles = ['solid', 'dashdot', 'dashed', 'dotted', (5, (10, 3)), (0, (3, 5, 1, 5))]
    print(L_models)



    for k, model in enumerate(L_models):
        for repulsion in L_repulsion:
            
            for j, base in enumerate(base_combinations):
                if number_tables == 3:
                    fig = plt.figure()
                    fig.set_figheight(4)
                    fig.set_figwidth(6)       
                    graph = plt.subplot2grid(shape=(3, 3), loc=(0, 0), colspan=2, rowspan=3)
                    # table1 = plt.subplot2grid(shape=(3, 3), loc=(0, 2), colspan=1)
                    table2 = plt.subplot2grid(shape=(3, 3), loc=(1, 2), colspan=1)
                    # table3 = plt.subplot2grid(shape=(3, 3), loc=(2, 2), colspan=1)
                L_delta_g_bb_graph = []
                L_mean_pur_pur_incorrect = []
                L_mean_py_py_incorrect = []
                L_mean_pur_py_incorrect = []
                L_mean_pur_py_correct = []
                for i,delta_g_pol in enumerate(L_delta_g_pol_number): 
                    rep = float(repulsion)
                    label = L_delta_g_pol_label[i]
                    for DIR in L_DIR:
                        for date in date_options:
                            
                            string = '/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/' + DIR + '/' +date + 'model_' + model + base +ratio + '_delta_g_pol_' + label+ "_repulsion" + repulsion +'_allsteps.csv'
                            file_exists = os.path.exists(string)
                            # print(string)
                            if file_exists == True:
                                print(string)
                                data = pd.read_csv(string)
                                epsilon_mean = data['error probability'].mean()
                                epsilon_sem = data['error probability'].std()#/(np.sqrt(len(data['error probability'])))
                                print(epsilon_mean, epsilon_sem)
                                L_delta_g_bb_graph.append(delta_g_pol)
                                
                            else:
                                print("file does not exist")
                                print(string)
                            
                            string2 = '/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/' + DIR + '/' +date + 'polymer_strings_model_' + model + base +ratio + '_delta_g_pol_' + label+ "_repulsion" + repulsion +'_allsteps.csv'
                            file_exists2 = os.path.exists(string2)
                            # print(string2)
                            if file_exists2 == True:
                                data2 = pd.read_csv(string2,on_bad_lines='skip')
                                prob_mean, prob_std = probability_base_base(data2)
                                array_pur_py_probs = cumulative_correct_incorrect_probs(prob_mean) #[pur_pur_incorrect,py_py_incorrect, pur_py_incorrect, pur_py_correct]
                                L_mean_pur_pur_incorrect.append(array_pur_py_probs[0])
                                L_mean_py_py_incorrect.append(array_pur_py_probs[1])
                                L_mean_pur_py_incorrect.append(array_pur_py_probs[2])
                                L_mean_pur_py_correct.append(array_pur_py_probs[3])
                                print("append", delta_g_pol)
                                
                                matrix_prob= matrix_purpy_values(prob_mean,bases_pur1py1pur2py2[j])
                                matrix_sem = matrix_purpy_values(prob_std,bases_pur1py1pur2py2[j])
                                print(matrix_prob)
                                print(matrix_sem)
                                # print(matrix_purpy_values(prob_std,bases_pur1py1pur2py2[j]))
                                
                                repulsion_int = float(repulsion)
                                print(bases_pur1py1pur2py2[j], "repulsion:", repulsion_int, "G_pol:",delta_g_pol)
                                # print(matrix_purpy(bases_pur1py1pur2py2[j]))
                                OH_matrix, matrix = make_E_OH_matrix(bases_pur1py1pur2py2[j], repulsion_int)
                                print(OH_matrix)
                            else:
                                print("polymer string file does not exist")
                            

                    

                #define error as all mismatches and plot the error as function of delta G_bb
                if number_tables ==3:
                    print(L_delta_g_bb_graph, L_mean_pur_pur_incorrect)
                    graph_cummulative_3tables(L_delta_g_bb_graph, L_mean_pur_pur_incorrect,L_mean_py_py_incorrect,L_mean_pur_py_incorrect,L_mean_pur_py_correct, L_models_titles[k], base, repulsion, OH_matrix, bases_pur1py1pur2py2[j],graph,table2)
                plt.show()

















if plot3==True:
    L_delta_g_pol_label = ["-ln4","-ln2", "1", "10"]#["-ln4","-ln2", "0", "1", "2","3","4", "5","6", "7", "8" ,"9", "10"]
    all_steps = True
    # DIR = "230331output"
    L_DIR = ["230414output","230414output_test"]
    model = "muriel_big_abs_BF_x1"
    model_title = "BF x=1"
    bases_pur1py1pur2py2=[["iG", "iC", "X", "K"],["aA","T","G","C"],["aA","T","X","K"]]#["delta","beta","X","K"]]
    base_combinations = ["iGiCXK","aATGC", "aATXK"] #"GCXK","aATXK","aATGC", "deltabetaXK",
    ratio ='1111'
    # number_steps = 0
    # MX_models = ["x=$e^{\Delta G_{TT} +5}$, Model B","x=$e^{\Delta G_{TT} +5}$, Model F", "x=1, Model B", "x=1, Model F" ]

    # L_models, L_models_titles = make_model_names(base_combinations, model, model_title)
    L_models = [model]
    L_models_titles = [model_title]
    date = date_now(datetime.datetime.now())
    date_options = ["2023413","2023414", "2023415","2023416","2023417","2023418",date]
    L_repulsion = ["-1","0"]#,"1","2","3","4", "5", "6","20"]


    #define colors
    colormap = plt.cm.plasma #nipy_spectral, Set1,Paired   
    color_list = [colormap(i) for i in np.linspace(0, 1,100)]
    L_color_bases=["red","blue", "green"]

    print(L_models)



    # line_styles = ['dashed', 'solid']

    # L_models = L_models[i:i+3] #for same M and X diffent models
    line_styles = ['solid', 'dashdot', 'dashed', 'dotted', (5, (10, 3)), (0, (3, 5, 1, 5))]
    print(L_models)

    
    for k, model in enumerate(L_models):
        for i,delta_g_pol in enumerate(L_delta_g_pol_number):
        
            for repulsion in L_repulsion:
                L_epsilon_mean = []
                L_epsilon_std_high = []
                L_epsilon_std_low = []
                                     
                for j, base in enumerate(base_combinations):
                    for date in date_options:
                        for DIR in L_DIR:                      
                            string2 = '/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/' + DIR + '/' +date + 'polymer_strings_model_' + model + base +ratio + '_delta_g_pol_' + L_delta_g_pol_label[i]+ "_repulsion" + repulsion +'_allsteps.csv'
                            file_exists2 = os.path.exists(string2)
                            print(string2)
                            if file_exists2 == True:
                                data2 = pd.read_csv(string2,on_bad_lines='skip')
                                prob_mean, prob_std = probability_base_base(data2)
                                matrix_prob= matrix_purpy_values(prob_mean,bases_pur1py1pur2py2[j])
                                matrix_sem = matrix_purpy_values(prob_std,bases_pur1py1pur2py2[j])
                                print(matrix_prob,matrix_sem)
                                # print(matrix_purpy_values(prob_std,bases_pur1py1pur2py2[j]))
                                repulsion_int = float(repulsion)
                                print(bases_pur1py1pur2py2[j], "repulsion:", repulsion_int, "G_pol:",delta_g_pol)
                                # print(matrix_purpy(bases_pur1py1pur2py2[j]))
                                OH_matrix,matrix = make_E_OH_matrix(bases_pur1py1pur2py2[j], repulsion_int)
                                print(OH_matrix)
                            else:
                                print("polymer string file does not exist")
                            

                    #put base prob versus bond strength graph per delta g_pol and rep strength
                    if repulsion =="-1" and delta_g_pol==-np.log(4):
                        break
                    else:
                        graph_prob_vs_OHstrength(bases_pur1py1pur2py2[j],prob_mean, prob_std, OH_matrix, L_delta_g_pol_label[i], repulsion,L_color_bases[j])
                    # plt.savefig('/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/' + DIR + '/prob_OHstrength_rep' +repulsion + 'Gpol' +L_delta_g_pol_label[i] +'_'+ base_combinations[j] )
                plt.show()
   


