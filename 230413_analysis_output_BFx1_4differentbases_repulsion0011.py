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

def graph_error_g_bb_3tables(L_delta_g_bb, L_epsilon, model_title, bases, L_std_high, L_std_low, color_def, repulsion, matrix_table1, colnames,iteration,graph,table1, table2, table3):
    # matrix_table1 = matrix_table1.round(decimals=2)
    if iteration ==0:
        table1.axis('tight')
        table1.axis('off')
        the_table1 = table1.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        the_table1.auto_set_font_size(False)
        the_table1.set_fontsize(9)
        table2.plot()
        table3.plot()
    if iteration ==1:
        table2.axis('tight')
        table2.axis('off')
        table2.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        # table2.set_fontsize(14)
        table1.plot()
        table3.plot()
    if iteration ==2:
        table3.axis('tight')
        table3.axis('off')
        table3.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        # table3.set_fontsize(14)
        table1.plot()
        table2.plot()
        
    graph.plot(L_delta_g_bb, L_epsilon, label = "%s" %(bases), color= color_def)
    graph.fill_between(L_delta_g_bb, L_std_low, L_std_high, color= color_def, alpha=0.1)
    graph.set_ylim([0,0.6])
    graph.set_xlabel("Backbone strength $\Delta G_{pol}$")
    graph.set_ylabel("Error probability $\epsilon$")
    box = graph.get_position()
    graph.set_position([box.x0, box.y0+box.height*0.1, box.width, box.height*0.9])
    graph.legend(loc='lower center', bbox_to_anchor=(0.5, -0.2), ncol=3)

    if repulsion =="-1":
        graph.set_title("%s: Effect $\Delta G_{pol}$ on error probability of the copy polymer, not allowing pur-pur or py-py binding"%(model_title))
    else:
        graph.set_title("%s: Effect $\Delta G_{pol}$ on error probability of the copy polymer, with repulsion %s"%(model_title,repulsion))
    
    return

def graph_error_g_bb_6tables(L_delta_g_bb, L_epsilon, model_title, bases, L_std_high, L_std_low, color_def, repulsion, matrix_table1, colnames,iteration,graph,table1, table2, table3, table4,table5,table6,OH_matrix):
    matrix_table1 = matrix_table1.round(decimals=2)
    if iteration ==0:
        table1.axis('tight')
        table1.axis('off')
        the_table1 = table1.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        # the_table1.auto_set_font_size(False)
        # the_table1.set_fontsize(9)
        table4.axis('tight')
        table4.axis('off')  
        the_table4 = table4.table(cellText=OH_matrix,colLabels=colnames,rowLabels = colnames, loc='center')
        # the_table4.auto_set_font_size(False)
        # the_table4.set_fontsize(9)
        table2.plot()
        table3.plot()
    if iteration ==1:
        table2.axis('tight')
        table2.axis('off')
        table2.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        # table2.set_fontsize(9)
        table5.axis('tight')
        table5.axis('off')  
        the_table5 = table5.table(cellText=OH_matrix,colLabels=colnames,rowLabels = colnames, loc='center')
        # the_table5.auto_set_font_size(False)
        # the_table5.set_fontsize(9)
        table1.plot()
        table3.plot()
    if iteration ==2:
        table3.axis('tight')
        table3.axis('off')
        table3.table(cellText=matrix_table1,colLabels=colnames,rowLabels = colnames, loc='center')
        # table3.set_fontsize(9)
        table6.axis('tight')
        table6.axis('off')  
        the_table6 = table6.table(cellText=OH_matrix,colLabels=colnames,rowLabels = colnames, loc='center')
        # the_table6.auto_set_font_size(False)
        # the_table6.set_fontsize(9)
        table1.plot()
        table2.plot()

    

    graph.plot(L_delta_g_bb, L_epsilon, label = "%s" %(bases), color= color_def)
    graph.fill_between(L_delta_g_bb, L_std_low, L_std_high, color= color_def, alpha=0.1)
    graph.set_ylim([0,0.6])
    graph.set_xlabel("Backbone strength $\Delta G_{pol}$")
    graph.set_ylabel("Error probability $\epsilon$")
    box = graph.get_position()
    graph.set_position([box.x0, box.y0+box.height*0.1, box.width, box.height*0.9])
    graph.legend(loc='lower center', bbox_to_anchor=(0.5, -0.2), ncol=3)

    if repulsion =="-1":
        graph.set_title("%s: Effect $\Delta G_{pol}$ on error probability of the copy polymer, not allowing pur-pur or py-py binding"%(model_title))
    else:
        graph.set_title("%s: Effect $\Delta G_{pol}$ on error probability of the copy polymer, with repulsion %s"%(model_title,repulsion))
    
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

def change_base_string_to_int(base):
    base_int = 0
    if base == "aA":
        base_int =1
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
    for i in range(0,len(base_i)):
        if base_i[i] != base_j[i]:
            match += 1
    return -1*(2 * match - 3) +rep

def make_E_OH_matrix(bases, rep):
    
    OH_matrix=np.zeros((len(bases), len(bases)))
    strength = 0

    for i, basei in enumerate (bases):
        for j, basej in enumerate (bases):
            base_i_int = change_base_string_to_int(basei)
            base_i = base_OH_groups(base_i_int)
            base_j_int = change_base_string_to_int(basej)
            base_j = base_OH_groups(base_j_int)
            base_i_dev = base_i_int%2
            base_j_dev = base_j_int%2
            tot_dev = base_i_dev +base_j_dev 
            if tot_dev ==1:
                strength = OH_matches_strength(base_i,base_j,0)
            else: 
                strength = OH_matches_strength(base_i,base_j,rep)
                if rep ==-1:
                    strength = 0
            OH_matrix[i,j] = strength
    return OH_matrix

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

plot1 =True
plot2 =False

if plot1 ==True:
    L_delta_g_pol_label = ["-ln4","-ln2", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    L_delta_g_pol_number = [-np.log(4),-np.log(2), 0,1,2,3,4,5,6,7,8,9,10]

    all_steps = True
    DIR = "230414output"
    model = "muriel_big_abs_BF_x1"
    model_title = "BF x=1"
    bases_pur1py1pur2py2=[["iG", "iC", "X", "K"],["aA","T","G","C"],["aA","T","X","K"]]#["delta","beta","X","K"]]
    base_combinations = ["iGiCXK","aATGC", "aATXK"] #"GCXK","aATXK","aATGC", "deltabetaXK",
    ratio ='1111'
    # number_steps = 0
    # MX_models = ["x=$e^{\Delta G_{TT} +5}$, Model B","x=$e^{\Delta G_{TT} +5}$, Model F", "x=1, Model B", "x=1, Model F" ]

    # L_models, L_models_titles = make_model_names(base_combinations, model, model_title)exit()\e
    L_models = [model]
    L_models_titles = [model_title]
    date = date_now(datetime.datetime.now())
    date_options = [date]
    L_repulsion = ["-1","0"]#,"20"]#"0","1","2","3","4", "5", "6","20"]
    number_tables=3

    #define colors
    colormap = plt.cm.plasma #nipy_spectral, Set1,Paired   
    color_list = [colormap(i) for i in np.linspace(0, 1,100)]

    print(L_models)



    # line_styles = ['dashed', 'solid']

    # L_models = L_models[i:i+3] #for same M and X diffent models
    line_styles = ['solid', 'dashdot', 'dashed', 'dotted', (5, (10, 3)), (0, (3, 5, 1, 5))]
    print(L_models)



    for k, model in enumerate(L_models):
        for repulsion in L_repulsion:
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
            for j, base in enumerate(base_combinations):
                L_epsilon_mean = []
                L_epsilon_std_high = []
                L_epsilon_std_low = []
                L_delta_g_bb_graph = []
                for i,delta_g_pol in enumerate(L_delta_g_pol_number):     

                    for date in date_options:
                        string = '/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/' + DIR + '/' +date + 'model_' + model + base +ratio + '_delta_g_pol_' + L_delta_g_pol_label[i]+ "_repulsion" + repulsion +'_allsteps.csv'
                        file_exists = os.path.exists(string)
                        print(string)
                        if file_exists == True:
                            data = pd.read_csv(string)
                            epsilon_mean = data['error probability'].mean()
                            epsilon_sem = data['error probability'].std()#/(np.sqrt(len(data['error probability'])))
                            print(epsilon_mean, epsilon_sem)
                            L_epsilon_mean.append(epsilon_mean)
                            L_epsilon_std_high.append(epsilon_mean + epsilon_sem)
                            L_epsilon_std_low.append(epsilon_mean - epsilon_sem)
                            L_delta_g_bb_graph.append(delta_g_pol)
                            
                        else:
                            print("file does not exist")
                        
                        string2 = '/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/' + DIR + '/' +date + 'polymer_strings_model_' + model + base +ratio + '_delta_g_pol_' + L_delta_g_pol_label[i]+ "_repulsion" + repulsion +'_allsteps.csv'
                        file_exists2 = os.path.exists(string2)
                        print(string2)
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
                            OH_matrix = make_E_OH_matrix(bases_pur1py1pur2py2[j], repulsion_int)
                            print(OH_matrix)
                        # else:
                            print("polymer string file does not exist")
                        

                    

                #define error as all mismatches and plot the error as function of delta G_bb
                # if number_tables ==3:
                #     graph_error_g_bb_3tables(L_delta_g_bb_graph, L_epsilon_mean, L_models_titles[k], base, L_epsilon_std_high, L_epsilon_std_low, color_list[j*15],repulsion, OH_matrix, bases_pur1py1pur2py2[j],j,graph,table1,table2,table3)
                # elif number_tables ==6:
                #     graph_error_g_bb_6tables(L_delta_g_bb_graph, L_epsilon_mean, L_models_titles[k], base, L_epsilon_std_high, L_epsilon_std_low, color_list[j*15],repulsion, matrix_prob, bases_pur1py1pur2py2[j],j,graph,table1,table2,table3,table4,table5,table6,OH_matrix)
            # plt.show()


if plot2==True:
    L_delta_g_pol_label = ["-ln2", "0", "1", "2","3","4", "5","6", "7", "8" ,"9", "10"]
    L_delta_g_pol_number = [-np.log(2), 0,1,2,3,4,5,6,7,8,9,10]

    all_steps = True
    DIR = "230331output"
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
    date_options = ["2023331", "202344","202343", "202345","202346", date]
    L_repulsion = ["-1","0","1","2","3","4", "5", "6","20"]


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
                            OH_matrix = make_E_OH_matrix(bases_pur1py1pur2py2[j], repulsion_int)
                            print(OH_matrix)
                        else:
                            print("polymer string file does not exist")
                        

                    #put base prob versus bond strength graph per delta g_pol and rep strength
                    graph_prob_vs_OHstrength(bases_pur1py1pur2py2[j],prob_mean, prob_std, OH_matrix, L_delta_g_pol_label[i], repulsion,L_color_bases[j])
                    # plt.savefig('/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/' + DIR + '/prob_OHstrength_rep' +repulsion + 'Gpol' +L_delta_g_pol_label[i] +'_'+ base_combinations[j] )
                plt.show()
   


