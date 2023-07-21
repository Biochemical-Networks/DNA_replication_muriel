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

# def graph_equilibration_check(length_polymer, error_prob,delta_g_bb, delta_g_tt, model):
#     plt.plot(length_polymer, error_prob, label = "$\Delta G_{bb}$ = %s, $\Delta G_{tt}$ = %s"%(delta_g_bb, delta_g_tt))
#     plt.title("Energetic: model: %s, Check when error is equilibrated as function of the length of the mRNA grown for different $\Delta G_{bb}$ and $\Delta G_{tt}$" %model)
#     plt.xlabel("length mRNA")
#     plt.ylabel("$\epsilon$")
#     plt.legend(bbox_to_anchor=(1.04, 1))
#     plt.subplots_adjust(right=0.7)
#     return

def graph_error_g_bb(L_delta_g_bb, L_epsilon, model, delta_g_tt, model_title, L_std_high, L_std_low,gtt_theory, color_def,fontsizee,widthh,ax):
    
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0+box.height*0.05, box.width*0.95 , box.height*0.95])
    ax.plot(L_delta_g_bb, L_epsilon, label = "$\Delta G_r - \Delta G_w$ = %s" %(delta_g_tt), color= color_def, linewidth=1*line_widths)
    if model_title[-1] =="1":
        plt.scatter([10], [gtt_theory], color= color_def)
    ax.fill_between(L_delta_g_bb, L_std_low, L_std_high, color= color_def, alpha=0.1)
    ax.set_title("%s"%(model_title), fontsize =fontsizee)
    ax.set_xlabel("Backbone strength $\Delta G_{bb}$",fontsize =fontsizee)
    ax.set_ylabel("Error probability $\epsilon$",fontsize =fontsizee)
    ax.legend(loc='lower center', bbox_to_anchor=(1.18, 0.6), ncol=1, fontsize =fontsizee)
    ax.set_ylim([0,0.6])
    ax.spines['left'].set_linewidth(widthh*0.5)
    ax.spines['right'].set_linewidth(widthh*0.5)
    ax.spines['top'].set_linewidth(widthh*0.5)
    ax.spines['bottom'].set_linewidth(widthh*0.5)
    ax.xaxis.set_tick_params(labelsize=fontsizee)
    ax.yaxis.set_tick_params(labelsize=fontsizee)
    return



def graph_error_g_bb_multiple_models_per_model(L_delta_g_bb, L_epsilon, model, delta_g_tt, model_title, L_std_high, L_std_low, color_def, line_style,fontsizee,widthh):
    
    plt.plot(L_delta_g_bb, L_epsilon, label = "$\Delta G_{tt}$ = %s, %s" %(delta_g_tt, model_title[4:12]), linestyle=line_style, color= color_def)
    plt.fill_between(L_delta_g_bb, L_std_low, L_std_high, color= color_def, alpha=0.1)
    plt.title("%s" %(model_title[0:4]+model_title[12:]))
    plt.xlabel("Backbone strength $\Delta G_{pol}$")
    plt.ylabel("Error probability $\epsilon$")
    plt.legend(bbox_to_anchor=(1.04, 1),fontsize =fontsizee)
    plt.subplots_adjust(right=0.7)
    ax.set_position([box.x0, box.y0+box.height*0.05, box.width*0.95 , box.height*0.95])
    return

def graph_error_g_bb_multiple_models_per_model2(L_delta_g_bb, L_epsilon, model, delta_g_tt, model_title, L_std_high, L_std_low, color_def, line_style,fontsizee,widthh,ax,f,k):
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0+box.height*0.05, box.width*0.95 , box.height*0.95])
    if line_style == 'solid':
        ax.plot(L_delta_g_bb, L_epsilon, label = "$\Delta G_r - \Delta G_w$ = %s" %(delta_g_tt), linestyle=line_style, color= color_def, linewidth=1*line_widths)
    else: 
        ax.plot(L_delta_g_bb, L_epsilon, linestyle=line_style, color= color_def, linewidth=1*line_widths)
    ax.fill_between(L_delta_g_bb, L_std_low, L_std_high, color= color_def, alpha=0.1)
    ax.set_title("%s" %(model_title[27:]),fontsize =fontsizee)
    ax.set_xlabel("Backbone strength $\Delta G_{bb}$",fontsize =fontsizee)
    ax.set_ylabel("Error probability $\epsilon$",fontsize =fontsizee)
    ax.legend(loc='lower center', bbox_to_anchor=(1.18, 0.6), ncol=1, fontsize =fontsizee)
    ax.set_ylim([0,0.6])
    ax.spines['left'].set_linewidth(widthh*0.5)
    ax.spines['right'].set_linewidth(widthh*0.5)
    ax.spines['top'].set_linewidth(widthh*0.5)
    ax.spines['bottom'].set_linewidth(widthh*0.5)
    ax.xaxis.set_tick_params(labelsize=fontsizee)
    ax.yaxis.set_tick_params(labelsize=fontsizee)
    return

def graph_error_g_bb_multiple_models_per_model3(L_delta_g_bb, L_epsilon, model, delta_g_tt, model_title, L_std_high, L_std_low, color_def, line_style,fontsizee,widthh):
    plt.plot(L_delta_g_bb, L_epsilon, label = "$\Delta G_{tt}$ = %s, %s" %(delta_g_tt, model_title), linestyle=line_style, color= color_def)
    plt.fill_between(L_delta_g_bb, L_std_low, L_std_high, color= color_def, alpha=0.1)
    # plt.title("Comparison models")
    plt.xlabel("Backbone strength $\Delta G_{pol}$")
    plt.ylim(0,0.6)
    plt.ylabel("Error probability $\epsilon$")
    plt.legend(bbox_to_anchor=(1.3, 1),fontsize =fontsizee)
    plt.subplots_adjust(right=0.7)
    return

def graph_error_g_bb_multiple_models_per_model4(L_delta_g_bb, L_epsilon, model, delta_g_tt, model_title, L_std_high, L_std_low, color_def, line_style,fontsizee,widthh,ax,f,k):
    ax = plt.gca()
    box = ax.get_position()
    labels = ["Fine grained BF x=1", "Fine grained BF x=$ e^{\Delta G_r - \Delta G_w + 5}$", "Coarse grained"]
    ax.set_position([box.x0, box.y0+box.height*0.05, box.width*0.95 , box.height*0.95])
    ax.plot(L_delta_g_bb, L_epsilon, label = labels[k], linestyle=line_style, color= color_def, linewidth=1*line_widths)
    ax.fill_between(L_delta_g_bb, L_std_low, L_std_high, color= color_def, alpha=0.1)
    ax.set_title("$\Delta G_r - \Delta G_w$ = %s" %delta_g_tt,fontsize =fontsizee)
    ax.set_xlabel("Backbone strength $\Delta G_{bb}$",fontsize =fontsizee)
    ax.set_ylabel("Error probability $\epsilon$",fontsize =fontsizee)
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3, fontsize =fontsizee)
    ax.set_ylim([0,0.6])
    ax.spines['left'].set_linewidth(widthh*0.5)
    ax.spines['right'].set_linewidth(widthh*0.5)
    ax.spines['top'].set_linewidth(widthh*0.5)
    ax.spines['bottom'].set_linewidth(widthh*0.5)
    ax.xaxis.set_tick_params(labelsize=fontsizee)
    ax.yaxis.set_tick_params(labelsize=fontsizee)
    # ax.set_xticklabels(fontsize =fontsizee)
    # ax.subplots_adjust(right=0.7)
    # f.tight_layout()
    return

def make_model_names(x_def, M_model, models, models_lose, models_titles):
    L_models = []
    L_models_titles = []
    for x in x_def:
        if x == "x1":
            x_title = "1"
        else:
            x_title = "$e^{\Delta G_r - \Delta G_w +5}$"
        for M in (M_model):
            if M =="MB1":
                M_title = "model BB"
            elif M == "MB2":
                M_title = "model BF"
            elif M =="MF1":
                M_title = "model FB"
            else: 
                M_title = "model FF"
            for i,model in enumerate(models):
                model = model + M + "_" + x
                L_models.append(model)
                L_models_titles.append(models_titles[i] +" " + M_title+ " x=" + x_title)
    for mod in models_lose:
        L_models.append(mod)
        L_models_titles.append(models_titles[-1])
    return L_models, L_models_titles

def graph_error_g_bb_inf(L_delta_g_bb, L_epsilon, model, delta_g_tt, model_title, L_std_high, L_std_low,gtt_theory, color_def,fontsizee,widthh,ax,ax2,theory):
    L_delta_g_bb.append(12)
    L_epsilon.append(theory)
    # ax = plt.gca()
    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0+box.height*0.05, box.width*0.95 , box.height*0.95])
    ax.plot(L_delta_g_bb[:-1], L_epsilon[:-1], label = "$\Delta G_r - \Delta G_w$ = %s" %(delta_g_tt), color= color_def, linewidth=1*line_widths)
    ax2.scatter(L_delta_g_bb, L_epsilon,  label = "$\Delta G_r - \Delta G_w$ = %s" %(delta_g_tt), color= color_def, linewidth=1*line_widths)
    if model_title[-1] =="1":
        plt.scatter([10], [gtt_theory], color= color_def)
    ax.fill_between(L_delta_g_bb[:-1], L_std_low, L_std_high, color= color_def, alpha=0.1)
    ax.set_title("%s"%(model_title), fontsize =fontsizee)
    ax.set_xlabel("Backbone strength $\Delta G_{bb}$",fontsize =fontsizee)
    ax.set_ylabel("Error probability $\epsilon$",fontsize =fontsizee)
    
    # ax.legend(loc='lower center', bbox_to_anchor=(1.18, 0.6), ncol=1, fontsize =fontsizee)
    ax.set_ylim([0,0.6])
    ax.set_xlim([-1,11])
    ax2.set_xlim([11,13])
    ax.yaxis.tick_left()
    ax2.yaxis.tick_right()
    ax2.tick_params(labelright='off')
    ax.tick_params(axis='y', labelsize=fontsizee, width=line_widths*0.5)
    ax2.tick_params(axis='y', labelsize=fontsizee,width=line_widths*0.5)
    ax.tick_params(axis='x', labelsize=fontsizee, width=line_widths*0.5)
    ax2.tick_params(axis='x', labelsize=fontsizee,width=line_widths*0.5)
    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.spines['left'].set_linewidth(line_widths*0.5)
    ax2.spines['right'].set_linewidth(line_widths*0.5)
    ax.spines['bottom'].set_linewidth(line_widths*0.5)
    ax2.spines['bottom'].set_linewidth(line_widths*0.5)
    ax.spines['top'].set_linewidth(line_widths*0.5)
    ax2.spines['top'].set_linewidth(line_widths*0.5)
    ax.spines['left'].set_linewidth(widthh*0.5)
    ax.spines['right'].set_linewidth(widthh*0.5)
    ax.spines['top'].set_linewidth(widthh*0.5)
    ax.spines['bottom'].set_linewidth(widthh*0.5)
    ax.xaxis.set_tick_params(labelsize=fontsizee)
    ax.yaxis.set_tick_params(labelsize=fontsizee)
    ax2.set_xticks([12])
    ax2.set_xticklabels(["$\infty$"], fontsize =fontsizee)
    # ax2.set_ylabel("Error probability $\epsilon$",fontsize =fontsizee)
    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1-d, 1+d), (-d, +d), **kwargs)
    ax.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    # ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
    # ax2.plot((-d, +d), (-d, +d), **kwargs)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0+box.height*0.1, box.width, box.height*0.9])
    ax2.set_position([box.x0+box.width +2*d, box.y0+box.height*0.1, box.width*0.2, box.height*0.9])
    ax.legend(loc='lower center', bbox_to_anchor=(0.62, -0.33), ncol=4, fontsize =fontsizee)
    return



L_delta_g_pol_label = ["-ln2", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
L_delta_g_pol_number = [-np.log(2), 0,1,2,3,4,5,6,7,8,9,10]
L_delta_g_tt = [0,2,4,6]

models = ["muriel_big_abs_", "no_abs_matrix_"]
models_title = ["Fine grained", "Fine grained computational", "Coarse grained model"]
all_steps = True
L_DIR = ["230418_equi_output", "230418output"]#, "230714output"]
x_def = ["xexp5", "x1"] #should be xdelta_g_tt_pow2
M_model = ["MB1","MB2", "MF1", "MF2"]
# number_steps = 0
# MX_models = ["x=$e^{\Delta G_{TT} +5}$, Model B","x=$e^{\Delta G_{TT} +5}$, Model F", "x=1, Model B", "x=1, Model F" ]

L_models, L_models_titles = make_model_names(x_def, M_model, models, ["Jenny_full_steps"], models_title)
date = date_now(datetime.datetime.now())
date_options = ["2023421", "2023422", "2023423", "2023424","202351","2023210","2023211",date]
L_length = ["300"]



#define colors
colormap = plt.cm.plasma #nipy_spectral, Set1,Paired   
color_list = [colormap(i) for i in np.linspace(0, 1,100)]
fontsizee = 23
line_widths = 3

print(L_models)

# i =10
# j= 8
# k=12
#graph multiple 1
# L_models = [L_models[i], L_models[i+3], L_models[i+6], L_models[i+9]] #per model all different M and x
# L_models_titles = [L_models_titles[i], L_models_titles[i+3], L_models_titles[i+6], L_models_titles[i+9]]
# L_models = [L_models[i], L_models[i+2], L_models[i+4], L_models[i+4]] #per model all different M and x = exp (i=0)
# L_models_titles = [L_models_titles[i], L_models_titles[i+2], L_models_titles[i+4], L_models_titles[i+6]] #per model all different M and x = exp
# L_models = [L_models[i], L_models[i+2], L_models[i+4], L_models[i+4]]#per model all different M and x = 1 (i=8)
# L_models_titles =  [L_models_titles[i], L_models_titles[i+2], L_models_titles[i+4], L_models_titles[i+6]]#per model all different M and x = 1

# graph multiple 2
# L_models = [L_models[i], L_models[i+1]] #FGA FGC for one x and for one M 
# L_models_titles = [L_models_titles[i], L_models_titles[i+1]] #(i=0,2,4 6, 8,10,12) , BF:2,10 


#graph multiple 3
# L_models = [L_models[i], L_models[i+2], L_models[j], L_models[j+2], L_models[k], L_models[k+2]]#, L_models[16]] #FGA for one x and for two M and CG
# L_models_titles = [L_models_titles[i], L_models_titles[i+2], L_models_titles[j], L_models_titles[j+2],L_models_titles[k], L_models_titles[k+2]]#,  L_models_titles[16]] ##FGA for one x and for two M and CG (i=0,4,8,12)
# L_delta_g_tt = [2,4]

# graph multiple 4
L_models = [L_models[2], L_models[10], L_models[16]] #FGA BF x =1 and x=exp
L_models_titles = [L_models_titles[2], L_models_titles[10], L_models_titles[16]]#
L_delta_g_tt = [4]


# line_styles = ['dashed', 'solid']

# L_models = L_models[i:i+3] #for same M and X diffent models
line_styles = ['solid', 'dashdot', 'dashed', 'dotted', (5, (10, 3)), (0, (3, 5, 1, 5))]
print(L_models)

L_gtt_theorie_xexp_B = [0.5000000000000001, 0.4985490335641511, 0.4983520193927259, 0.4983253445301224, 0.4983217342620933, 0.49832124566145497]
L_gtt_theorie_F = [0.5, 0.11920292202211756, 0.017986209962091555, 0.0024726231566347743, 0.00033535013046647816,4.5397868702434395e-05]
L_gtt_theorie_B = [0.5, 0.19251027051493744, 0.034723337448322934, 0.004920911195322951, 0.0006702507235977487, 9.078749428737782e-05]
L_gtt_theorie = L_gtt_theorie_xexp_B

f, ax = plt.subplots()
# f, (ax, ax2) = plt.subplots(1, 2, sharey=True, facecolor='w', gridspec_kw={'width_ratios': (9,1)})

for k, model in enumerate(L_models):


    for j, delta_g_tt in enumerate(L_delta_g_tt):
        L_epsilon_mean = []
        L_epsilon_std_high = []
        L_epsilon_std_low = []
        L_delta_g_bb_graph = []
        for i,delta_g_pol in enumerate(L_delta_g_pol_number):     
            for length in L_length:
                for DIR in L_DIR:
                    for date in date_options:
                        # print('/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/%s/%smodel_%s_delta_g_pol_%s_delta_g_tt_%s_allsteps.csv'%(DIR, date, model, L_delta_g_pol_label[i], delta_g_tt))                                                
                        file_exists = os.path.exists('/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/%s/%smodel_%s_delta_g_pol_%s_delta_g_tt_%slength_copy%s_allsteps.csv'%(DIR, date, model, L_delta_g_pol_label[i], delta_g_tt,length))
                        # print('/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/%s/%smodel_%s_delta_g_pol_%s_delta_g_tt_%slength_copy%s_allsteps.csv'%(DIR, date, model, L_delta_g_pol_label[i], delta_g_tt,length))
                        if file_exists == True:
                            print(DIR, date)
                            print(model)
                            data = pd.read_csv('/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/%s/%smodel_%s_delta_g_pol_%s_delta_g_tt_%slength_copy%s_allsteps.csv'%(DIR, date, model, L_delta_g_pol_label[i], delta_g_tt,length))
                            epsilon_mean = data['error probability'].mean()
                            epsilon_std = data['error probability'].std()
                            print(model, length)
                            print(delta_g_pol,delta_g_tt)
                            print(epsilon_mean, epsilon_std)#/(np.sqrt(len(data['error probability']))))
                            L_epsilon_mean.append(epsilon_mean)
                            L_epsilon_std_high.append(epsilon_mean + epsilon_std)
                            L_epsilon_std_low.append(epsilon_mean - epsilon_std)
                            L_delta_g_bb_graph.append(delta_g_pol)
                            
                        else: 
                            file_exists = os.path.exists('/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/%s/%smodel_%s_delta_g_pol_%s_delta_g_tt_%s_allsteps.csv'%(DIR, date, model, L_delta_g_pol_label[i], delta_g_tt))
                            # print('/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/%s/%smodel_%s_delta_g_pol_%s_delta_g_tt_%slength_copy%s_allsteps.csv'%(DIR, date, model, L_delta_g_pol_label[i], delta_g_tt,length))
                            if file_exists == True:
                                print(DIR, date)
                                print(model)
                                data = pd.read_csv('/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/%s/%smodel_%s_delta_g_pol_%s_delta_g_tt_%s_allsteps.csv'%(DIR, date, model, L_delta_g_pol_label[i], delta_g_tt))
                                epsilon_mean = data['error probability'].mean()
                                epsilon_std = data['error probability'].std()
                                print(model, length)
                                print(delta_g_pol,delta_g_tt)
                                print(epsilon_mean, epsilon_std)#/(np.sqrt(len(data['error probability']))))
                                L_epsilon_mean.append(epsilon_mean)
                                L_epsilon_std_high.append(epsilon_mean + epsilon_std)
                                L_epsilon_std_low.append(epsilon_mean - epsilon_std)
                                L_delta_g_bb_graph.append(delta_g_pol)


                        # else:
                        #     print("file does not exist")
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
        # graph_error_g_bb(L_delta_g_bb_graph, L_epsilon_mean, model, delta_g_tt, L_models_titles[k], L_epsilon_std_high, L_epsilon_std_low, L_gtt_theorie[j], color_list[j*20], fontsizee,line_widths,ax)
        # graph_error_g_bb_multiple_models_per_model(L_delta_g_bb_graph, L_epsilon_mean, model, delta_g_tt, L_models_titles[k], L_epsilon_std_high, L_epsilon_std_low, color_list[j*10], line_styles[k], fontsizee,line_widths)
        # graph_error_g_bb_multiple_models_per_model2(L_delta_g_bb_graph, L_epsilon_mean, model, delta_g_tt, L_models_titles[k], L_epsilon_std_high, L_epsilon_std_low, color_list[j*20], line_styles[k], fontsizee,line_widths,ax,f,k)
        # graph_error_g_bb_multiple_models_per_model3(L_delta_g_bb_graph, L_epsilon_mean, model, delta_g_tt, L_models_titles[k], L_epsilon_std_high, L_epsilon_std_low, color_list[j*40], line_styles[k], fontsizee,line_widths)
        graph_error_g_bb_multiple_models_per_model4(L_delta_g_bb_graph, L_epsilon_mean, model, delta_g_tt, L_models_titles[k], L_epsilon_std_high, L_epsilon_std_low, color_list[40], line_styles[k], fontsizee,line_widths,ax,f,k)
        # graph_error_g_bb_inf(L_delta_g_bb_graph, L_epsilon_mean, model, delta_g_tt, L_models_titles[k], L_epsilon_std_high, L_epsilon_std_low, L_gtt_theorie[j], color_list[j*20], fontsizee,line_widths,ax,ax2,L_gtt_theorie[j])

    # print(model, delta_g_tt, L_epsilon_mean)
plt.show()


