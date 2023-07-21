"""
input: the compiled code you want to run=code_name, file with the parameters=parameter_file"""

import subprocess
#lines will be the list with the string defining what your are going to run
lines = ["./Documents/programming/DNA_replication_muriel/4_bases"] #name of compiled code you want to run, lines = ["./{code_name}"]


#read parameter file, with open('{parameter_file}' as f:)
with open('Documents/programming/DNA_replication_muriel/parameter_files/4bases/4_bases_ATiGiC_rep-1_equi') as f: 

    for i,line in enumerate(f): 
        if line.rstrip() == "":
            #read first part of the file with the parameters
            break
        if line.rstrip().split('=',1)[0] == "--L_g_pol" or line.rstrip().split('=',1)[0] == "--L_g_tt" or line.rstrip().split('=',1)[0] == "--ratio_template":
            #make {} of the parameter variabels defining the list of G_pol and G_tt
            linee= line.rstrip().split('=',1)[1].split('.',1)[0].replace(' ', ',')[:-1]
            lines.append(line.rstrip().split('=',1)[0] + "={" + linee + "}")
        else: 
            lines.append(line.rstrip().split('.',1)[0])
            #add the line you need for defining the parameters you want to run
        
        
lines = [' '.join(lines)]
subprocess.run(lines, shell=True) #run compiled program with parameters





