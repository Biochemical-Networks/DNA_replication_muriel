"""
input: the compiled code you want to run=code_name, file with the parameters=parameter_file"""

import subprocess
#lines will be the list with the string defining what your are going to run
lines = ["./Big_abs"] #name of compiled code you want to run, lines = ["./{code_name}"]

#read parameter file, with open('{parameter_file}' as f:)
with open('parameter_files/Big_abs_M1_x1_parameters') as f: 
    for i,line in enumerate(f): 
        if line.rstrip() == "":
            #read first part of the file with the parameters
            break
        if line.rstrip().split('=',1)[0] == "--L_g_pol" or line.rstrip().split('=',1)[0] == "--L_g_tt":
            #make {} of the parameter variabels defining the list of G_pol and G_tt
            linee= line.rstrip().split('=',1)[1].split('.',1)[0].replace(' ', ',')[:-1]
            lines.append(line.rstrip().split('=',1)[0] + "={" + linee + "}")
        else: 
            lines.append(line.rstrip().split('.',1)[0])
            #add the line you need for defining the parameters you want to run
        
        
lines = [' '.join(lines)]
subprocess.run(lines, shell=True) #run compiled program with parameters





# --multiple_DNA=1.
# --DIR=test_output.
# --number_DNA=3.
# --x_def=n.
# --M_model=n.
# --mRNA_length=30.
# --L_g_pol=1 2 .
# --L_g_tt=1 3 .

# we are computing multiple DNA strings: 1. #1=true, #0=false
# output file is set to: test_output.
# number of DNA computed: 3.
# x value: n.
# M model: n.
# length of the mRNA formed: 30.
# delta G_pol list definition: 1 2 .
# delta G_tt list definition: 1 2 .