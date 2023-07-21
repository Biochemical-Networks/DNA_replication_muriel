# DNA_replication_muriel

### 230412_model_big_abs_matrix_differentX_grow_mRNA_wr_energetic_beginww.cpp: 

Code for making chapter 2 results of my thesis making the fine grained systems using an absorption matrix

2 monomer system of right(0) and wrong(1) 

To run on the cluster :

1) g++ -lboost_program_options 230412_model_big_abs_matrix_differentX_grow_mRNA_wr_energetic_beginww.cpp -o Muriel
2) change 230131_read_parameters_models.py 
3) change run_model.sh
4) qsub run_model.sh


Paramaters:

    `bool all_steps = false`;    -> can be removed I think   
    
    `bool multiple_DNA = true`;   
    
            -> if false: Makes files for every polymer made describing all transitions, e.g did the system add a monomer (0 or 1) or did the ystem move back removing the current monomer tip.  
            
            -> if true: makes one file of every polymer made describing only the end polymer and the error probability of that polymer.
            
    `string DIR = "230714output"`;- -> directory where it puts the files
    
    `string x_def = "1"`; --> x definition of the polymerisation step
    
    `string M_model = "MB2"`; -> fine grained model used (MB2 = BF, MB1 = BB, MF2 = FF, MF1 = FB)
    
    `int number_steps = 30`; -> number of polymers simulated
    
    `int length_mRNA = 300`; -> length of ervey polymer that is siumlated
    
    `bool variables_in_main = false`; 
    
            -> if false: run script on your computer
            
            -> if true run the script on cluster Then aslo comment out the upper section of the code and comment in the section where i define parameters using external files (needed for cluster)
            
    `vector<double> L_delta_g_tt = {3,5}`;   --> delta g_tt = Delta G_r - Delta G_w: parameter defining copy/template interactions
    
    `vector<double> L_delta_g_pol = {0,1,2,3,4,5,6,7,8,9,10}`;   --> Delta G_bb: backbone strength variation


### 230412_grow_mRNA_wr_energetic_without_absorption_matrix_beginw.cpp

Code for makring chapter 2 of my thesis making the fully computational model

2 monomer system: right(0) wrong(1)

1) g++ -lboost_program_options 230412_grow_mRNA_wr_energetic_without_absorption_matrix_beginw.cpp -o no_abs
2) change 230131_read_parameters_models.py 
3) change run_model.sh
4) qsub run_model.sh

Paramaters te same as for 230412_model_big_abs_matrix_differentX_grow_mRNA_wr_energetic_beginww.cpp


### 230526_Big_abs_BF_x1_4differentbases_with_pupu_pypy_coststep13.cpp
Code for results chapter 3 + 4

4 nucleotide system of polymer copying

To run on the cluster :

1) g++ -lboost_program_options 230526_Big_abs_BF_x1_4differentbases_with_pupu_pypy_coststep13.cpp -o 4_bases
2) change 230413_read_parameters_4bases_forloop.py
3) change run_model_4bases.sh
4) qsub run_model_4bases.sh


Paramaters:

    `bool variables_in_main =true`;  

                -> if false: run script on your computer
            
            -> if true run the script on cluster Then aslo comment out the upper section of the code and comment in the section where i define parameters using external files (needed for cluster)
             
    `string DIR = "";` = ("DIR", po::value<string>(), "output directory name") --> directory for output files
    
    `bool multiple_DNA`; = ("multiple_DNA", po::value<bool>(), "compute multiple DNA string") 

    
            -> if false: Makes files for every polymer made describing all transitions, e.g did the system add a monomer (0 or 1) or did the ystem move back removing the current monomer tip.  
            
            -> if true: makes one file of every polymer made describing only the end polymer and the error probability of that polymer.
    
    `int number_steps`;  = ("number_DNA", po::value<int>(), "number of DNA computed") --> number of polymer string made
    
    `string x_def = ""`;  =   ("x_def", po::value<string>(), "strength x, backbone")  --> x definition of polymerisation step, b+ and b- 
    
    `int length_mRNA`; = ("mRNA_length", po::value<int>(), "length of the mRNA string formed") -> length of polymers computed 
    
    `vector<double> L_delta_g_pol`; = ("L_g_pol", po::value< vector<double> >(), "delta G_pol definitions list")  -> Delta G_bb: backbone strength variation
    
    `vector<double> ratio_bases_template`; =  ("ratio_template", po::value< vector<double> >(), "the ratio of the four bases we are using") --> ratio of purine1,purine2, pyrimidn1 and pyrimdine2 on template
    
    `string purine_1` = ""; = ("purine_1", po::value<string>(), "purine 1 base definition") --> which purine used 
    
    `string purine_2` = ""; =  ("purine_2", po::value<string>(), "purine 2 base definition") -->which purine used
    
    `string pyrimidines_1` = ""; =  ("pyrimidines_1", po::value<string>(), "pyrimidines 1 base definition") --> which pyrimidine used, complementary to purine1
    
    `string pyrimidines_2` = ""; =  ("pyrimidines_2", po::value<string>(), "pyrimidines 2 base definition") --> which pyrimidine used, complementary to purine2
    
    `int repulsion_00_11`; = ("repulsion_11_00", po::value<int>(), "strength repulsion purine-pruine or pyrimidine-pyrimidine binding") --> \Delta G_{cost}

### 230501_Big_abs_BF_x1_4differentbases_with_pupu_pypy_editcost

4 nucleotide system of polymer copying putting the cost on b+ 


### 230630_limit_G_bb_abs_matrix_6_bases.py

Code for theoretical result 6 nucleotide alphabet: conclusion and outlook chapter of thesis

Note that all fine grained models of 6 nucleotides do not work --> they do not give the right results yes


### run_model.sh
For running or 2 monomer system codes

### run_models_4bases.sh
For running all codes with 4 or 6 bases 

### mathematica/230714_absoption_matrix_example.nb
absoption matrix of the simplified system of figure 8

### mathematic/230317_absoption_matrix_BFx1_4bases.nb
Big absoption matrix(figure 21, thesis) for 4 bases system

### mathematica/230327_absoption_matrix_BFx1_6bases.nb
Big absoption matrix(figure 21, thesis) but then for 6 bases instead of 4

### mathematica/230526_absoption_matrix_4_bases_g_bb_limit.nb
absoption matrix of system of 4 bases in the limit of inifite backbone strength (figure 12 of thesis) with x=1

### mathematica/230630_absoption_matrix_6_bases_g_bb_limit.nb
absoption matrix of system of 6 bases in the limit of inifite backbone strength with x=1

### mathematica/230714_absoption_matrix_xexp_simplified_system.nv               
Absorption matrix of system chapter 2 in the limit of infinite backbone strength((figure 12 of thesis) with x can be varied

### all other codes are old or not working correctly!              
                
               
                
               
               
               
                
    

