/*
27-10-2022, Muriel Louman, AMOLF

grow mRNA from a random template string using 4 bases, aA, C, G, T. Assuming a energetic case, delta G_k = 0 always. 
delta_G is defined by the number of right and wrong matches. When there is a right match +1 (attraction), when there \
is a wrong match -1 (repulsion). 

model muriel:
Using a absorption matrix consisting of the probabilities to move between the absorbing states and the transition 
states (defined by the rates). Use of a probability to move forward to t1r or t1w and backward to t2r or t2w.  
*/

#include <iostream>
#include <cmath> // math functions
#include <vector> // vector container
#include <string> // text strings
#include <fstream> // output to file
#include <sstream> 
#include <random> // random number generator
#include <bits/stdc++.h> //for count function
using namespace std;

typedef struct {
    //define the part of the OH group: does it have a O- or a H at a surtain position.
    int position_1;
    int position_2;
    int position_3;
    } OH_groups ;

std::vector<int> base_OH_groups(int base){
    std::vector<int> OH_groups;
    if(base == 1){
        //base = 1 means, base is aA.
        OH_groups = {1, 0, 1};
    }
    else if(base == 2){
        //base = 1 means, base is aA.
        OH_groups = {0, 1, 0};
    }
    else if(base == 3){
        //base = 1 means, base is aA.
        OH_groups = {0, 1, 1};
    }
    else if(base == 4){
        //base = 1 means, base is aA.
        OH_groups = {1, 0, 0};
    }
    return OH_groups;
}


typedef struct {
    //typedef gives a way to define the type of struct somewhere without typing struct in front of it
    // define variables based on how good the match is between the OH groups
    double concentration;
    double concentration_star;
    double delta_G;
    double delta_G_bb;
    double delta_G_k;
    } variable_definition_t ;

int matches(std::vector<int> base_i, std::vector<int> base_j){
    int match = 0;
    for(int i=0; i< base_i.size(); ++i){
        if(base_i[i] != base_j[i]){
            match += 1;
        }
    }
    return match;
}

int strength_of_matches(int matches){
    return 2 * matches - 3;
}

// define needed free energies and concentrations dependent on template en mRNA binding
variable_definition_t def_variables_given_ij(int i,int j, double G_bb, double G){
    variable_definition_t variable_definition;
    // i and j can be 1,2,3,4
    std::vector<int> base_i = base_OH_groups(i);
    std::vector<int> base_j = base_OH_groups(j);
    int OH_matches = matches(base_i, base_j);

    variable_definition.concentration = 1;
    variable_definition.concentration_star= 1.5;
    variable_definition.delta_G = strength_of_matches(OH_matches) * G;
    variable_definition.delta_G_bb = G_bb;
    variable_definition.delta_G_k = 0;

    return variable_definition;
}

// define rates
double a_2(){
    return 1;
}

double a_1(double k_on, int i, int j, double G_bb, double G){
    variable_definition_t variables = def_variables_given_ij(i,j,G_bb, G);
    return exp(-1 * variables.delta_G);
}

double b_2(double k_bb, int j_M, int j_M1, double G_bb, double G){
    //j_M is of s[M] j_M1 is of s[M+1]
    variable_definition_t variables = def_variables_given_ij(j_M, j_M1, G_bb, G);
    return exp(variables.delta_G_bb);
}

double b_1(){
    return 1;
}

double c_2(){
    return 1;
}

double c_1(double k_on, int i, int j, double G_bb, double G){
    variable_definition_t variables = def_variables_given_ij(i,j,G_bb, G);
    return exp(variables.delta_G);
}



void print_matrix(std::vector<int> input){
    int n = sizeof(input)/sizeof(input[0]);
 
    // loop through the array elements
    for(int i=0; i< input.size(); ++i){
        std::cout << input.at(i) << ' ';
        }
}

//Create a structure variable called mRNA_move
struct structure_mRNA_move{
    std::vector<int> template_string;
    std::vector<int> polymer;
    int length;
    string transition_state;
    double error;
    };



void print_matrix_of_structures(std::vector<structure_mRNA_move> mRNA_movess){
    // loop through the array elements, which are all structures defined by structure_mRNA_move
    for(int i=0; i< mRNA_movess.size(); ++i){
        cout << "template DNA" << endl;
        print_matrix(mRNA_movess[i].template_string);
        cout << "creates new polymer" << endl;
        print_matrix(mRNA_movess[i].polymer);
        cout << "of length" << endl;
        cout << mRNA_movess[i].length << endl;
        cout << "transition states used" << endl;
        cout << mRNA_movess[i].transition_state << endl;
        cout << "   " << endl;
        cout << "with error probability"  << endl;
        cout << mRNA_movess[i].error << endl;
        }
}


// void save_matrix_of_structures(std::vector<structure_mRNA_move> mRNA_movess, std::string filename){
    // //takes to much memory
//     //saves the structures of different computations in a matrix
//     ofstream fss;
//     fss.open(filename.c_str());
//     // write the file headers
//     fss << "template string" << ","  << "start polymer" << "," << "length polymer" << "," << "transition state used" <<"," << "error probability"  << std::endl;
//     // loop through the array elements, which are all structures defined by structure_mRNA_move
//     for(int i=0; i< mRNA_movess.size(); ++i){
//         std::vector<int> in = mRNA_movess[i].template_string;
//         // loop through the array elements
//         for(int j=0; j< in.size(); ++j){
//             fss << in.at(j);
//             }
//         fss << ",";
//         std::vector<int> input = mRNA_movess[i].polymer;
 
//         // loop through the array elements
//         for(int j=0; j< input.size(); ++j){
//             fss << input.at(j);
//             }
//         fss << "," << mRNA_movess[i].length <<"," << mRNA_movess[i].transition_state << ","  << mRNA_movess[i].error << std ::endl;
//     }
//     fss.close();
// }




string date(time_t now){
    //define date in this way: 20221011, which is 11 October 2022
    tm *ltm = localtime(&now);
    string datum = std::to_string(1900 + ltm->tm_year) + std::to_string(1 + ltm->tm_mon) + std::to_string(ltm->tm_mday);
    return datum;
}

std::vector<int> template_string(int string_length, double aA_ratio, double T_ratio, double G_ratio, double C_ratio){
    // create a (uniformly) randomized template string 
    std::vector<int> string;
    std::default_random_engine generator;
    double ratios = aA_ratio + T_ratio + C_ratio + G_ratio;
    std::uniform_real_distribution<double> distribution(0,ratios);
    for(int j=0; j<string_length; j++){
        double random_number = distribution(generator);
        if(random_number < aA_ratio){
            string.push_back(1);
        }
        else if(random_number < aA_ratio + T_ratio){
            string.push_back(2);
        }
        else if(random_number < aA_ratio + T_ratio + G_ratio){
            string.push_back(3);
        }
        else{
            string.push_back(4);
        }
    }
    return string;
}

int error(int i, int j){
        //error vector with -5 is incorrect, 5 is correct
        int E = 0;
        if(i == 1 && j==2){
            E = 5;
            }
        else if(i == 2 && j==1){
            E = 5;
            }
        else if(i == 3 && j==4){
            E = 5;}
        else if(i == 4 && j==3){
            E = 5;}
        else{
            E = -5;}
    return E;
    }


int main(){
    // define variables
    double k_on = 1.0;
    double k_bb = 1.0;
    int length_mRNA = 2000;

    // random number generator for uniform distr between 0 and 1
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    double L_delta_g_bb[] = {0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10};  //{-1,-0.7,-0.6, -0.5,-0.4,-0.3,-0.2,-0.1};//we would like to have more but takes computer time
    double delta_g = {1}; //strength of a OH bond/repulsion
    // string L_delta_g_bb_string[] = {"min1", "min0.7", "min0.6", "min0.5", "min0.4", "min0.3", "min0.2", "min0.1"};

    // define template string with ratios 
    int rows = 5 ;
    int columns = 4;
    double L_proportions_template[rows][columns] ={{1,1,1,1}, {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}};


    for (int k =0; k<rows;k++){
        std::vector<int> mRNA_template = template_string(length_mRNA + 1, L_proportions_template[k][0], L_proportions_template[k][1], L_proportions_template[k][2], L_proportions_template[k][3]);
        
        int iteration_delta_g_bb = 0;
        //loop over different g_r = g_w and g_bb
        for (double delta_g_bb : L_delta_g_bb) {
            //error vector with -5 is incorrect, 5 is correct
            std::vector<int> mRNA_error(length_mRNA,0);

            // mRNA string we are going to grow, start with a 1(=aA) on the first place
            std::vector<int> mRNA_string;
            mRNA_string.push_back(1); 


            //define variables later used in code
            int i = 0;
            int transition_state = 0;
            int iteration = 0;
            // std::vector<structure_mRNA_move> mRNA_moves;
            structure_mRNA_move mRNA_moves;
            string transition_state_string = "0";

            // saving the information
            //save the polymer and the path it took to get there
            int rounded_delta_g_bb = round(delta_g_bb);
            int rounded_delta_g = round(delta_g);
            int rounded_ratio_aA = round(L_proportions_template[k][0]);
            int rounded_ratio_T = round(L_proportions_template[k][1]);
            int rounded_ratio_G = round(L_proportions_template[k][2]);
            int rounded_ratio_C = round(L_proportions_template[k][3]);
            string date_now = date(time(0));
            string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/221124output/";
            string filename_output = path + date_now + "model_muriel_delta_g_bb_" + std::to_string(rounded_delta_g_bb) + "strength_OH_bond" + std::to_string(rounded_delta_g) + "_aATGC_ratio" +  std::to_string(rounded_ratio_aA) + std::to_string(rounded_ratio_T)+ std::to_string(rounded_ratio_G) + std::to_string(rounded_ratio_C) + ".csv";
            fstream fss;
            fss.open(filename_output.c_str(), ios::out | ios::app);
            // write the file headers
            fss << "start polymer" << "," << "length polymer" << "," << "transition state used" <<"," << "error probability"  << "\n";
   
            
            while(mRNA_string.size()< length_mRNA){
            // for(int j =0;j<20;j++){
                // create a random number betwen 0 and 1
                double random_number1 = distribution(generator);

                // keep track of the current position. S(M)
                i = mRNA_string.size() -1;

                                                        
                // if string is shorter than 2 monomers, always go add 1(25%), 2(25%), 3(25%) or 4(25%) to the mRNA
                if(mRNA_string.size() == 1){              
                    if(random_number1 <0.25){
                        //add a 1 (=aA) to the mRNA
                        transition_state = -1;
                        transition_state_string = "t1_aA";
                    }
                    else if(random_number1 <0.5){
                        //add a 2(=T) to the mRNA
                        transition_state = -2;
                        transition_state_string = "t1_T";
                    }
                    else if(random_number1 <0.75){
                        //add a 3(=G) to the mRNA
                        transition_state = -3;
                        transition_state_string = "t1_G";
                    }
                    else{
                        //add a 4(=C) to the mRNA
                        transition_state = -4;
                        transition_state_string = "t1_C";
                    }
                }

                else{
                    // define probabilities needed for going to t1 or t2
                    double alpha_2_add_aA = a_2()/(a_2() + a_2() + a_2() + a_2() + c_1(k_on, mRNA_template[i-1], mRNA_string[i-1], delta_g_bb, delta_g));
                    double alpha_2_add_T  = a_2()/(a_2() + a_2() + a_2() + a_2() + c_1(k_on, mRNA_template[i-1], mRNA_string[i-1], delta_g_bb, delta_g));
                    double alpha_2_add_G  = a_2()/(a_2() + a_2() + a_2() + a_2() + c_1(k_on, mRNA_template[i-1], mRNA_string[i-1], delta_g_bb, delta_g));
                    double alpha_2_add_C  = a_2()/(a_2() + a_2() + a_2() + a_2() + c_1(k_on, mRNA_template[i-1], mRNA_string[i-1], delta_g_bb, delta_g));
                    double gamma_1 = c_1(k_on, mRNA_template[i-1], mRNA_string[i-1], delta_g_bb, delta_g)/(a_2() + a_2() + a_2() + a_2() + c_1(k_on, mRNA_template[i-1], mRNA_string[i-1], delta_g_bb, delta_g));
                    
                    if(random_number1 < alpha_2_add_aA){
                        //go to transition state 1, adding a 1(=aA) to the mRNA
                        transition_state = -1;
                        transition_state_string = "t1_aA";
                    }
                    else if(random_number1 < (alpha_2_add_aA + alpha_2_add_T)){
                        //go to transition state 1, adding a 2(=T) to the mRNA
                            transition_state = -2;
                            transition_state_string = "t1_T";
                    }
                    else if(random_number1 < (alpha_2_add_aA + alpha_2_add_T + alpha_2_add_G)){
                        //go to transition state 1, adding a 3(=G) to the mRNA
                            transition_state = -3;
                            transition_state_string = "t1_G";
                    }
                    else if(random_number1 < (alpha_2_add_aA + alpha_2_add_T + alpha_2_add_G + alpha_2_add_C)){
                        //go to transition state 1, adding a 4(=C) to the mRNA
                            transition_state = -4;
                            transition_state_string = "t1_C";
                    }
                    else{
                        //go back to transition state 2
                        transition_state = 2;
                        transition_state_string = "t2";
                    }
                }


                double random_number2 = distribution(generator);
                // if in trasition state -1, adding a 1(=aA) to the mRNA string
                if(transition_state == -1){
                    //define the probabilities for changing between states --> together the prob to add monomer 1 
                    double beta_2_add_1  = b_2(k_bb, mRNA_string[i], 1, delta_g_bb, delta_g)/(a_1(k_on, mRNA_template[i+1], 1, delta_g_bb, delta_g) + b_2(k_bb, mRNA_string[i], 1, delta_g_bb, delta_g));
                    double beta_1_add_1  = b_1()/(b_1() + c_2());
                    double gamma_2_add_1 = c_2()/(b_1() + c_2());
                    double prob_add_monomer_1 = (beta_2_add_1 * gamma_2_add_1)/ (1 - beta_1_add_1 *beta_2_add_1);
                    if(prob_add_monomer_1 > random_number2){
                        //s(M+1) = 1 --> add 1 to the mRNA string
                        mRNA_string.push_back(1);
                    }
                    else{
                        // go back to state s(M)
                    } 
                }

                // if in trasition state -2, adding a 2(=T) to the mRNA string
                else if(transition_state == -2){
                    //define the probabilities for changing between states --> together the prob to add monomer 1 
                    double beta_2_add_1  = b_2(k_bb, mRNA_string[i], 2, delta_g_bb, delta_g)/(a_1(k_on, mRNA_template[i+1], 2, delta_g_bb, delta_g) + b_2(k_bb, mRNA_string[i], 2, delta_g_bb, delta_g));
                    double beta_1_add_1  = b_1()/(b_1() + c_2());
                    double gamma_2_add_1 = c_2()/(b_1() + c_2());
                    double prob_add_monomer_1 = (beta_2_add_1 * gamma_2_add_1)/ (1 - beta_1_add_1 *beta_2_add_1);
                    if(prob_add_monomer_1 > random_number2){
                        //s(M+1) = 1 --> add 1 to the mRNA string
                        mRNA_string.push_back(2);
                    }
                    else{
                        // go back to state s(M)
                    } 
                }

                // if in trasition state -3, adding a 3(=G) to the mRNA string
                else if(transition_state == -3){
                    //define the probabilities for changing between states --> together the prob to add monomer 1 
                    double beta_2_add_1  = b_2(k_bb, mRNA_string[i], 3, delta_g_bb, delta_g)/(a_1(k_on, mRNA_template[i+1], 3, delta_g_bb, delta_g) + b_2(k_bb, mRNA_string[i], 3, delta_g_bb, delta_g));
                    double beta_1_add_1  = b_1()/(b_1() + c_2());
                    double gamma_2_add_1 = c_2()/(b_1() + c_2());
                    double prob_add_monomer_1 = (beta_2_add_1 * gamma_2_add_1)/ (1 - beta_1_add_1 *beta_2_add_1);
                    if(prob_add_monomer_1 > random_number2){
                        //s(M+1) = 1 --> add 1 to the mRNA string
                        mRNA_string.push_back(3);
                    }
                    else{
                        // go back to state s(M)
                    } 
                }

                // if in trasition state -4, adding a 1(=C) to the mRNA string
                else if(transition_state == -4){
                    //define the probabilities for changing between states --> together the prob to add monomer 1 
                    double beta_2_add_1  = b_2(k_bb, mRNA_string[i], 4, delta_g_bb, delta_g)/(a_1(k_on, mRNA_template[i+1], 4, delta_g_bb, delta_g) + b_2(k_bb, mRNA_string[i], 4, delta_g_bb, delta_g));
                    double beta_1_add_1  = b_1()/(b_1() + c_2());
                    double gamma_2_add_1 = c_2()/(b_1() + c_2());
                    double prob_add_monomer_1 = (beta_2_add_1 * gamma_2_add_1)/ (1 - beta_1_add_1 *beta_2_add_1);
                    if(prob_add_monomer_1 > random_number2){
                        //s(M+1) = 1 --> add 1 to the mRNA string
                        mRNA_string.push_back(4);
                    }
                    else{
                        // go back to state s(M)
                    } 
                }
                

                // if in trasition state 2, remove a monomer or stay in s(M)
                else if(transition_state == 2){
                    //define the probabilities for changing between states --> together the prob to add monomer 1 
                    double beta_2_remove  = b_2(k_bb, mRNA_string[i-1], mRNA_string[i], delta_g_bb, delta_g)/(a_1(k_on, mRNA_template[i], mRNA_string[i], delta_g_bb, delta_g) + b_2(k_bb, mRNA_string[i-1], mRNA_string[i], delta_g_bb, delta_g));
                    double beta_1_remove  = b_1()/(b_1() + c_2());
                    double gamma_2_remove = c_2()/(b_1() + c_2());
                    double prob_stay_s_M = gamma_2_remove/ (1 - beta_1_remove *beta_2_remove);
                    if(prob_stay_s_M > random_number2){
                        //stay in s(M)
                    }
                    else{
                        // remove s(M) monomer and move back to S(M-1)
                        mRNA_string.pop_back();
                        mRNA_error[i] = 0;
                    } 
                }


                //define for the string until the current position, S(M) the error.
                int error_true = error(mRNA_template[i], mRNA_string[i]);
                mRNA_error[i] = error_true;
            
                 
                mRNA_moves.template_string = mRNA_template;
                mRNA_moves.polymer = mRNA_string;
                mRNA_moves.length= mRNA_string.size();
                mRNA_moves.transition_state = transition_state_string;
                mRNA_moves.error = count(mRNA_error.begin(), mRNA_error.end(), -5)/(double)mRNA_string.size(); //error prob is #w/length polymer
                
                //save information to CSV file 
                std::vector<int> input = mRNA_moves.polymer;
                // loop through the array elements
                for(int j=0; j< input.size(); ++j){
                    fss << input.at(j);
                    }
                fss << "," << mRNA_moves.length <<"," << mRNA_moves.transition_state << ","  << mRNA_moves.error;
                fss << "\n";
                iteration += 1;
            } // end while loop --> make mRNa
        fss << "template string" << ",";
        std::vector<int> in = mRNA_moves.template_string;
                // loop through the array elements
                for(int j=0; j< in.size(); ++j){
                    fss << in.at(j) ;
                    }
        fss.close();

        // print the polymer computed
        // std::cout << "end polymer" << std::endl;
        // print_matrix(mRNA_string);
        // std::cout << "   " << std::endl;
        // std::cout << "template polymer:" << std::endl;
        // print_matrix(mRNA_template);
        // std::cout << mRNA_moves[iteration].error_old << std::endl;
        // std::cout << "   " << std::endl;

        // // // print path it took to come to this polymer (which states it passed)
        // std::cout << "path taken to get to this polymer:" << std::endl;
        // print_matrix_of_structures(mRNA_moves);

        // iteration_delta_g_bb +=1;
        } //end delta_g_bb loop
    } //end proption loop
}//end main


