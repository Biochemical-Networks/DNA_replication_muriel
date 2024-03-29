/* 
07-10-22, Muriel Louman, AMOLF

Grow mRNA, assumingh a randomly distributed template string of zeros and ones. Reproduicing the kinetic case of Jenny's paper. 
delta G_TT = 0, delta G_k = delta G. 

mRNA only consisting of zeros and ones.
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
    //typedef gives a way to define the type of struct somewhere without typing struct in front of it
    double concentration;
    double concentration_star;
    double delta_G;
    double delta_G_bb;
    double delta_G_k;
    } variable_definition_t ;

// define needed free energies and concentrations dependent on template en mRNA binding
variable_definition_t def_variables_given_ij(int i,int j, double G_bb, double G_k){
    variable_definition_t variable_definition;
    
    // G_bb = 0;
    if(i ==0 && j==0){
        variable_definition.concentration = 1;
        variable_definition.concentration_star= 1.5;
        variable_definition.delta_G = 5;
        variable_definition.delta_G_bb = G_bb;
        variable_definition.delta_G_k = G_k;
    }
    else if(i == 0 && j==1){
        variable_definition.concentration = 1;
        variable_definition.concentration_star= 1.5;
        variable_definition.delta_G = 5;
        variable_definition.delta_G_bb = G_bb;
        variable_definition.delta_G_k = 0;
    }
    else if(i == 1 && j==0){
        variable_definition.concentration = 1;
        variable_definition.concentration_star= 1.5;
        variable_definition.delta_G = 5;
        variable_definition.delta_G_bb = G_bb;
        variable_definition.delta_G_k = 0;
    }
    else {
        variable_definition.concentration = 1;
        variable_definition.concentration_star= 1.5;
        variable_definition.delta_G = 5;
        variable_definition.delta_G_bb = G_bb;
        variable_definition.delta_G_k = G_k;
    }
    return variable_definition;
}

// define rates
double a_2(double k_on, int i, int j, double G_bb, double G_k){
    variable_definition_t variables = def_variables_given_ij(i,j,G_bb, G_k);
    return k_on * variables.concentration * exp(-1 * variables.delta_G_k);
}

double a_1(double k_on, int i, int j, double G_bb, double G_k){
    variable_definition_t variables = def_variables_given_ij(i,j,G_bb, G_k);
    return k_on * exp(-1 * variables.delta_G - 1 * variables.delta_G_k);
}

double b_2(double k_bb){
    return k_bb;
}

double b_1(double k_bb, int j_M, int j_M1, double G_bb, double G_k){
    variable_definition_t variables = def_variables_given_ij(j_M, j_M1, G_bb, G_k);
    return k_bb * exp(-1 * variables.delta_G_bb);
}

double c_2(double k_on, int i, int j, double G_bb, double G_k){
    variable_definition_t variables = def_variables_given_ij(i,j,G_bb, G_k);
    return k_on * exp(-1 * variables.delta_G);
    }

double c_1(double k_on, int i, int j, double G_bb, double G_k){
    variable_definition_t variables = def_variables_given_ij(i,j,G_bb, G_k);
    return k_on * variables.concentration_star;
}


// define probabilities between states instead of rates--> done in loop


void print_matrix(std::vector<int> input){
    int n = sizeof(input)/sizeof(input[0]);
 
    // loop through the array elements
    for(int i=0; i< input.size(); ++i){
        std::cout << input.at(i) << ' ';
        }
}

//Create a structure variable called mRNA_move
struct structure_mRNA_move{
    std::vector<int> polymer_old;
    int length_old;
    double error_old;
    string transition_state;
    };


void print_matrix_of_structures(std::vector<structure_mRNA_move> mRNA_movess){
    // loop through the array elements, which are all structures defined by structure_mRNA_move
    for(int i=0; i< mRNA_movess.size(); ++i){
        cout << "creates new polymer" << endl;
        print_matrix(mRNA_movess[i].polymer_old);
        cout << "   " << endl;
        cout << "of length" << endl;
        cout << mRNA_movess[i].length_old << endl;
        cout << "with error probability"  << endl;
        cout << mRNA_movess[i].error_old << endl;
        cout << "moves to transition state" << endl;
        cout << mRNA_movess[i].transition_state << endl;
        }
}

void save_matrix_of_structures(std::vector<structure_mRNA_move> mRNA_movess, std::string filename){
    //saves the structures of different computations in a matrix
    ofstream fss;
    fss.open(filename.c_str());
    // write the file headers
    fss << "start polymer" << "," << "length polymer" <<"," << "error probability" << "," << "transition state"  << std::endl;
    // loop through the array elements, which are all structures defined by structure_mRNA_move
    for(int i=0; i< mRNA_movess.size(); ++i){
        std::vector<int> input = mRNA_movess[i].polymer_old;
        int n = sizeof(input)/sizeof(input[0]);
 
        // loop through the array elements
        for(int j=0; j< input.size(); ++j){
            fss << input.at(j);
            }
        fss << "," << mRNA_movess[i].length_old <<"," << mRNA_movess[i].error_old << "," << mRNA_movess[i].transition_state << std ::endl;
    }
    fss.close();
}

string date(time_t now){
    //define date in this way: 20221011, which is 11 October 2022
    tm *ltm = localtime(&now);
    string datum = std::to_string(1900 + ltm->tm_year) + std::to_string(1 + ltm->tm_mon) + std::to_string(ltm->tm_mday);
    return datum;
}

std::vector<int> template_string(int string_length){
    // create a (uniformly) randomized template string 
    std::vector<int> string;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    for(int j=0; j<string_length; j++){
        double random_number = distribution(generator);
        int rounded_random_number = round(random_number);
        string.push_back(rounded_random_number);
        }
    return string;
}

int error(int i, int j){
        //error vector with -5 is incorrect, 5 is correct
        int E = 0;
        if(i==j){
            E = -5;}
        else{
            E = 5;}
    return E;
    }


int main(){
    // define variables
    double k_on = 1.0;
    double k_bb = 1.0;
    int length_mRNA = 150;

    // random number generator for uniform distr between 0 and 1
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    double L_delta_g_bb[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; //we would like to have more but takes computer time
    double L_delta_g_k[] = {0, 2, 4, 6, 8, 10};

    // define template string
    std::vector<int> mRNA_template = template_string(length_mRNA +1);
    
    
    //loop over different g_r = g_w and g_bb
    for (double delta_g_bb : L_delta_g_bb) {
        for (double delta_g_k : L_delta_g_k) {
            //error vector with -5 is incorrect, 5 is correct
            std::vector<int> mRNA_error(length_mRNA,0);

            // mRNA string we are going to grow, start with a zero on the first place
            std::vector<int> mRNA_string;
            mRNA_string.push_back(0); 

            // keep track of the transition states and the length of mRNA trough the code. This is also done in mRNA moves. start with a mRNA length of 1
            std::vector<int> mRNA_string_states_length;
            mRNA_string_states_length.push_back(1);

            //define variables later used in code
            int i = 0;
            int transition_state = 0;
            int iteration = 0;
            std::vector<structure_mRNA_move> mRNA_moves;
            
            while(mRNA_string.size()< length_mRNA){
            // for(int j =0;j<20;j++){
                // create a random number betwen 0 and 1
                double random_number1 = distribution(generator);

                // keep track of the current position. S(M)
                i = mRNA_string.size() -1;
                
                //define for the current position, S(M) the error.
                int error_true = error(mRNA_template[i], mRNA_string[i]);
                mRNA_error[i] = error_true;
            
                // create first element of vector which is the structure we defined before
                mRNA_moves.push_back(structure_mRNA_move());
                // // add elements in the structure of this iteration
                mRNA_moves[iteration].polymer_old = mRNA_string;
                mRNA_moves[iteration].length_old= mRNA_string.size();
                mRNA_moves[iteration].error_old = count(mRNA_error.begin(), mRNA_error.end(), -5)/(double)mRNA_string.size(); // error prob is #w/length polymer
                

                                       
                // if string is shorter than 2 monomers, always go add 1(50%) or 0(50%) to the mRNA
                if(mRNA_string.size() == 1){              
                    if(random_number1 <0.5){
                        //add a 0 to the mRNA
                        transition_state = -1;
                        mRNA_string_states_length.push_back(transition_state);
                        mRNA_moves[iteration].transition_state = "t1_0";
                    }
                    else{
                        //add a 1 to the mRNA
                        transition_state = 1;
                        mRNA_string_states_length.push_back(transition_state);
                        mRNA_moves[iteration].transition_state = "t1_1";
                    }
                }

                else{
                    // define probabilities needed for going to t1 or t2
                    double alpha_2_add_0 = a_2(k_on, mRNA_template[i+1], 0, delta_g_bb, delta_g_k)/(a_2(k_on, mRNA_template[i+1], 0, delta_g_bb, delta_g_k) + a_2(k_on, mRNA_template[i+1], 1, delta_g_bb, delta_g_k) + c_1(k_on, mRNA_template[i-1], mRNA_string[i-1], delta_g_bb, delta_g_k));
                    double alpha_2_add_1 = a_2(k_on, mRNA_template[i+1], 1, delta_g_bb, delta_g_k)/(a_2(k_on, mRNA_template[i+1], 0, delta_g_bb, delta_g_k) + a_2(k_on, mRNA_template[i+1], 1, delta_g_bb, delta_g_k) + c_1(k_on, mRNA_template[i-1], mRNA_string[i-1], delta_g_bb, delta_g_k));
                    double gamma_1 = c_1(k_on, mRNA_template[i-1], mRNA_string[i-1], delta_g_bb, delta_g_k)/(a_2(k_on, mRNA_template[i+1], 0, delta_g_bb, delta_g_k) + a_2(k_on, mRNA_template[i+1], 1, delta_g_bb, delta_g_k) + c_1(k_on, mRNA_template[i-1], mRNA_string[i-1], delta_g_bb, delta_g_k));
                    
                    if(random_number1 < alpha_2_add_0){
                        //go to transition state 1, adding a 0 to the mRNA
                        transition_state = -1;
                        mRNA_string_states_length.push_back(transition_state);
                        mRNA_moves[iteration].transition_state = "t1_0";
                    }

                    else if(random_number1 < (alpha_2_add_0 + alpha_2_add_1)){
                        //go to transition state 1, adding a 1 to the mRNA
                            transition_state = 1;
                            mRNA_string_states_length.push_back(transition_state);
                            mRNA_moves[iteration].transition_state = "t1_1";
                    }
                    else{
                        //go back to transition state 2
                        transition_state = 2;
                        mRNA_string_states_length.push_back(transition_state);
                        mRNA_moves[iteration].transition_state = "t2";
                    }
                }


                double random_number2 = distribution(generator);
                // if in trasition state 1, adding a 1 to the mRNA string
                if(transition_state == 1){
                    //define the probabilities for changing between states --> together the prob to add monomer 1 
                    double beta_2_add_1 = b_2(k_bb)/(a_1(k_on, mRNA_template[i+1], 1, delta_g_bb, delta_g_k) + b_2(k_bb));
                    double beta_1_add_1 = b_1(k_bb, mRNA_string[i], 1, delta_g_bb, delta_g_k)/(b_1(k_bb, mRNA_string[i], 1, delta_g_bb, delta_g_k) + c_2(k_on, mRNA_template[i], mRNA_string[i], delta_g_bb, delta_g_k));
                    double gamma_2_add_1 = c_2(k_on, mRNA_template[i], mRNA_string[i], delta_g_bb, delta_g_k)/(b_1(k_bb, mRNA_string[i], 1, delta_g_bb, delta_g_k) + c_2(k_on, mRNA_template[i], mRNA_string[i], delta_g_bb, delta_g_k));
                    double prob_add_monomer_1 = (beta_2_add_1 * gamma_2_add_1)/ (1 - beta_1_add_1 *beta_2_add_1);
                    if(prob_add_monomer_1 > random_number2){
                        //s(M+1) = 1 --> add 1 to the mRNA string
                        mRNA_string.push_back(1);
                    }
                    else{
                        // go back to state s(M)
                    } 
                }
                
                // if in trasition state 1, adding a 0 to the mRNA string
                else if(transition_state == -1){
                    //define the probabilities for changing between states --> together the prob to add monomer 1 
                    double beta_2_add_0 = b_2(k_bb)/(a_1(k_on, mRNA_template[i+1], 0, delta_g_bb, delta_g_k) + b_2(k_bb));
                    double beta_1_add_0 = b_1(k_bb, mRNA_string[i], 0, delta_g_bb, delta_g_k)/(b_1(k_bb, mRNA_string[i], 0, delta_g_bb, delta_g_k) + c_2(k_on, mRNA_template[i], mRNA_string[i], delta_g_bb, delta_g_k));
                    double gamma_2_add_0 = c_2(k_on, mRNA_template[i], mRNA_string[i], delta_g_bb, delta_g_k)/(b_1(k_bb, mRNA_string[i], 0, delta_g_bb, delta_g_k) + c_2(k_on, mRNA_template[i], mRNA_string[i], delta_g_bb, delta_g_k));
                    double prob_add_monomer_0 = (beta_2_add_0 * gamma_2_add_0)/ (1 - beta_1_add_0 *beta_2_add_0);
                    if(prob_add_monomer_0 > random_number2){
                        //s(M+1) = 0 --> add 0 to the mRNA string
                        mRNA_string.push_back(0);
                    }
                    else{
                        // go back to state s(M)
                    } 
                }

                // if in trasition state 2, remove a monomer or stay in s(M)
                else if(transition_state == 2){
                    //define the probabilities for changing between states --> together the prob to add monomer 1 
                    double beta_2_remove = b_2(k_bb)/(a_1(k_on, mRNA_template[i], mRNA_string[i], delta_g_bb, delta_g_k) + b_2(k_bb));
                    double beta_1_remove = b_1(k_bb, mRNA_string[i-1], mRNA_string[i], delta_g_bb, delta_g_k)/(b_1(k_bb, mRNA_string[i-1], mRNA_string[i], delta_g_bb, delta_g_k) + c_2(k_on, mRNA_template[i-1], mRNA_string[i-1], delta_g_bb, delta_g_k));
                    double gamma_2_remove = c_2(k_on, mRNA_template[i-1], mRNA_string[i-1], delta_g_bb, delta_g_k)/(b_1(k_bb, mRNA_string[i-1], mRNA_string[i], delta_g_bb, delta_g_k) + c_2(k_on, mRNA_template[i-1], mRNA_string[i-1], delta_g_bb, delta_g_k));
                    double prob_stay_s_M = gamma_2_remove/ (1 - beta_1_remove *beta_2_remove);
                    if(prob_stay_s_M > random_number2){
                        //stay in s(M)
                    }
                    else{
                        // remove s(M) monomer and move back to S(M-1)
                        mRNA_string.pop_back();
                    } 
                }

                mRNA_string_states_length.push_back(mRNA_string.size());

                // //prints to check if program works right. from up to down: size polymer, polymer computed, 
                // //polymer length and transition states it passed up intill this iteration
                // std::cout << mRNA_string.size() << std::endl;
                // print_matrix(mRNA_string);
                // std::cout << " " << std::endl;
                // print_matrix(mRNA_template);
                // std::cout << " " << std::endl;
                // print_matrix(mRNA_error);
                // std::cout << mRNA_moves[iteration].error_old << std::endl;
                // print_matrix(mRNA_string_states_length);
                // std::cout << "new computation"  << std::endl;
                iteration += 1;
            }

        // add end polymer to array of structures
        mRNA_moves.push_back(structure_mRNA_move());    
        mRNA_moves[iteration].polymer_old = mRNA_string;
        mRNA_moves[iteration].length_old= mRNA_string.size();
        mRNA_moves[iteration].error_old = count(mRNA_error.begin(), mRNA_error.end(), -5)/(double)mRNA_string.size();
        mRNA_moves[iteration].transition_state = "end computation";

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

        //save the polymer and the path it took to get there
        int rounded_delta_g_bb = round(delta_g_bb);
        int rounded_delta_g_k = round(delta_g_k);
        std:cout << " done" << std::endl;
        string date_now = date(time(0));
        string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/221007output/";
        string filename_output = path + date_now + "_delta_g_bb_" + std::to_string(rounded_delta_g_bb) + "_delta_g_k_" + std::to_string(rounded_delta_g_k) + ".csv";
        save_matrix_of_structures(mRNA_moves, filename_output);

        }
    }
}

