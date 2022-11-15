/* 
06-10-22, Muriel Louman, AMOLF

Grow mRNA, assumingh template string of only ones. The energetic case of Jenny's paper. 
delta G_k = 0, delta G_TT = delta G

mRNA only consisting of zeros and ones.
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
    tm *ltm = localtime(&now);
    string datum = std::to_string(1900 + ltm->tm_year) + std::to_string(1 + ltm->tm_mon) + std::to_string(ltm->tm_mday);
    return datum;
}

std::vector<int> template_string(int string_length){
    std::vector<int> string;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    for(int j; j<string_length; j++){
        double random_number = distribution(generator);
        int rounded_random_number = round(random_number);
        string.push_back(rounded_random_number);
        }
    return string;
}


int main(){
    // define variables
    double k_on = 1.0;
    double k_bb = 1;
    double r_con = 1.0;
    double w_con = 1.0;
    double r_con_star = 1.5;
    double w_con_star = 1.5;
    int length_mRNA = 300;
    

    // random number generator for uniform distr between 0 and 1
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    double L_delta_g_bb[] = {0,1,2,3,4,5,6,7,8,9,10}; //we would like to have more but takes computer time
    double L_delta_g_rw[] = {0,2,4,6,8,10};
    
    

    // define template string
    // std::vector<int> mRNA_template = template_string(length_mRNA);

    //loop over different g_r = g_w and g_bb
    for (double delta_g_bb : L_delta_g_bb) {
        for (double delta_g_rw : L_delta_g_rw) {
            
            // // G_TT defenition
            // double delta_g_w = -1 * (delta_g_rw/2);
            // double delta_g_r = delta_g_rw/2;
            // string model1 = "G2";
            
            // double x = -2;
            // double delta_g_w = x;
            // double delta_g_r = x + delta_g_rw;
            // int rounded_x = round(x);
            // string model1 = "Gx" + std::to_string(rounded_x);

            //define rates
            // double a_1r = k_on * exp(-1 * delta_g_r);
            // double a_1w = k_on * exp(-1 * delta_g_w);
            // double a_2r = k_on * r_con;
            // double a_2w = k_on * w_con;

            // double b_1 = k_bb * exp(-1 * delta_g_bb);
            // double b_2 = k_bb;
            // string model2 = "b1exp";

            // // double b_1 = k_bb;
            // // double b_2 = k_bb * exp(delta_g_bb);
            // // string model2 = "b2exp";

            // double c_1r = k_on * r_con_star; 
            // double c_1w = k_on * w_con_star; 
            // double c_2r = k_on * exp(-1* delta_g_r); 
            // double c_2w = k_on * exp(-1* delta_g_w); 

            double a_1r = 1;
            double a_1w = exp(delta_g_rw);
            double a_2r = 1;
            double a_2w = 1;

            double b_1 = exp(-1 * delta_g_bb);
            double b_2 = 1;
            string model1 = "b1exp";

            // double b_1 = 1;
            // double b_2 = exp(delta_g_bb);
            // string model1 = "b2exp";

            double c_1r = 1; 
            double c_1w = exp(-1 * delta_g_rw); 
            double c_2r = 1; 
            double c_2w = 1; 
            string model2 = "Jenny";

            // define probabilities between states instead of rates
            double alpha_1r = a_1r/(a_1r + b_2);   
            double alpha_1w = a_1w/(a_1w + b_2); 
            double beta_2r = b_2/(a_1r + b_2);
            double beta_2w = b_2/(a_1w + b_2);
            double gamma_1r = c_1r/(a_2r + a_2w + c_1r);
            double gamma_1w = c_1w/(a_2r + a_2w + c_1w); 
            double beta_1r = b_1/(c_2r + b_1);
            double beta_1w = b_1/(c_2w + b_1);
            double gamma_2r = c_2r/(c_2r + b_1);
            double gamma_2w = c_2w/(c_2w + b_1);
            //NB: alpha_x2y means: given s(M-1) is x what is the prob for adding y in t1
            double alpha_r2r = a_2r/(a_2r + a_2w + c_1r);
            double alpha_r2w = a_2w/(a_2r + a_2w + c_1r);
            double alpha_w2r = a_2r/(a_2r + a_2w + c_1w);   
            double alpha_w2w = a_2w/(a_2r + a_2w + c_1w);

            std::vector<int> mRNA_string;
            mRNA_string.push_back(0);
            std::vector<int> mRNA_string_states_length;
            mRNA_string_states_length.push_back(1);
            int i = 0;
            int transition_state = 0;
            int iteration = 0;
            std::vector<structure_mRNA_move> mRNA_moves;
            
            while(mRNA_string.size()< length_mRNA){
            // for(int j =0;j<20;j++){
                // Go to transition state 1r or 1w or 2, first MC step
                // create a random number betwen 0 and 1
                double random_number1 = distribution(generator);
                i = mRNA_string.size() -1;


                // create first element of vector which is the structure we defined before
                mRNA_moves.push_back(structure_mRNA_move());
                // // add elements in the structure of this iteration
                mRNA_moves[iteration].polymer_old = mRNA_string;
                mRNA_moves[iteration].length_old= mRNA_string.size();
                mRNA_moves[iteration].error_old = count(mRNA_string.begin(), mRNA_string.end(), 1)/(double)mRNA_string.size(); // error prob is #w/length polymer
                
                            
                // if string is shorter than 2 monomers, always go to t1r(50%) or t1w(50%)
                if(mRNA_string.size() == 1){
                    if(random_number1 <0.5){
                        transition_state = 1;
                        mRNA_string_states_length.push_back(transition_state);
                        mRNA_moves[iteration].transition_state = "t1r";
                    }
                    else{
                        transition_state = -1;
                        mRNA_string_states_length.push_back(transition_state);
                        mRNA_moves[iteration].transition_state = "t1w";
                    }
                }
                else{
                    // if s(M-1) = r
                    if(mRNA_string[i-1] == 0){
                        if(random_number1 < alpha_r2r){
                            transition_state = 1;
                            mRNA_string_states_length.push_back(transition_state);
                            mRNA_moves[iteration].transition_state = "t1r";
                        }
                        else if(random_number1 < (alpha_r2r + alpha_r2w)){
                            transition_state = -1;
                            mRNA_string_states_length.push_back(transition_state);
                            mRNA_moves[iteration].transition_state = "t1w";
                        }
                        else{
                            transition_state = 2;
                            mRNA_string_states_length.push_back(transition_state);
                            mRNA_moves[iteration].transition_state = "t2r";
                        }
                    }
                    // if s(M-1) = w
                    if(mRNA_string[i-1]==1){
                        if(random_number1 < alpha_w2r){
                            transition_state = 1;
                            mRNA_string_states_length.push_back(transition_state);
                            mRNA_moves[iteration].transition_state = "t1r";
                        }
                        else if(random_number1 < (alpha_w2r + alpha_w2w)){
                            transition_state = -1;
                            mRNA_string_states_length.push_back(transition_state);
                            mRNA_moves[iteration].transition_state = "t1w";
                        }
                        else{
                            transition_state = -2;
                            mRNA_string_states_length.push_back(transition_state);
                            mRNA_moves[iteration].transition_state = "t2w";
                        }
                    }         
                }

                double random_number2 = distribution(generator);
                // if in transition state 1r
                if(transition_state == 1){
                    // if s(M) = r
                    if(mRNA_string[i] == 0){
                        double prob_mr_m_1r = (beta_2r * gamma_2r)/ (1 - beta_1r *beta_2r);
                        if(prob_mr_m_1r > random_number2){
                            //s(M+1) = r --> add r
                            mRNA_string.push_back(0);
                        }
                        else{
                            // go back to state s(M) = r
                        } 
                    }
                    // if s(M) = w
                    if(mRNA_string[i] == 1){
                        double prob_mw_m_1r = (beta_2r * gamma_2w)/ (1 - beta_1w *beta_2r);
                        if(prob_mw_m_1r > random_number2){
                            //s(M+1) = r --> add r
                            mRNA_string.push_back(0);
                        }
                        else{
                            // go back to state s(M) = w
                        } 
                    }
                }

                // if in transition state 1w
                if(transition_state == -1){
                    // if s(M) = r
                    if(mRNA_string[i] == 0){
                        double prob_mr_m_1w = (beta_2w * gamma_2r)/ (1 - beta_1r *beta_2w);
                        if(prob_mr_m_1w > random_number2){
                            //s(M+1) = w --> add w
                            mRNA_string.push_back(1);
                        }
                        else{
                            // go back to state s(M) = r
                        } 
                    }
                    // if s(M) = w
                    if(mRNA_string[i] == 1){
                        double prob_mw_m_1w = (beta_2w * gamma_2w)/ (1 - beta_1w *beta_2w);
                        if(prob_mw_m_1w > random_number2){
                            //s(M+1) = w --> add w
                            mRNA_string.push_back(1);
                        }
                        else{
                            // go back to state s(M) = w
                        } 
                    }
                }

                // if in transition state 2r
                if(transition_state == 2){
                    // if s(M) = r
                    if(mRNA_string[i] == 0){
                        double prob_2mr_m_1r = (alpha_1r * beta_1r)/(1 - beta_1r * beta_2r);
                        if(prob_2mr_m_1r > random_number2){
                            //remove s(M) = r --> go to s(M-1) = r
                            mRNA_string.pop_back();
                        }
                        else{
                            // go back to state s(M) = r
                        } 
                    }
                    // if s(M) = w
                    if(mRNA_string[i] == 1){
                        double prob_2mw_m_1r = (alpha_1w * beta_1r)/(1 - beta_1r * beta_2w);
                        if(prob_2mw_m_1r > random_number2){
                            //remove s(M) = w --> go to s(M-1) = r
                            mRNA_string.pop_back();
                        }
                        else{
                            // go back to state s(M) = w
                        } 
                    }
                }

                // if in transition state 2w
                if(transition_state == -2){
                    // if s(M) = r
                    if(mRNA_string[i] == 0){
                        double prob_2mr_m_1w = (alpha_1r * beta_1w)/(1 - beta_1w * beta_2r);
                        if(prob_2mr_m_1w > random_number2){
                            //remove s(M) = r --> go to s(M-1) = w
                            mRNA_string.pop_back();
                        }
                        else{
                            // go back to state s(M) = r
                        } 
                    }
                    // if s(M) = w
                    if(mRNA_string[i] == 1){
                        double prob_2mw_m_1w = (alpha_1w * beta_1w)/(1 - beta_1w * beta_2w);
                        if(prob_2mw_m_1w > random_number2){
                            //remove s(M) = w --> go to s(M-1) = w
                            mRNA_string.pop_back();
                        }
                        else{
                            // go back to state s(M) = w
                        } 
                    }
                }

                
                mRNA_string_states_length.push_back(mRNA_string.size());

                // //prints to check if program works right. from up to down: size polymer, polymer computed, 
                // //polymer length and transition states it passed up intill this iteration
                // std::cout << mRNA_string.size() << std::endl;
                // print_matrix(mRNA_string);
                // print_matrix(mRNA_string_states_length);
                // std::cout << "new computation"  << std::endl;
                iteration += 1;
            }

        // add end polymer to array of structures
        mRNA_moves.push_back(structure_mRNA_move());    
        mRNA_moves[iteration].polymer_old = mRNA_string;
        mRNA_moves[iteration].length_old= mRNA_string.size();
        mRNA_moves[iteration].error_old = count(mRNA_string.begin(), mRNA_string.end(), 1)/(double)mRNA_string.size();
        mRNA_moves[iteration].transition_state = "end computation";

        // print the polymer computed
        std::cout << "end polymer" << std::endl;
        print_matrix(mRNA_string);
        // std::cout << "   " << std::endl;
        // std::cout << "   " << std::endl;

        // // // print path it took to come to this polymer (which states it passed)
        // std::cout << "path taken to get to this polymer:" << std::endl;
        // print_matrix_of_structures(mRNA_moves);

        //save the polymer and the path it took to get there
        int rounded_delta_g_bb = round(delta_g_bb);
        int rounded_delta_g_k = round(delta_g_rw);
        int rounded_k_bb = round(k_bb);
        std:cout << " done" << std::endl;
        string date_now = date(time(0));
        string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/221109output/";
        string filename_output = path + date_now + "model" + model1 + "_" + model2 + "_delta_g_bb_" + std::to_string(rounded_delta_g_bb) + "_delta_g_k_" + std::to_string(rounded_delta_g_k) + "_k_bb_"+ std::to_string(rounded_k_bb) + ".csv";
        save_matrix_of_structures(mRNA_moves, filename_output);

        }
    }
}

