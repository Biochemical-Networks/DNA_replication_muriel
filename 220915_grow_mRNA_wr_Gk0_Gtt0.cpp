#include <iostream>
#include <cmath> // math functions
#include <vector> // vector container
#include <string> // text strings
#include <fstream> // output to file
#include <sstream> 
#include <random> // random number generator
#include <bits/stdc++.h> //for count function
using namespace std;

// define rates
double a_1rw_func(double k_on, double delta_g_rw){
    return k_on * exp(-1 * delta_g_rw);}

double a_2rw_func(double k_on, double rw_con){
    return k_on * rw_con;}

double b_1_func(double k_bb, double delta_g_bb){
    return k_bb * exp(-1 * delta_g_bb);}

double b_2_func(double k_bb){
    return k_bb;}

double c_1rw_func(double k_on, double rw_con_star){
    return k_on * rw_con_star;}

double c_2rw_func(double k_on, double delta_g_rw){
    return k_on * exp(-1* delta_g_rw);}



// define probabilities between states instead of rates
double prob_alpha_1rw(double a_1rw, double b_2){
    // r if we add a right one (in t1), w if we are adding a wrong one (in t1)
    return a_1rw/(a_1rw + b_2);}

double prob_beta_2rw(double a_1rw, double b_2){
    // r if we add a right one (in t1), w if we are adding a wrong one (in t1)
    return b_2/(a_1rw + b_2);}

double prob_alpha_2rw_given_sM_1_r(double a_2rw, double a_2r, double a_2w, double c_1r){
    //given s(M-1) = r, rw beacomes: r if we are adding a right one (in t1), w if we are adding a wrong one (in t1)
    return a_2rw/(a_2r + a_2w + c_1r);}

double prob_alpha_2rw_given_sM_1_w(double a_2rw, double a_2r, double a_2w, double c_1w){
    //given s(M-1) = w, w beacomes: r if we are adding a right one (in t1), w if we are adding a wrong one (in t1)
    return a_2rw/(a_2r + a_2w + c_1w);}

double prob_gamma_1rw(double a_2r, double a_2w, double c_1rw){
    // rw becomes: r if we let go of a right one(in t2), w if we let go of a wrong one (in t2)
    return c_1rw/(a_2r + a_2w + c_1rw);}

double prob_beta_1rw(double b_1, double c_2rw){
    // rw becomes: r if we let go of a right one(in t2), w if we let go of a wrong one (in t2)
    return b_1/(c_2rw + b_1);}

double prob_gamma_2rw(double b_1, double c_2rw){
    // rw becomes: r if we let go of a right one(in t2), w if we let go of a wrong one (in t2)
    return c_2rw/(c_2rw + b_1);}


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


int main(){
    // define variables
    double k_on = 1.0;
    double k_bb = 1.0;
    double r_con = 1.0;
    double w_con = 1.0;
    double r_con_star = 1.5;
    double w_con_star = 1.5;
    int length_mRNA = 300;

    // random number generator for uniform distr between 0 and 1
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    double L_delta_g_bb[] = {0,1,2,3,4,5,6,7,8,9,10};
    double L_delta_g_r[] = {0,2,4,6,8,10};
    
    //loop over different g_r = g_w and g_bb
    for (double delta_g_bb : L_delta_g_bb) {
        for (double delta_g_r : L_delta_g_r) {
            
            double delta_g_w = delta_g_r;
            double a_1r = a_1rw_func(k_on, delta_g_r);
            double a_1w = a_1rw_func(k_on, delta_g_w);
            double a_2r = a_2rw_func(k_on, r_con);
            double a_2w = a_2rw_func(k_on, w_con);
            double b_1 = b_1_func(k_bb, delta_g_bb);
            double b_2 = b_2_func(k_bb);
            double c_1r = c_1rw_func(k_on,r_con_star);
            double c_1w = c_1rw_func(k_on,w_con_star);
            double c_2r = c_2rw_func(k_on, delta_g_r);
            double c_2w = c_2rw_func(k_on, delta_g_w);

            double alpha_1r = prob_alpha_1rw(a_1r, b_2);
            double alpha_1w = prob_alpha_1rw(a_1w, b_2);
            double beta_2r = prob_beta_2rw(a_1r, b_2);
            double beta_2w = prob_beta_2rw(a_1w, b_2);
            double gamma_1r = prob_gamma_1rw(a_2r, a_2w, c_1r);
            double gamma_1w = prob_gamma_1rw(a_2r, a_2w, c_1w);
            double beta_1r = prob_beta_1rw(b_1, c_2r);
            double beta_1w = prob_beta_1rw(b_1, c_2w);
            double gamma_2r = prob_gamma_2rw(b_1, c_2r);
            double gamma_2w = prob_gamma_2rw(b_1, c_2w);
            //NB: alpha_x2y means: given s(M-1) is x what is the prob for adding y in t1
            double alpha_r2r = prob_alpha_2rw_given_sM_1_r(a_2r, a_2r, a_2w, c_1r);
            double alpha_r2w = prob_alpha_2rw_given_sM_1_r(a_2w, a_2r, a_2w, c_1r);
            double alpha_w2r = prob_alpha_2rw_given_sM_1_w(a_2r, a_2r, a_2w, c_1w);
            double alpha_w2w = prob_alpha_2rw_given_sM_1_w(a_2w, a_2r, a_2w, c_1w);

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
                        else if(random_number1>=alpha_r2r && random_number1 < (alpha_r2r + alpha_r2w)){
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
                        else if(random_number1>=alpha_w2r && random_number1 < (alpha_w2r + alpha_w2w)){
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
                        double prob_mr_m_1r = (alpha_1r * beta_1r)/(1 - beta_1r * beta_2r);
                        if(prob_mr_m_1r > random_number2){
                            //remove s(M) = r --> go to s(M-1) = r
                            mRNA_string.pop_back();
                        }
                        else{
                            // go back to state s(M) = r
                        } 
                    }
                    // if s(M) = w
                    if(mRNA_string[i] == 1){
                        double prob_mw_m_1r = (alpha_1w * beta_1r)/(1 - beta_1r * beta_2w);
                        if(prob_mw_m_1r > random_number2){
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
                        double prob_mr_m_1w = (alpha_1r * beta_1w)/(1 - beta_1w * beta_2r);
                        if(prob_mr_m_1w > random_number2){
                            //remove s(M) = r --> go to s(M-1) = w
                            mRNA_string.pop_back();
                        }
                        else{
                            // go back to state s(M) = r
                        } 
                    }
                    // if s(M) = w
                    if(mRNA_string[i] == 1){
                        double prob_mw_m_1w = (alpha_1w * beta_1w)/(1 - beta_1w * beta_2w);
                        if(prob_mw_m_1w > random_number2){
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
        std::cout << "   " << std::endl;
        std::cout << "   " << std::endl;

        // // print path it took to come to this polymer (which states it passed)
        std::cout << "path taken to get to this polymer:" << std::endl;
        print_matrix_of_structures(mRNA_moves);

        //save the polymer and the path it took to get there
        int rounded_delta_g_bb = round(delta_g_bb);
        int rounded_delta_g_r = round(delta_g_r);
        std:cout << " done" << std::endl;
        string date_now = date(time(0));
        string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/220915output/";
        string filename_output = path + date_now + "_delta_g_bb_" + std::to_string(rounded_delta_g_bb) + "_delta_g_r_" + std::to_string(rounded_delta_g_r) + ".csv";
        save_matrix_of_structures(mRNA_moves, filename_output);

        }
    }
}

