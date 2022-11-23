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

string date(time_t now){
    tm *ltm = localtime(&now);
    string datum = std::to_string(1900 + ltm->tm_year) + std::to_string(1 + ltm->tm_mon) + std::to_string(ltm->tm_mday);
    return datum;
}

//Create a structure variable called mRNA_move
struct structure_mRNA_move{
    std::vector<int> polymer;
    string transition_states;
    int length;
    double error;
    
    };

void print_matrix(std::vector<int> input){
    int n = sizeof(input)/sizeof(input[0]);
 
    // loop through the array elements
    for(int i=0; i< input.size(); ++i){
        std::cout << input.at(i) << ' ';
        }
}

void print_matrix_of_structures(std::vector<structure_mRNA_move> mRNA_movess){
    // loop through the array elements, which are all structures defined by structure_mRNA_move
    for(int i=0; i< mRNA_movess.size(); ++i){
        cout << "creates new polymer" << endl;
        print_matrix(mRNA_movess[i].polymer);
        cout << "   " << endl;
        cout << "transition states used" << endl;
        cout << mRNA_movess[i].transition_states << endl;
        cout << "of length" << endl;
        cout << mRNA_movess[i].length << endl;
        cout << "with error probability"  << endl;
        cout << mRNA_movess[i].error << endl;
        }
}

void save_matrix_of_structures(std::vector<structure_mRNA_move> mRNA_movess, std::string filename){
    ofstream fss;
    fss.open(filename.c_str());
    // write the file headers
    fss << "polymer" << "," << "transition states used" << "," << "length polymer" <<"," << "error probability"  << std::endl;
    // loop through the array elements, which are all structures defined by structure_mRNA_move
    for(int i=0; i< mRNA_movess.size(); ++i){
        std::vector<int> input = mRNA_movess[i].polymer;
        int n = sizeof(input)/sizeof(input[0]);
 
        // loop through the array elements
        for(int j=0; j< input.size(); ++j){
            fss << input.at(j);
            }

        fss << "," << mRNA_movess[i].transition_states ;
        fss << "," << mRNA_movess[i].length <<"," << mRNA_movess[i].error << std ::endl;
    }
    fss.close();
}
int main(){
    // define variables
    double k_on = 1.0;
    double k_bb = 1;
    double r_con = 1.0;
    double w_con = 1.0;
    double r_con_star = 1.5;
    double w_con_star = 1.5;
    

    // random number generator for uniform distr between 0 and 1
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    double L_delta_g_pol[] ={50};//{1,2,3,4,5,6,7,8,9,10}; //we would like to have more but takes computer time
    double L_delta_g_tt[] = {2};//{2,4,6,8,10};
    double error_prob = 0;
    int length_mRNA = 3000;
    

    // define template string
    // std::vector<int> mRNA_template = template_string(length_mRNA);

    //loop over different g_r = g_w and g_bb
    for (double delta_g_pol : L_delta_g_pol) {
        for (double delta_g_tt : L_delta_g_tt) {
            
            // define rates 
            // string model_1 = "r3";
            // double delta_g_w = 2;
            // double delta_g_r = delta_g_w + delta_g_tt;

            string model_1 = "r2";
            double delta_g_r = delta_g_tt/2;
            double delta_g_w = -1 * (delta_g_tt/2);

            double a_2r = k_on * r_con;
            double a_1r = k_on * exp(-1 * delta_g_r);
            double a_2w = k_on * w_con;
            double a_1w = k_on * exp(-1 * delta_g_w);
            double b2 = k_bb;
            double b1 = k_bb * exp(-delta_g_pol);
            double c_2r = k_on * exp(-1 * delta_g_r);
            double c_1r = k_on * r_con_star;
            double c_2w = k_on * exp(-delta_g_w);
            double c_1w = k_on * w_con_star;
            
            // double a_2r = 1;
            // double a_1r = 1;
            // double a_2w = 1;
            // double a_1w = exp(delta_g_tt);
            // double b2 = 1;
            // double b1 = exp(-delta_g_pol);
            // double c_2r = 1;
            // double c_1r = 1;
            // double c_2w = 1;
            // double c_1w = exp(-delta_g_tt);
            // string model_1 = "r1";


            // define probabilities between states instead of rates
            double alpha_1r = a_1r/(a_1r + b2);   
            double alpha_1w = a_1w/(a_1w + b2); 
            double beta_2r = b2/(a_1r + b2);
            double beta_2w = b2/(a_1w + b2);
            // double gamma_1r = c_1r/(a_2r + a_2w + c_1r);
            // double gamma_1w = c_1w/(a_2r + a_2w + c_1w); 
            double beta_1r = b1/(c_2r + b1);
            double beta_1w = b1/(c_2w + b1);
            double gamma_2r = c_2r/(c_2r + b1);
            double gamma_2w = c_2w/(c_2w + b1);
            //NB: alpha_x2y means: given s(M-1) is x what is the prob for adding y in t1
            // double alpha_r2r = a_2r/(a_2r + a_2w + c_1r);
            // double alpha_r2w = a_2w/(a_2r + a_2w + c_1r);
            // double alpha_w2r = a_2r/(a_2r + a_2w + c_1w);   
            // double alpha_w2w = a_2w/(a_2r + a_2w + c_1w);

            //probabilities following absorbtion matrix 
            double prob_mr_m_1r = (beta_2r * gamma_2r) / (1 - beta_1r * beta_2r); //s[M]=r --> s[M+1] = r
            double prob_2mr_m_1r = (alpha_1r * beta_1r)/(1 - beta_1r * beta_2r); //s[M]=r --> s[M-1] = r
            double prob_mr_m_1w = (beta_2w * gamma_2r) / (1 - beta_1r * beta_2w); //s[M]=r --> s[M+1] = w
            double prob_2mw_m_1r = (alpha_1w * beta_1r)/(1 - beta_1r * beta_2w); //s[M]= w --> s[M-1] = r
            double prob_mw_m_1w = (beta_2w * gamma_2w) / (1 - beta_1w * beta_2w); //s[M]=w --> s[M+1] = w
            double prob_2mw_m_1w = (alpha_1w * beta_1w)/(1 - beta_1w * beta_2w); //s[M]=w --> s[M-1] = w
            double prob_mw_m_1r = (beta_2r * gamma_2w) / (1 - beta_1w * beta_2r); //s[M]=w --> s[M+1] = r
            double prob_2mr_m_1w = (alpha_1r * beta_1w)/(1 - beta_1w * beta_2r); //s[M]=r --> s[M-1] = w
            string model_2 = "_abs_matrix_prob";

            // double prob_mr_m_1r = (b2 * c_2r) / (1 - b1 * b2); //s[M]=r --> s[M+1] = r
            // double prob_2mr_m_1r = (a_1r * b1)/(1 - b1 * b2); //s[M]=r --> s[M-1] = r
            // double prob_mr_m_1w = (b2 * c_2r) / (1 - b1 * b2); //s[M]=r --> s[M+1] = w
            // double prob_2mw_m_1r = (a_1w * b1)/(1 - b1 * b2); //s[M]= w --> s[M-1] = r
            // double prob_mw_m_1w = (b2* c_2w) / (1 - b1 * b2); //s[M]=w --> s[M+1] = w
            // double prob_2mw_m_1w = (a_1w * b1)/(1 - b1 * b2); //s[M]=w --> s[M-1] = w
            // double prob_mw_m_1r = (b2 * c_2w) / (1 - b1 * b2); //s[M]=w --> s[M+1] = r
            // double prob_2mr_m_1w = (a_1r * b1)/(1 - b1 * b2); //s[M]=r --> s[M-1] = w
            // string model_2 = "_abs_matrix_rates";

            //define rates for moving forward or backward
            double psi_rr_plus = a_2r * prob_mr_m_1r;
            double psi_rr_min = c_1r * prob_2mr_m_1r;
            double psi_rw_plus = a_2w * prob_mr_m_1w;
            double psi_rw_min = c_1r * prob_2mw_m_1r;
            double psi_ww_plus = a_2w * prob_mw_m_1w;
            double psi_ww_min = c_1w * prob_2mr_m_1w;
            double psi_wr_plus = a_2r * prob_mw_m_1r;
            double psi_wr_min = c_1w * prob_2mr_m_1w;
                  
            std::vector<int> mRNA_string;
            mRNA_string.push_back(0); //we begin with a right monomer
            int transition_state = 0;
            std::vector<structure_mRNA_move> mRNA_moves;
            int iteration = 0;
            string transition_state_string;
            int adding_monomer = 0;
            
            while(mRNA_string.size()< length_mRNA){
            // for(int j =0;j<20;j++){
                // Go to transition state 1r or 1w or 2, first MC step
                // create a random number betwen 0 and 1
                int M = mRNA_string.size() -1;


                double random_number1 = distribution(generator);      
                // if string is shorter than 2 monomers, always add a r or w
                if(mRNA_string.size() == 1){
                    if (mRNA_string[M] == 0){
                        double prob_rr = psi_rr_plus/(psi_rr_plus + psi_rw_plus);
                        if (random_number1 < prob_rr){
                            transition_state_string = "rr_plus";
                            mRNA_string.push_back(0);
                        }
                        else{
                            transition_state_string = "rw_plus";
                            mRNA_string.push_back(1);
                        }
                    }
                    else if (mRNA_string[M] == 1){
                        double prob_ww = psi_ww_plus/(psi_ww_plus + psi_wr_plus);
                        if (random_number1 < prob_ww){
                            transition_state_string = "ww_plus";
                            mRNA_string.push_back(1);
                        }
                        else{
                            transition_state_string = "wr_plus";
                            mRNA_string.push_back(0);
                        }
                    }
                }
                // if string longer than 1, you can go back 
                else{
                    // double random_number2 = distribution(generator);
                    // if s(M-1) = r
                    if(mRNA_string[M-1] == 0){
                        //if s[M] = r
                        if (mRNA_string[M] == 0){
                            //&rr is current string 
                            double prob_rr = psi_rr_plus/(psi_rr_plus + psi_rw_plus + psi_rr_min);
                            double prob_rw = psi_rw_plus/(psi_rr_plus + psi_rw_plus + psi_rr_min);
                            if (random_number1 < prob_rr){
                                transition_state_string = "rr_plus";
                                mRNA_string.push_back(0);
                            }
                            else if (random_number1 < (prob_rr + prob_rw)){
                                transition_state_string = "rw_plus";
                                mRNA_string.push_back(1);
                            }
                            else{
                                // monomer is removed from string
                                transition_state_string = "rr_min";
                                mRNA_string.pop_back();
                            }
                        }
                        //if s[M] = w
                        else if (mRNA_string[M] == 1){
                            //&rw is current string 
                            double prob_wr = psi_wr_plus/(psi_wr_plus + psi_ww_plus + psi_rw_min);
                            double prob_ww = psi_ww_plus/(psi_wr_plus + psi_ww_plus + psi_rw_min);
                            if (random_number1 < prob_wr){
                                transition_state_string = "wr_plus";
                                mRNA_string.push_back(0);
                            }
                            else if (random_number1 < (prob_wr + prob_ww)){
                                transition_state_string = "ww_plus";
                                mRNA_string.push_back(1);
                            }
                            else{
                                // monomer is removed from string
                                transition_state_string = "rw_min";
                                mRNA_string.pop_back();
                            }
                        }
                    }
                    // if s(M-1) = w
                    else if(mRNA_string[M-1] == 1){
                        //if s[M] = r
                        if (mRNA_string[M] == 0){
                            //&wr is current string 
                            double prob_rr = psi_rr_plus/(psi_rr_plus + psi_rw_plus + psi_wr_min);
                            double prob_rw = psi_rw_plus/(psi_rr_plus + psi_rw_plus + psi_wr_min);
                            if (random_number1 < prob_rr){
                                transition_state_string = "rr_plus";
                                mRNA_string.push_back(0);
                            }
                            else if (random_number1 < (prob_rr + prob_rw)){
                                transition_state_string = "rw_plus";
                                mRNA_string.push_back(1);
                            }
                            else{
                                // monomer is removed from string
                                transition_state_string = "wr_min";
                                mRNA_string.pop_back();
                            }
                        }
                        //if s[M] = w
                        else if (mRNA_string[M] == 1){
                            //&ww is current string 
                            double prob_wr = psi_wr_plus/(psi_wr_plus + psi_ww_plus + psi_ww_min);
                            double prob_ww = psi_ww_plus/(psi_wr_plus + psi_ww_plus + psi_ww_min);
                            if (random_number1 < prob_wr){
                                transition_state_string = "wr_plus";
                                mRNA_string.push_back(0);
                            }
                            else if (random_number1 < (prob_wr + prob_ww)){
                                transition_state_string = "ww_plus";
                                mRNA_string.push_back(1);
                            }
                            else{
                                // monomer is removed from string
                                transition_state_string = "ww_min";
                                mRNA_string.pop_back();
                            }
                        }
                    }
                }
            mRNA_moves.push_back(structure_mRNA_move());
            error_prob = count(mRNA_string.begin(), mRNA_string.end(), 1)/(double)mRNA_string.size();
            mRNA_moves[iteration].polymer = mRNA_string;
            mRNA_moves[iteration].length= mRNA_string.size();
            mRNA_moves[iteration].error = error_prob;
            mRNA_moves[iteration].transition_states = transition_state_string;
            iteration += 1;
            }

        // print the polymer computed
        // std::cout << "end polymer" << std::endl;
        // print_matrix(mRNA_string);
        // std::cout << "   " << std::endl;
        // std::cout << "   " << std::endl;

        // // // print path it took to come to this polymer (which states it passed)
        // std::cout << "path taken to get to this polymer:" << std::endl;
        // print_matrix_of_structures(mRNA_moves);

        //save the polymer and the path it took to get there
        int rounded_delta_g_pol = round(delta_g_pol);
        int rounded_delta_g_tt = round(delta_g_tt);
        std:cout << " done" << std::endl;
        string date_now = date(time(0));
        string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/221121output/";
        string filename_output = path + date_now + "model_jennys_way_" + model_1 + model_2 + "_delta_g_pol_" + std::to_string(rounded_delta_g_pol) + "_delta_g_tt_" + std::to_string(rounded_delta_g_tt) + ".csv";
        save_matrix_of_structures(mRNA_moves, filename_output);
        }
    }
}

