/* 
11-11-22, Muriel Louman, AMOLF

Grow mRNA, assumingh template string of only ones. The energetic case of Jenny's paper. 
delta G_k = 0, delta G_TT = delta G. Testing th effect of different rates.

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


string date(time_t now){
    tm *ltm = localtime(&now);
    string datum = std::to_string(1900 + ltm->tm_year) + std::to_string(1 + ltm->tm_mon) + std::to_string(ltm->tm_mday);
    return datum;
}

//Create a structure variable called mRNA_move
struct structure_mRNA_move{
    // std::vector<int> polymer;
    string monomer;
    int length;
    double error;
    string transition_state;
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
        // cout << "creates new polymer" << endl;
        // print_matrix(mRNA_movess[i].polymer);
        cout << "monomer added or removed" << endl;
        cout << mRNA_movess[i].monomer << endl;
        cout << "   " << endl;
        cout << "of length" << endl;
        cout << mRNA_movess[i].length << endl;
        cout << "with error probability"  << endl;
        cout << mRNA_movess[i].error << endl;
        cout << "moves to transition state" << endl;
        cout << mRNA_movess[i].transition_state << endl;
        }
}

double alpha_2(double a_2r, double a_2w, double c_1, double a_2){
    return a_2/(a_2r + a_2w + c_1);
}

double alpha_1(double b_2, double a_1){
    return a_1/(a_1 + b_2);
}

double beta_2(double b_2, double a_1){
    return b_2/(b_2 + a_1);
}

double beta_1(double b_1, double c_2){
    return b_1/(b_1 + c_2);
}

double gamma_2(double b_1, double c_2){
    return c_2/(b_1 + c_2);
}

double gamma_1(double a_2r, double a_2w, double c_1){
    return c_1/(a_2r + a_2w + c_1);
}

string declare_filename(bool all_steps, string DIR, string model_1, int delta_g_pol, int delta_g_tt, int step){
    string date_now = date(time(0));
    int rounded_delta_g_pol = round(delta_g_pol);
    int rounded_delta_g_tt = round(delta_g_tt);
    string filename_output = "";
    if(all_steps ==true){
        string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/"+ DIR + "/";
        filename_output = path + date_now + "model_muriel_big_abs_" + model_1 + "_delta_g_pol_" + std::to_string(rounded_delta_g_pol) + "_delta_g_tt_" + std::to_string(rounded_delta_g_tt) + "_allsteps.csv";
    }
    else if(all_steps == false){
        int step_rounded = round(step);
        string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/" + DIR + "/";
        filename_output = path + date_now + "model_muriel_big_abs_step" + std::to_string(step_rounded) + "_" + model_1 + "_delta_g_pol_" + std::to_string(rounded_delta_g_pol) + "_delta_g_tt_" + std::to_string(rounded_delta_g_tt) + ".csv";            
    }
    return filename_output;
}

int main(){
    //in main
    bool all_steps = true;
    string DIR = "230110output";
    string x_def = "1"; //"exp"
    string M_model = "M1"; //"M2"

    // define variables
    double k_on = 1.0;
    double k_bb = 1;
    double r_con = 1.0;
    double w_con = 1.0;
    double r_con_star = 1.5;
    double w_con_star = 1.5;

    double a_2r, a_1r, a_2w, a_1w, b2, b1, c_2r, c_1r, c_2w, c_1w;
    a_2r = a_1r = a_2w = a_1w = b2 = b1 = c_2r = c_1r = c_2w = c_1w = 0;
    string model_1 = "0";
    double error_prob = 0;
    int number_steps = 20;
    int x = 0;
    int step = 0;

    // random number generator for uniform distr between 0 and 1
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    //for making the error landscape defined by different delta G_tt and delta G_pol
    double L_delta_g_tt[] = {6};//{0,2,4,6,8,10};
    double L_delta_g_pol[] = {10};//{0,1,2,3,4,5,6,7,8,9,10};
    int length_mRNA = 300;
       
    
    //loop over different g_tt and g_bb
    for (double delta_g_pol : L_delta_g_pol){
        for (double delta_g_tt : L_delta_g_tt){
            
            if(x_def =="exp"){
                x = exp(pow(delta_g_tt,2));
                string x_string = "delta_g_tt_pow2";}
            else if(x_def == "1"){
                x = 1;}

            int y = 1;         
            int rounded_x = round(x);
            int rounded_y = round(y);
            string x_string = std::to_string(rounded_x);
            
            if (M_model == "M1"){
            a_2r = a_1r = a_2w = 1;
            b2 = x;
            b1 = x * exp(-delta_g_pol);
            c_2r = c_1r = c_2w = y;
            c_1w = y * exp(-delta_g_tt);
            model_1 = "M1_x" + x_string;}
            
            if (M_model == "M2"){
            a_2r = a_1r = a_2w = 1;
            a_1w = exp(delta_g_tt);
            b2 = x;
            b1 = x * exp(-delta_g_pol);
            c_2r = c_1r = c_1w = y;
            c_2w = y * exp(delta_g_tt);
            model_1 = "M2_x" + x_string;}

            double A, B, C, D, E, F, G, H, I, J, L, M, N, O, P;
            A, B, C, D, E, F, G, H, I, J, L, M, N, O, P = 0;
            
            fstream fss;
            
            
            // for saving multiple steps end error 
            if(all_steps == true){
                string filename = declare_filename(all_steps, DIR, model_1, delta_g_pol, delta_g_tt, step);
                fss.open(filename.c_str(), ios::out | ios::app);
                // write the file headers
                fss << "mRNA polymer" << "," << "error probability"  << "\n";
                    }
    
            for (step = 0; step <= number_steps; step += 1) {
                // for saving all moves to make the polymer 
                if(all_steps == false){
                    string filename = declare_filename(all_steps, DIR, model_1, delta_g_pol, delta_g_tt, step);
                    fss.open(filename.c_str(), ios::out | ios::app);
                    // write the file headers
                    fss << "monomer added removed" << "," << "length polymer" << "," << "transition state used" <<"," << "error probability"  << "\n";
                }

                std::vector<int> mRNA_string;
                mRNA_string.push_back(0); 
                mRNA_string.push_back(0); //we begin with two right monomers
                int transition_state = 0;
                structure_mRNA_move mRNA_moves;
                int iteration = 0;
                string transition_state_string = "t";
                string monomer_status = "0";   

                while (mRNA_string.size()<length_mRNA){
                    // Go to transition state 1r or 1w or 2, first MC step
                    // create a random number betwen 0 and 1
                    int Q = mRNA_string.size() -1;

                    
                    // define the variables used in the absorption matrix
                    if (mRNA_string[Q] == 0){
                        G = alpha_1(b2, a_1r);
                        H = beta_2(b2, a_1r);
                        I = beta_1(b1, c_2r);
                        J = gamma_2(b1, c_2r);
                        M = alpha_1(b2, a_1w);
                        N = beta_2(b2, a_1w);
                        O = beta_1(b1, c_2r);
                        P = gamma_2(b1, c_2r);
                        if(mRNA_string[Q-1] == 0){
                            A = alpha_1(b2,a_1r);
                            B = beta_2(b2, a_1r);
                            C = beta_1(b1, c_2r);
                            D = gamma_2(b1, c_2r);
                            E = gamma_1(a_2r, a_2w, c_1r);
                            F = alpha_2(a_2r,a_2w, c_1r,a_2r);
                            L = alpha_2(a_2r,a_2w, c_1r,a_2w);
                        }
                        else if(mRNA_string[Q-1] == 1){
                            A = alpha_1(b2, a_1r);
                            B = beta_2(b2, a_1r);
                            C = beta_1(b1, c_2w);
                            D = gamma_2(b1, c_2w);
                            E = gamma_1(a_2r, a_2w, c_1w);
                            F = alpha_2(a_2r, a_2w, c_1w, a_2r);
                            L = alpha_2(a_2r, a_2w, c_1w, a_2w);
                        }
                    }

                    else if (mRNA_string[Q] == 1){
                        G = alpha_1(b2, a_1r);
                        H = beta_2(b2, a_1r);
                        I = beta_1(b1, c_2w);
                        J = gamma_2(b1, c_2w);
                        M = alpha_1(b2, a_1w);
                        N = beta_2(b2, a_1w);
                        O = beta_1(b1, c_2w);
                        P = gamma_2(b1, c_2w);
                        if(mRNA_string[Q-1] == 0){
                            A = alpha_1(b2, a_1w);
                            B = beta_2(b2, a_1w);
                            C = beta_1(b1, c_2r);
                            D = gamma_2(b1, c_2r);
                            E = gamma_1(a_2r, a_2w, c_1r);
                            F = alpha_2(a_2r, a_2w, c_1r, a_2r);
                            L = alpha_2(a_2r, a_2w, c_1r, a_2w);
                        }
                        else if(mRNA_string[Q-1] == 1){
                            A = alpha_1(b2, a_1w);
                            B = beta_2(b2, a_1w);
                            C = beta_1(b1, c_2w);
                            D = gamma_2(b1, c_2w);
                            E = gamma_1(a_2r, a_2w, c_1w);
                            F = alpha_2(a_2r, a_2w, c_1w, a_2r);
                            L = alpha_2(a_2r, a_2w, c_1w, a_2w);
                        }
                    }

                    double Z = (1 - B*C - D*E - F*G + B*C*F*G - H*I + B*C*H*I + D*E*H*I - L*M + B*C*L*M + H*I*L*M - B*C*H*I*L*M - N*O + B*C*N*O + D*E*N*O + F*G*N*O - B*C*F*G*N*O + H*I*N*O - B*C*H*I*N*O - D*E*H*I*N*O);
                    double backward = (A*(C*E -C*E*H*I - C*E*N*O + C*E*H*I*N*O))/Z;
                    double forward_r = (J*(F*H - B*C*F*H - F*H*N*O + B*C*F*H*N*O))/Z;
                    double forward_w = (P*(L*N - B*C*L*N - H*I*L*N + B*C*H*I*L*N))/Z;
                    double random_number1 = distribution(generator);   
                    if (mRNA_string.size() ==2){
                        double P_forward_r = forward_r /(forward_r + forward_w);
                        double P_forward_w = forward_w /(forward_r + forward_w);
                        if (random_number1 < P_forward_r){
                                mRNA_string.push_back(0);
                                monomer_status = "0";
                            }
                        else{
                                mRNA_string.push_back(1);
                                monomer_status = "1";
                            }
                    }
                    else{
                        double P_backward = backward/(forward_r + forward_w + backward);
                        double P_forward_r = forward_r /(forward_r + forward_w + backward);
                        double P_forward_w = forward_w /(forward_r + forward_w + backward);
                        if (random_number1 < P_forward_r){
                                mRNA_string.push_back(0);
                                monomer_status = "0";
                            }
                        else if(random_number1 < (P_forward_r + P_forward_w)){
                                mRNA_string.push_back(1);
                                monomer_status = "1";
                            }
                        else{
                            mRNA_string.pop_back();
                            monomer_status = "-";
                        }
                    }
                
                error_prob = count(mRNA_string.begin(), mRNA_string.end(), 1)/(double)mRNA_string.size();
                mRNA_moves.monomer = monomer_status;
                mRNA_moves.length= mRNA_string.size();
                mRNA_moves.error = error_prob;
                mRNA_moves.transition_state = transition_state_string;

                //save information for every iteration to CSV file  
                if(all_steps == false){         
                    fss << mRNA_moves.monomer << "," << mRNA_moves.length <<"," << mRNA_moves.transition_state << ","  << mRNA_moves.error;
                    fss << "\n";
                    }

                iteration += 1;
                }

            //save information for every iteration to CSV file 
            if(all_steps == false){
                std::vector<int> in = mRNA_string;
                // loop through the array elements
                for(int j=0; j< in.size(); ++j){
                    fss << in.at(j) ;
                    }
                fss << std ::endl;
                fss.close();
            }

            //save information for one step, only error info
            if(all_steps == true){
                std::vector<int> in = mRNA_string;
                // loop through the array elements
                for(int j=0; j< in.size(); ++j){
                    fss << in.at(j) ;
                    }
                fss << "," << mRNA_moves.error << "\n";
            }
            
            }
        //save information for all steps, error info
        if(all_steps == true){
            fss << std ::endl;
            fss.close();
        }
        }
    }
}