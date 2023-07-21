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

#include <boost/program_options.hpp>
namespace po = boost::program_options;


string date(time_t now){
    tm *ltm = localtime(&now);
    string datum = std::to_string(1900 + ltm->tm_year) + std::to_string(1 + ltm->tm_mon) + std::to_string(ltm->tm_mday);
    return datum;
}

//Create a structure variable called mRNA_move
struct structure_mRNA_move{
    string monomer;
    string transition_state;
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

void print_matrix_double(std::vector<double> input){
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

string declare_filename(bool all_steps, string DIR, string model_1, double delta_g_pol, double delta_g_tt, int step){
    string date_now = date(time(0));
    int rounded_delta_g_pol = round(delta_g_pol);
    string str_rounded_delta_g_pol = std::to_string(rounded_delta_g_pol);
    if(delta_g_pol == -log(2)){
        str_rounded_delta_g_pol = "-ln2";
    }
    int rounded_delta_g_tt = round(delta_g_tt);
    string filename_output = "";
    if(all_steps ==true){
        string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/" + DIR + "/";
        filename_output = path + date_now + "model_muriel_" + model_1 + "_delta_g_pol_" + str_rounded_delta_g_pol + "_delta_g_tt_" + std::to_string(rounded_delta_g_tt) + "_allsteps.csv";
    }
    else if(all_steps == false){
        int step_rounded = round(step);
        string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/" + DIR + "/";
        filename_output = path + date_now + "model_muriel_step" +  std::to_string(step_rounded) + "_" + model_1 + "_delta_g_pol_" + str_rounded_delta_g_pol + "_delta_g_tt_" + std::to_string(rounded_delta_g_tt) + ".csv";
    }
    return filename_output;
}

string declare_filename_parameter(bool all_steps, string DIR, string M_model, string x_def){
    string date_now = date(time(0));
    string filename_output = "";
    string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/" + DIR + "/";
    if(all_steps ==true){
        filename_output = path + date_now + "parameter_file_model_muriel_" + M_model + "_x" + x_def + "_allsteps.csv";}
    else{
        filename_output = path + date_now + "parameter_file_model_muriel_" + M_model + "_x" + x_def + "_seperate_steps.csv";}
     
    return filename_output;
}

int main(){
    //in main
    bool all_steps = true;
    bool multiple_DNA = true;
    string DIR = "test_output";
    string x_def = "exp";
    string M_model = "M1"; //"M2"
    int number_steps = 100;
    int length_mRNA = 100;
    vector<double> L_delta_g_tt = {0,2,4,6,8,10};
    vector<double> L_delta_g_pol = {0,1,2,3,4,5,6,7,8,9,10};


    
fstream fss;
    string parameter_output_file = declare_filename_parameter(multiple_DNA, DIR, M_model, x_def);
    fss.open(parameter_output_file.c_str(), ios::out | ios::app);
    fss << "--multiple_DNA=" << multiple_DNA << ".\n" ;
    fss << "--DIR=" << DIR << ".\n";
    fss << "--number_DNA=" << number_steps << ".\n";
    fss << "--x_def=" << x_def << ".\n";
    fss << "--M_model=" << M_model << ".\n";
    fss << "--mRNA_length=" << length_mRNA << ".\n";
    fss << "--L_g_pol=";
    std::vector<double> in = L_delta_g_pol;
    // loop through the array elements
    for(int j=0; j< in.size(); ++j){
        fss << in.at(j);
        fss << " ";}
    fss <<".\n"  ;
    fss << "--L_g_tt=";
    std::vector<double> inn = L_delta_g_tt;
    // loop through the array elements
    for(int j=0; j< inn.size(); ++j){
        fss << inn.at(j);
        fss << " ";}       
    fss << ".\n";
    fss << "\n";
    fss << "model used: Muriel.\n";   
    fss.close();

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
    double x = 0;
    int x1= 0;
    int step = 0;
    string x_string = "";

    //for making the error landscape defined by different delta G_tt and delta G_pol
    // double L_delta_g_tt[] = {0,1,2,3,4,6,8,10};//{0,2,4,6,8,10}; //{2,4,6,8,10};
    // double L_delta_g_pol[] = {0,1,2,3,4,5,6,7,8,9,10};
    // int length_mRNA = 500;

    
    // random number generator for uniform distr between 0 and 1
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    
    for (double delta_g_pol : L_delta_g_pol){
        if(delta_g_pol == 3000){
                delta_g_pol = -log(2);
        }
        for (double delta_g_tt : L_delta_g_tt){
            // define rates 
            if(x_def =="exp"){
                // x = exp(delta_g_tt + 4);
                // x = exp(delta_g_tt + 20);
                // x1 = exp(delta_g_tt + 20);
                x = exp(delta_g_tt + 5);
                std::cout << x << std::endl;
                std::cout << x1 << std::endl;
                // x = x1;
                
                x_string = "delta_g_tt_pow2";}
            else if(x_def == "1"){
                x = 1;
                int rounded_x = round(x);
                x_string = std::to_string(rounded_x);}

            int y = 1;         
            int rounded_y = round(y);
            
            if (M_model == "M1"){
            a_2r = a_1r = a_2w = 1;
            a_1w = exp(delta_g_tt);
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


            //define probabilities 
            //if s[M-1] = r
            double alpha_r2r = a_2r/(a_2w + a_2r + c_1r);
            double alpha_r2w = a_2w/(a_2w + a_2r + c_1r);
            double gamma_1r =  c_1r/(a_2w + a_2r + c_1r);
            //if s[M-1] = w
            double alpha_w2r = a_2r/(a_2r + a_2w + c_1w);
            double alpha_w2w = a_2w/(a_2r + a_2w + c_1w);
            double gamma_1w =  c_1w/(a_2r + a_2w + c_1w);
            //if s[M+1] = r
            double alpha_1r = a_1r/(a_1r +b2);
            double beta_2r = b2/(a_1r +b2);
            //if s[M+1] = w
            double alpha_1w = a_1w/(a_1w +b2);
            double beta_2w = b2/(a_1w + b2);
            //if s[M] = r
            double beta_1r = b1/(c_2r + b1);
            double gamma_2r = c_2r/(c_2r + b1);
            //if s[M] = w
            double beta_1w = b1/(c_2w + b1);
            double gamma_2w = c_2w/(c_2w + b1);
            
            
            //for saving multiple steps end error 
            if(multiple_DNA == true){
                string filename = declare_filename(multiple_DNA, DIR, model_1, delta_g_pol, delta_g_tt, step);
                fss.open(filename.c_str(), ios::out | ios::app);
                // write the file headers
                fss << "mRNA polymer" << "," << "error probability"  << "\n";
            }


            for (step = 0; step <= number_steps; step += 1) {
                // for saving all moves to make the polymer 
                if(multiple_DNA == false){
                    string filename = declare_filename(multiple_DNA, DIR, model_1, delta_g_pol, delta_g_tt, step);
                    fss.open(filename.c_str(), ios::out | ios::app);
                    // write the file headers
                    fss << "monomer added removed" << "," << "length polymer" << "," << "transition state used" <<"," << "error probability"  << "\n";
                }
                std::vector<int> mRNA_string;
                mRNA_string.push_back(0); //we begin with a right monomer
                int transition_state = 0;
                structure_mRNA_move mRNA_moves;
                int iteration = 0;
                string transition_state_string = "t";
                string monomer_status = "0";  
        

            


                while (mRNA_string.size()<length_mRNA){
                    int M = mRNA_string.size() -1;

                    //define to which treansition state we go
                    //if string size = 1 you can only go forward
                    if(mRNA_string.size() == 1){
                        double random_number1 = distribution(generator);
                        if (random_number1 < 0.5){
                            //go to t_1r
                            transition_state = 1;
                            transition_state_string = "t1r";
                        }
                        else{
                            //go to t_1w
                            transition_state = -1;
                            transition_state_string = "t1w";
                        }
                    }
                    else{
                        double random_number2 = distribution(generator);
                        // if s[M-1] = r
                        if (mRNA_string[M-1] ==0){
                            if (random_number2 < alpha_r2r){
                                transition_state = 1;
                                transition_state_string = "t1r";
                            }
                            else if (random_number2 < (alpha_r2r + alpha_r2w)){
                                transition_state = -1;
                                transition_state_string = "t1w";
                            }
                            else{
                                transition_state = 2;
                                transition_state_string = "t2r";
                            }   
                        }
                        // if s[M-1] = w
                        double random_number7 = distribution(generator);
                        if (mRNA_string[M-1] ==1){
                            if (random_number7 < alpha_w2r){
                                transition_state = 1;
                                transition_state_string = "t1r";
                            }
                            else if (random_number7 < (alpha_w2r + alpha_w2w)){
                                transition_state = -1;
                                transition_state_string = "t1w";
                            }
                            else{
                                transition_state = -2;
                                transition_state_string = "t2w";
                            }   
                        }
                    }



                    // define from the transition state we are in, which monomer is added, not added or removed.
                    // if in trabnsition state t1r
                    if (transition_state == 1){
                        double random_number3 = distribution(generator);
                        // if s[M] = r
                        if (mRNA_string[M] == 0){
                            double prob_r_t1r = (beta_2r * gamma_2r)/(1- beta_1r * beta_2r);
                            if (random_number3 < prob_r_t1r){
                                //add r
                                mRNA_string.push_back(0);
                                monomer_status = "0";  
                            }
                            else{
                                // else nothing happens
                                monomer_status = "n";
                            }
                            
                        }
                        // if s[M] = w
                        else if (mRNA_string[M] ==1){
                            double prob_w_1tr = (beta_2r * gamma_2w)/(1 - beta_1w * beta_2r);
                            if (random_number3 < prob_w_1tr){
                                mRNA_string.push_back(0);
                                monomer_status = "0";
                            }
                            else{
                                // else nothing happens --> dont add a monomer
                                monomer_status = "n";
                            }
                            
                        }
                    }
                    // if in transition state 1tw
                    else if (transition_state == -1){
                        double random_number4 = distribution(generator);
                        if (mRNA_string[M] == 0){
                            double prob_r_1tw = (beta_2w *gamma_2r)/(1 - beta_1r * beta_2w);
                            if (random_number4 < prob_r_1tw){
                                mRNA_string.push_back(1);
                                monomer_status = "1";
                            }
                            else{
                                // else nothing happens --> dont add a monomer
                                monomer_status = "n";
                            }
                        }
                        else if (mRNA_string[M] == 1){
                            double prob_w_1tw = (beta_2w * gamma_2w)/(1 - beta_1w * beta_2w);
                            if (random_number4 < prob_w_1tw){
                                mRNA_string.push_back(1);
                                monomer_status = "1";
                            } 
                            else{
                                // else nothing happens --> dont add a monomer
                                monomer_status = "n";
                            }
                        }
                    }
                    //if in transition state 2tr
                    else if (transition_state == 2){
                        double random_number5 = distribution(generator);
                        if (mRNA_string[M] == 0){
                            double prob_r_2tr = (alpha_1r * beta_1r)/(1 - beta_1r * beta_2r);
                            if (random_number5 < prob_r_2tr){
                                mRNA_string.pop_back();
                                monomer_status = "-";
                            }
                            else{
                                // else nothing happens --> dont add a monomer
                                monomer_status = "n";
                            }
                        }
                        else if (mRNA_string[M] == 1){
                            double prob_w_2tr = (alpha_1w * beta_1r)/(1 - beta_1r * beta_2w);
                            if (random_number5 < prob_w_2tr){
                                mRNA_string.pop_back();
                                monomer_status = "-";
                            }
                            else{
                                // else nothing happens --> dont add a monomer
                                monomer_status = "n";
                            }
                        }
                    }
                    // if in transition state 2tw
                    else if (transition_state == -2){
                        double random_number6 = distribution(generator);
                        if (mRNA_string[M] == 0){
                            double prob_r_2tw = (alpha_1r * beta_1w) / (1 - beta_1w * beta_2r);
                            if (random_number6 < prob_r_2tw){
                                mRNA_string.pop_back();
                                monomer_status = "-";
                            } 
                            else{
                                // else nothing happens --> dont add a monomer
                                monomer_status = "n";
                            }
                        }
                        else if (mRNA_string[M] == 1){
                            double prob_w_2tw = (alpha_1w * beta_1w)/(1 - beta_1w * beta_2w);
                            if (random_number6 < prob_w_2tw){
                                mRNA_string.pop_back();
                                monomer_status = "-";
                            }
                            else{
                                // else nothing happens --> dont add a monomer
                                monomer_status = "n";
                            }
                        }
                    }
                
                error_prob = count(mRNA_string.begin(), mRNA_string.end(), 1)/(double)mRNA_string.size();
                mRNA_moves.monomer = monomer_status;
                mRNA_moves.length= mRNA_string.size();
                mRNA_moves.error = error_prob;
                mRNA_moves.transition_state = transition_state_string;

                                //save information for every iteration to CSV file  
                if(multiple_DNA == false){         
                    fss << mRNA_moves.monomer << "," << mRNA_moves.length <<"," << mRNA_moves.transition_state << ","  << mRNA_moves.error;
                    fss << "\n";
                    }
                iteration += 1;
                }
            //save information for every iteration to CSV file 
            if(multiple_DNA == false){
                std::vector<int> in = mRNA_string;
                // loop through the array elements
                for(int j=0; j< in.size(); ++j){
                    fss << in.at(j) ;
                    }
                fss << std ::endl;
                fss.close();
            }

            //save information for one step, only error info
            if(multiple_DNA == true){
                std::vector<int> in = mRNA_string;
                // loop through the array elements
                for(int j=0; j< in.size(); ++j){
                    fss << in.at(j) ;
                    }
                fss << "," << mRNA_moves.error << "\n";
            }
            }
        //save information for all steps, error info
        if(multiple_DNA == true){
            fss << std ::endl;
            fss.close();
        }
        }
    } 
}
