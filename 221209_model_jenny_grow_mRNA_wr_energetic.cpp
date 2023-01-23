/* 
21-11-22, Muriel Louman, AMOLF

Grow mRNA, assumingh template string of only ones. The energetic case of Jenny's paper. 
delta G_k = 0, delta G_TT = delta G. Testing th effect of different rates and different definition of the absorption 
matrix: consisting of the probabilities to move between the absorbing states and the transition 
states (defined by the rates) or consiting of rates to move between the states.

mRNA only consisting of zeros and ones.

model Jenny (not what she did in paper but what other students used):
Use psi_xy_plus or psi_xy_minus which define the rate to move forward(plus) or backward(min) between x and. using 
the rate t1 op t2 for the first step times the absorption matrix to move from there to the next absorbing state. In this
model you always make a move foward or backward.  
*/

#include <iostream>
#include <cmath> // math functions
#include <vector> // vector container
#include <string> // text strings
#include <fstream> // output to file
#include <sstream> 
#include <random> // random number generator
#include <bits/stdc++.h> //for count function

#include <boost/program_options.hpp>
namespace po = boost::program_options;


// #include <boost/variabels_def.hpp
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
    // std::vector<int> polymer;
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

void print_matrix_of_structures(std::vector<structure_mRNA_move> mRNA_movess){
    // loop through the array elements, which are all structures defined by structure_mRNA_move
    for(int i=0; i< mRNA_movess.size(); ++i){
        // cout << "creates new polymer" << endl;
        // print_matrix(mRNA_movess[i].polymer);
        cout << "adds or removed monomer" << endl;
        cout << mRNA_movess[i].monomer << endl;
        cout << "   " << endl;
        cout << "transition state used" << endl;
        cout << mRNA_movess[i].transition_state << endl;
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
    fss << "monomer" << "," << "transition state used" << "," << "length polymer" <<"," << "error probability"  << std::endl;
    // loop through the array elements, which are all structures defined by structure_mRNA_move
    for(int i=0; i< mRNA_movess.size(); ++i){
        fss << mRNA_movess[i].monomer << "," << mRNA_movess[i].transition_state ;
        fss << "," << mRNA_movess[i].length <<"," << mRNA_movess[i].error << std ::endl;
    }
    fss.close();
}

string declare_filename(bool all_steps, string DIR, string model_1, int delta_g_pol, int delta_g_tt, int step){
    string date_now = date(time(0));
    int rounded_delta_g_pol = round(delta_g_pol);
    int rounded_delta_g_tt = round(delta_g_tt);
    string filename_output = "";
    if(all_steps ==true){
        string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/" + DIR + "/";
        filename_output = path + date_now + "model_Jenny_" + model_1 + "_delta_g_pol_" + std::to_string(rounded_delta_g_pol) + "_delta_g_tt_" + std::to_string(rounded_delta_g_tt) + "_allsteps.csv";
    }
    else if(all_steps == false){
        int step_rounded = round(step);
        string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/" + DIR + "/";
        filename_output = path + date_now + "model_Jenny_step" + std::to_string(step_rounded) + "_" + model_1 + "_delta_g_pol_" + std::to_string(rounded_delta_g_pol) + "_delta_g_tt_" + std::to_string(rounded_delta_g_tt) + ".csv";
    }
    return filename_output;
}


int main(int ac, char* av[]){
    //in main
    // bool multiple_DNA = false;
    // string DIR = "230110output";
    // string x_def = "1"; //"exp"
    // string M_model = "M1"; //"M2"
    // int number_steps = 0;
    // bool variables_in_main = false;

    bool multiple_DNA = true;
    string DIR = "";
    int number_steps = 0;

    bool variables_in_main = true;
    if(variables_in_main == true){
        po::options_description desc("Allowed options");
        try {

            po::options_description desc("Allowed options");
            desc.add_options()
                ("help", "produce help message")
                ("multiple_DNA", po::value<bool>(), "compute multiple DNA string") //set declaration of values
                ("DIR", po::value<string>(), "output directory name") 
                ("number_DNA", po::value<int>(), "number of DNA computed")
            ;
            po::variables_map vm;
            po::store(po::parse_command_line(ac, av, desc), vm);
            po::notify(vm);



            if (vm.count("help")) {
                std::cout << desc << "\n";
                return 0;}
            if (vm.count("multiple_DNA")) {
                std::cout << "we are computing multiple DNA strings: "
                    << vm["multiple_DNA"].as<bool>() << ".\n"; //print value definition if set
                multiple_DNA = vm["multiple_DNA"].as<bool>();
                } 
            if (vm.count("DIR")) {
                std::cout << "output file is set to: "
                    << vm["DIR"].as<string>() << ".\n"; //print value definition if set
                DIR = vm["DIR"].as<string>();
                } 
            if (vm.count("number_DNA")) {
                std::cout << "number of DNA computed: "
                    << vm["number_DNA"].as<int>() << ".\n"; //print value definition if set
                number_steps = vm["number_DNA"].as<int>();
                }
            else {
                std::cout << "Not all variables are set \n";}
        }
        catch(exception& e) {
            cerr << "error: " << e.what() << "\n";
            return 1;}
        catch(...) {
            cerr << "Exception of unknown type!\n";}
    }




    fstream fss;


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
    int x = 0;
    int step = 0;

    double L_delta_g_pol[] = {0,1,2,3,4,5,6,7,8,9,10}; //we would like to have more but takes computer time
    double L_delta_g_tt[] = {0,2,4,6,8,10};
    int length_mRNA = 3000;

    // random number generator for uniform distr between 0 and 1
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
       
    

    // define template string
    // std::vector<int> mRNA_template = template_string(length_mRNA);

    //loop over different g_tt and g_bb
    for (double delta_g_pol : L_delta_g_pol) {
        for (double delta_g_tt : L_delta_g_tt) {
            
            //define rates for moving forward or backward
            double psi_rr_plus, psi_rw_plus, psi_ww_plus, psi_wr_plus; 
            psi_rr_plus = psi_rw_plus = psi_ww_plus = psi_wr_plus = 1;
            double psi_rr_min = exp(-delta_g_pol);
            double psi_rw_min = exp(-delta_g_pol+delta_g_tt);
            double psi_ww_min = exp(-delta_g_pol);
            double psi_wr_min = exp(-delta_g_pol-delta_g_tt);
            string model_1 ="full_steps";


            
            
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
                string transition_state_string;
                int adding_monomer = 0;
                string monomer_status = "0";

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
                                monomer_status = "0";
                            }
                            else{
                                transition_state_string = "rw_plus";
                                mRNA_string.push_back(1);
                                monomer_status = "1";
                            }
                        }
                        else if (mRNA_string[M] == 1){
                            double prob_ww = psi_ww_plus/(psi_ww_plus + psi_wr_plus);
                            if (random_number1 < prob_ww){
                                transition_state_string = "ww_plus";
                                mRNA_string.push_back(1);
                                monomer_status = "1";
                            }
                            else{
                                transition_state_string = "wr_plus";
                                mRNA_string.push_back(0);
                                monomer_status = "0";
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
                                    monomer_status = "0";
                                }
                                else if (random_number1 < (prob_rr + prob_rw)){
                                    transition_state_string = "rw_plus";
                                    mRNA_string.push_back(1);
                                    monomer_status = "1";
                                }
                                else{
                                    // monomer is removed from string
                                    transition_state_string = "rr_min";
                                    mRNA_string.pop_back();
                                    monomer_status = "-";
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
                                    monomer_status = "0";
                                }
                                else if (random_number1 < (prob_wr + prob_ww)){
                                    transition_state_string = "ww_plus";
                                    mRNA_string.push_back(1);
                                    monomer_status = "1";
                                }
                                else{
                                    // monomer is removed from string
                                    transition_state_string = "rw_min";
                                    mRNA_string.pop_back();
                                    monomer_status = "-";
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
                                    monomer_status = "0";
                                }
                                else if (random_number1 < (prob_rr + prob_rw)){
                                    transition_state_string = "rw_plus";
                                    mRNA_string.push_back(1);
                                    monomer_status = "1";
                                }
                                else{
                                    // monomer is removed from string
                                    transition_state_string = "wr_min";
                                    mRNA_string.pop_back();
                                    monomer_status = "-";
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
                                    monomer_status = "0";
                                }
                                else if (random_number1 < (prob_wr + prob_ww)){
                                    transition_state_string = "ww_plus";
                                    mRNA_string.push_back(1);
                                    monomer_status = "1";
                                }
                                else{
                                    // monomer is removed from string
                                    transition_state_string = "ww_min";
                                    mRNA_string.pop_back();
                                    monomer_status = "-";
                                }
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

