/* 
11-11-22, Muriel Louman, AMOLF

Grow mRNA, assumingh template string of only ones. The energetic case of Jenny's paper. 
delta G_k = 0, delta G_TT = delta G. Testing th effect of different rates.

mRNA only consisting of zeros and ones.
No_abs_matrix model:
Using no absorption matrix such that every step is defined by the probability to move between the next state or 
transition state 

Interstingly this gives the same results for the different defenition of the rates as when using model muriel
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
        cout << mRNA_movess[i].transition_states << endl;
        }
}



string declare_filename(bool all_steps, string DIR, string model_1, double delta_g_pol, double delta_g_tt, int step, string length){
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
        filename_output = path + date_now + "model_no_abs_matrix_" + model_1 + "_delta_g_pol_" + str_rounded_delta_g_pol + "_delta_g_tt_" + std::to_string(rounded_delta_g_tt) + "length_copy" + length +"_allsteps.csv";
    }
    else if(all_steps == false){
        int step_rounded = round(step);
        string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/" + DIR + "/";
        filename_output = path + date_now + "model_no_abs_matrix_step" +  std::to_string(step_rounded) + "_" + model_1 + "_delta_g_pol_" + str_rounded_delta_g_pol + "_delta_g_tt_" + std::to_string(rounded_delta_g_tt) + ".csv";
    }
    return filename_output;
}

string declare_filename_parameter(bool all_steps, string DIR, string M_model, string x_def, string length){
    string date_now = date(time(0));
    string filename_output = "";
    string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/" + DIR + "/";
    if(all_steps ==true){
        filename_output = path + date_now + "parameter_file_model_no_abs_" + M_model + "_x" + x_def + "length_copy" + length +"_allsteps.csv";}
    else{
        filename_output = path + date_now + "parameter_file_model_no_abs_" + M_model + "_x" + x_def + "_seperate_steps.csv";}
     
    return filename_output;
}


int main(int ac, char* av[]){
    // in main
    bool multiple_DNA=true;
    string DIR ="230714output";
    int number_steps =30;
    string x_def= "exp";
    string M_model= "MB2";
    int length_mRNA=300;
    vector<double> L_delta_g_pol = {0,1,2,3,4,5,6,7,8,9,10};
    vector<double> L_delta_g_tt={3};
    bool variables_in_main = false;

    // bool multiple_DNA = true;
    // string DIR = "";
    // int number_steps = 0;
    // string x_def = "exp";
    // string M_model = "M2";
    // int length_mRNA = 1;
    // vector<double> L_delta_g_tt;
    // vector<double> L_delta_g_pol;


    // bool variables_in_main = true;
    // if(variables_in_main == true){
    //     po::options_description desc("Allowed options");
    //     try {

    //         po::options_description desc("Allowed options");
    //         desc.add_options()
    //             ("help", "produce help message")
    //             ("multiple_DNA", po::value<bool>(), "compute multiple DNA string") //set declaration of values
    //             ("DIR", po::value<string>(), "output directory name") 
    //             ("number_DNA", po::value<int>(), "number of DNA computed")
    //             ("x_def", po::value<string>(), "strength x, backbone") 
    //             ("M_model", po::value<string>(), "rate definition M1(backward correction) or M2(forward correction)") 
    //             ("mRNA_length", po::value<int>(), "length of the mRNA string formed")
    //             ("L_g_pol", po::value< vector<double> >(), "delta G_pol definitions list") 
    //             ("L_g_tt", po::value< vector<double> >(), "delta G_tt definitions list") 
    //         ;
    //         po::variables_map vm;
    //         po::store(po::parse_command_line(ac, av, desc), vm);
    //         po::notify(vm);

    //         std::cout << "--multiple_DNA=" << vm["multiple_DNA"].as<bool>() << ".\n" ;
    //         std::cout << "--DIR=" << vm["DIR"].as<string>() << ".\n";
    //         std::cout << "--number_DNA=" << vm["number_DNA"].as<int>() << ".\n";
    //         std::cout << "--x_def=" << vm["x_def"].as<string>() << ".\n";
    //         std::cout << "--M_model=" << vm["M_model"].as<string>() << ".\n";
    //         std::cout << "--mRNA_length=" << vm["mRNA_length"].as<int>() << ".\n";
    //         std::cout << "--L_g_pol=";
    //         print_matrix_double(vm["L_g_pol"].as< vector<double> >()) ;
    //         std::cout <<".\n"  ;
    //         std::cout << "--L_g_tt=";
    //         print_matrix_double(vm["L_g_tt"].as< vector<double> >());
    //         std::cout << ".\n";
    //         std::cout << "\n";
    //         std::cout << "model used: No_abs.\n";

    //         if (vm.count("help")) {
    //             std::cout << desc << "\n";
    //             return 0;}
    //         if (vm.count("multiple_DNA")) {
    //             std::cout << "we are computing multiple DNA strings: "
    //                 << vm["multiple_DNA"].as<bool>() << ".\n"; //print value definition if set
    //             multiple_DNA = vm["multiple_DNA"].as<bool>();
    //             } 
    //         if (vm.count("DIR")) {
    //             std::cout << "output file is set to: "
    //                 << vm["DIR"].as<string>() << ".\n"; //print value definition if set
    //             DIR = vm["DIR"].as<string>();
    //             } 
    //         if (vm.count("number_DNA")) {
    //             std::cout << "number of DNA computed: "
    //                 << vm["number_DNA"].as<int>() << ".\n"; //print value definition if set
    //             number_steps = vm["number_DNA"].as<int>();
    //             }
    //         if (vm.count("x_def")) {
    //             std::cout << "x value: "
    //                 << vm["x_def"].as<string>() << ".\n"; //print value definition if set
    //             x_def = vm["x_def"].as<string>();
    //             }
    //         if (vm.count("M_model")) {
    //             std::cout << "M model: "
    //                 << vm["M_model"].as<string>() << ".\n"; //print value definition if set
    //             M_model = vm["M_model"].as<string>();
    //             }
    //         if (vm.count("mRNA_length")) {
    //             std::cout << "length of the mRNA formed: "
    //                 << vm["mRNA_length"].as<int>() << ".\n"; //print value definition if set
    //             length_mRNA = vm["mRNA_length"].as<int>();
    //             }
    //         if (vm.count("L_g_pol")){
    //             std::cout << "delta G_pol list definition: " ;
    //             print_matrix_double(vm["L_g_pol"].as< vector<double> >());
    //             std::cout << ".\n";
    //             L_delta_g_pol = vm["L_g_pol"].as< vector<double> >();
    //         }
    //         if (vm.count("L_g_tt")){
    //             std::cout << "delta G_tt list definition: " ;
    //             print_matrix_double(vm["L_g_tt"].as< vector<double> >());
    //             std::cout << ".\n";
    //             L_delta_g_tt = vm["L_g_tt"].as< vector<double> >();
    //         }
    //         else {
    //             std::cout << "Not all variables are set \n";}
            
    //     }
    //     catch(exception& e) {
    //         cerr << "error: " << e.what() << "\n";
    //         return 1;}
    //     catch(...) {
    //         cerr << "Exception of unknown type!\n";}
    // }

        // random number generator for uniform distr between 0 and 1
    // std::default_random_engine generator;
    std::random_device seed; 
    // uint32_t seed = device();
    std::mt19937 generator(seed());
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    

    string length_str = std::to_string(length_mRNA);
    fstream fss;
    string parameter_output_file = declare_filename_parameter(multiple_DNA, DIR, M_model, x_def, length_str);
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
    fss << "model used: Big_abs.\n";   
    fss << "the seed:" << seed()<< ". \n";
    fss.close();
    //for making the error landscape defined by different delta G_tt and delta G_pol
    // double L_delta_g_tt[] = {0,1,2,3,4,6,8,10};
    // double L_delta_g_pol[] = {0,1,2,3,4,5,6,7,8,9,10};

    // int length_mRNA = 300;
    // if(x_def == "exp"){
    //     length_mRNA = 300;}
    // else{
    //     length_mRNA = 3000;}
    

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
    int step = 0;
    string x_string = "";
    
    // random number generator for uniform distr between 0 and 1
    // std::default_random_engine generator;
    // std::uniform_real_distribution<double> distribution(0.0,1.0);
    int time_before_loop_begins = time(NULL);
    int number_moves = 0;
    
    for (double delta_g_pol : L_delta_g_pol){
        int rounded_delta_g_pol = round(delta_g_pol);
        string str_rounded_delta_g_pol = std::to_string(rounded_delta_g_pol);
        if(delta_g_pol == 3000){
            delta_g_pol = -log(2);
            str_rounded_delta_g_pol = "-ln(2)";
        }
        for (double delta_g_tt : L_delta_g_tt){
            
            
            
            int rounded_delta_g_tt = round(delta_g_tt);

            // // define rates 
            if(x_def =="exp"){
                x = exp(delta_g_tt + 5);
                x_string = "exp5";}
            else if(x_def == "1"){
                x = 1;
                int rounded_x = round(x);
                x_string = std::to_string(rounded_x);}

            int y = 1;         
            int rounded_y = round(y);
            
            if (M_model == "MF1"){
            a_2r = a_1r = a_1w = 1;
            a_2w = exp(-delta_g_tt);
            b2 = x;
            b1 = x * exp(-delta_g_pol);
            c_2r = c_1r = c_2w = y;
            c_1w = y * exp(-delta_g_tt);
            model_1 = "MF1_x" + x_string;}
            
            if (M_model == "MF2"){
            a_2r = a_1r = a_1w = 1;
            a_2w = exp(-delta_g_tt);
            b2 = x;
            b1 = x * exp(-delta_g_pol);
            c_2r = c_1r = c_1w = y;
            c_2w = y * exp(delta_g_tt);
            model_1 = "MF2_x" + x_string;}

            if (M_model == "MB1"){
            a_2r = a_1r = a_2w = 1;
            a_1w = exp(delta_g_tt);
            b2 = x;
            b1 = x * exp(-delta_g_pol);
            c_2r = c_1r = c_2w = y;
            c_1w = y * exp(-delta_g_tt);
            model_1 = "MB1_x" + x_string;}
            
            if (M_model == "MB2"){
            a_2r = a_1r = a_2w = 1;
            a_1w = exp(delta_g_tt);
            b2 = x;
            b1 = x * exp(-delta_g_pol);
            c_2r = c_1r = c_1w = y;
            c_2w = y * exp(delta_g_tt);
            model_1 = "MB2_x" + x_string;}


            // define probabilities 
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

            fstream fss;
            //for saving multiple steps end error 
            if(multiple_DNA == true){
                string filename = declare_filename(multiple_DNA, DIR, model_1, delta_g_pol, delta_g_tt, step,length_str);
                fss.open(filename.c_str(), ios::out | ios::app);
                // write the file headers
                fss << "mRNA polymer" << "," << "error probability"  << "\n";
            }
          
            for (step = 0; step <= number_steps; step += 1) {
                // for saving all moves to make the polymer 
                if(multiple_DNA == false){
                    string filename = declare_filename(multiple_DNA, DIR, model_1, delta_g_pol, delta_g_tt, step, length_str);
                    fss.open(filename.c_str(), ios::out | ios::app);
                    // write the file headers
                    fss << "monomer added removed" << "," << "length polymer" << "," << "transition state used" <<"," << "error probability"  << "\n";
                }
                std::vector<int> mRNA_string;
                mRNA_string.push_back(1); //we begin with a wrong monomer
                int transition_state = 0;
                structure_mRNA_move mRNA_moves;
                int iteration = 0;
                string transition_state_string;
                int adding_monomer = 0;
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
                            transition_state_string = "t1r ";
                            adding_monomer = 0;
                            number_moves +=1;
                        }
                        else{
                            //go to t_1w
                            transition_state = -1;
                            transition_state_string = "t1w ";
                            adding_monomer = 1;
                            number_moves +=1;
                        }
                    }
                    else{
                        double random_number2 = distribution(generator);
                        // if s[M-1] = r
                        if (mRNA_string[M-1] ==0){
                            if (random_number2 < alpha_r2r){
                                transition_state = 1;
                                transition_state_string = "t1r ";
                                adding_monomer = 0;
                                number_moves +=1;
                            }
                            else if (random_number2 < (alpha_r2r + alpha_r2w)){
                                transition_state = -1;
                                transition_state_string = "t1w ";
                                adding_monomer = 1;
                                number_moves +=1;
                            }
                            else{
                                transition_state = 2;
                                transition_state_string = "- t2r ";
                                adding_monomer = mRNA_string[M];
                                mRNA_string.pop_back();
                                monomer_status = "-";
                                number_moves +=1;
                            }   
                        }
                        // if s[M-1] = w
                        double random_number7 = distribution(generator);
                        if (mRNA_string[M-1] ==1){
                            if (random_number7 < alpha_w2r){
                                transition_state = 1;
                                transition_state_string = "t1r ";
                                adding_monomer = 0;
                                number_moves +=1;
                            }
                            else if (random_number7 < (alpha_w2r + alpha_w2w)){
                                transition_state = -1;
                                transition_state_string = "t1w ";
                                adding_monomer = 1;
                                number_moves +=1;
                            }
                            else{
                                transition_state = -2;
                                transition_state_string = "- t2w ";
                                adding_monomer = mRNA_string[M];
                                mRNA_string.pop_back();
                                monomer_status = "-";
                                number_moves +=1;
                            }   
                        }
                    }

                    // you have to define M again since you could have popped back when you went to transition state 2
                    M = mRNA_string.size() -1;

                    // define from the transition state we are in, which monomer is added or removed.
                    while (transition_state != 0){
                        // if in trabnsition state t1r
                        if (transition_state == 1){
                            double random_number3 = distribution(generator);
                            if (random_number3 < beta_2r){
                                if(mRNA_string[M] == 0){
                                    transition_state = 2;
                                    transition_state_string += "t2r ";
                                    number_moves +=1;
                                }
                                else{ //s[M] = 1
                                    transition_state = -2;
                                    transition_state_string += "t2w ";
                                    number_moves +=1;
                                }
                            }
                            else{
                                // else go back to s[M]
                                transition_state = 0;
                                monomer_status += "n"; 
                                number_moves +=1;
                            }    
                        }

                        // if in transition state 1tw
                        else if (transition_state == -1){
                            double random_number4 = distribution(generator);
                            if(random_number4 < beta_2w){
                                if(mRNA_string[M] == 0){
                                    transition_state = 2;
                                    transition_state_string = transition_state_string + "t2r ";
                                    number_moves +=1;
                                }
                                else{ //s[M] = 1
                                    transition_state = -2;
                                    transition_state_string = transition_state_string + "t2w ";
                                    number_moves +=1;
                                }
                            }
                            else{
                                // else go back to s[M]
                                transition_state = 0;
                                monomer_status += "n"; 
                                number_moves +=1;
                            }
                        }

                        //if in transition state 2tr
                        else if (transition_state == 2){
                            double random_number5 = distribution(generator);
                            if (random_number5 < gamma_2r){
                                if(adding_monomer == 0){
                                    mRNA_string.push_back(0);
                                    transition_state = 0;
                                    monomer_status += "0"; 
                                    number_moves +=1;
                                }
                                else{
                                    mRNA_string.push_back(1);
                                    transition_state = 0;
                                    monomer_status += "1"; 
                                    number_moves +=1;
                                }
                            }
                            else{
                                if(adding_monomer == 0){
                                    transition_state = 1;
                                    transition_state_string = transition_state_string + "t1r ";
                                    number_moves +=1;
                                }
                                else{
                                    transition_state = -1;
                                    transition_state_string = transition_state_string + "t1w ";
                                    number_moves +=1;
                                }   
                            }
                        }


                        //if in transition state 2tw
                        else if (transition_state == -2){
                            double random_number5 = distribution(generator);
                            if (random_number5 < gamma_2w){
                                if(adding_monomer == 0){
                                    mRNA_string.push_back(0);
                                    transition_state = 0;
                                    monomer_status += "0"; 
                                    number_moves +=1;
                                }
                                else{
                                    mRNA_string.push_back(1);
                                    transition_state = 0;
                                    monomer_status += "1"; 
                                    number_moves +=1;
                                }
                            }
                            else{
                                if(adding_monomer == 0){
                                    transition_state = 1;
                                    transition_state_string = transition_state_string + "t1r ";
                                    number_moves +=1;
                                }
                                else{
                                    transition_state = -1;
                                    transition_state_string = transition_state_string + "t1w ";
                                    number_moves +=1;
                                }   
                            }
                        }
                }

                error_prob = count(mRNA_string.begin(), mRNA_string.end(), 1)/(double)mRNA_string.size();
                mRNA_moves.monomer = monomer_status;
                mRNA_moves.length= mRNA_string.size();
                mRNA_moves.error = error_prob;
                mRNA_moves.transition_states = transition_state_string;
            
                //save information for every iteration to CSV file  
                if(multiple_DNA == false){         
                    fss << mRNA_moves.monomer << "," << mRNA_moves.length <<"," << mRNA_moves.transition_states << ","  << mRNA_moves.error;
                    fss << "\n";
                    }
                iteration += 1;
                monomer_status = "";
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
                error_prob = count(mRNA_string.begin() +1, mRNA_string.end()-15, 1)/((double)mRNA_string.size()-16);
                std::vector<int> in = mRNA_string;
                // loop through the array elements
                for(int j=0; j< in.size(); ++j){
                    fss << in.at(j) ;
                    }
                fss << "," << error_prob << "\n";
            }
            
            }
        //save information for all steps, error info
        if(multiple_DNA == true){
            fss << std ::endl;
            fss.close();
        }
        double time_end_loop= time(NULL);
        double time_difference = time_end_loop - time_before_loop_begins;
        double speed = number_moves/time_difference;
        std::cout << "g_pol: " << str_rounded_delta_g_pol << ", g_tt: " << rounded_delta_g_tt << ", speed of simulation (moves/sec): " << speed << ".\n";

        }
    } 
}
