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

int change_base_string_to_int(string base){
    int base_int = 0;
    if(base == "aA"){
        base_int =1;
    }
    else if(base == "T"){
        base_int = 2;
    }
    else if(base == "G"){
        base_int = 3;
    }
    else if(base == "C"){
        base_int = 4;
    }
    else if(base == "iG"){
        base_int = 5;
    }
    else if(base == "iC"){
        base_int = 6;
    }
    else if(base == "S"){
        base_int = 7;
    }
    else if(base == "L"){
        base_int = 8;
    }
    else if(base == "Lambda"){
        base_int = 9;
    }
    else if(base == "Sigma"){
        base_int = 10;
    }
    else if(base == "delta"){
        base_int = 11;
    }
    else if(base == "beta"){
        base_int = 12;
    }
    else if(base == "X"){
        base_int = 13;
    }
    else if(base == "K"){
        base_int = 14;
    }
    return base_int;
}


std::vector<int> base_OH_groups(int base){
    std::vector<int> OH_groups;
    if(base == 1){
        //base = 1 means, base is aA.
        OH_groups = {1, 0, 1};
    }
    else if(base == 2){
        //base = 2 means, base is T.
        OH_groups = {0, 1, 0};
    }
    else if(base == 3){
        //base = 1 means, base is G.
        OH_groups = {0, 1, 1};
    }
    else if(base == 4){
        //base = 1 means, base is C.
        OH_groups = {1, 0, 0};
    }
    else if(base == 5){
        //base = 2 means, base is T.
        OH_groups = {1, 1, 0};
    }
    else if(base == 6){
        //base = 1 means, base is G.
        OH_groups = {0, 0, 1};
    }
    else if(base == 7){
        //base = 1 means, base is C.
        OH_groups = {0, 0, 0};
    }
    else if(base == 8){
        //base = 2 means, base is T.
        OH_groups = {1, 1, 1};
    }
    else if(base == 9){
        //base = 1 means, base is G.
        OH_groups = {1, 1, 1};
    }
    else if(base == 10){
        //base = 1 means, base is C.
        OH_groups = {0, 0, 0};
    }
    else if(base == 11){
        //base = 1 means, base is C.
        OH_groups = {1, 0, 0};
    }
    else if(base == 12){
        //base = 1 means, base is C.
        OH_groups = {0, 1, 1};
    }
    else if(base == 13){
        //base = 1 means, base is C.
        OH_groups = {0, 1, 0};
    }
    else if(base == 14){
        //base = 1 means, base is C.
        OH_groups = {1, 0, 1};
    }
    return OH_groups;
}

double OH_matches_strength(std::vector<int> base_i, std::vector<int> base_j){
    int match = 0;
    for(int i=0; i< base_i.size(); ++i){
        if(base_i[i] != base_j[i]){
            match += 1;
        }
    }
    return -1*(2 * match - 3);
}

double a_plus(int base_template, int base_copy,double repulsion_same_base){
    double strength;
    string templated;
    string copy;
    if (base_template % 2 == 0){
        templated ="even";
    }
    else{
        templated = "odd";
    }
    if (base_copy % 2 == 0){
        copy = "even";
    }
    else{
        copy = "odd";
    }
    if(repulsion_same_base>=0){
        if(copy =="even" && templated=="even"){
            strength = 1;
        }
        else if(copy =="odd" && templated=="odd"){
            strength = 1;
        }
        else{
            strength = 1;
        }
    }
    else{
        if(copy =="even" && templated=="even"){
            strength = 0;
        }
        else if(copy =="odd" && templated=="odd"){
            strength = 0;
        }
        else{
            strength = 1;
        }
    }
    
    return strength;
}

double a_min(int base_template, int base_copy, double repulsion_same_base){
    double strength;
    string templated;
    string copy;
    if (base_template % 2 == 0){
        templated ="even";
    }
    else{
        templated = "odd";
    }
    if (base_copy % 2 == 0){
        copy = "even";
    }
    else{
        copy = "odd";
    }
    std::vector<int> base_i = base_OH_groups(base_template);
    std::vector<int> base_j = base_OH_groups(base_copy);
    double OH_strength = OH_matches_strength(base_i, base_j);
    if(repulsion_same_base >=0){
        if(copy =="even" && templated=="even"){
            OH_strength += repulsion_same_base;
            strength = exp(OH_strength);
        }
        else if(copy =="odd" && templated=="odd"){
            OH_strength += repulsion_same_base;
            strength = exp(OH_strength);
        }
        else{
            strength = exp(OH_strength);
        }
    }
    else{
        strength = exp(OH_strength);
        // if(copy =="even" && templated=="even"){
        //     strength = 0;
        // }
        // else if(copy =="odd" && templated=="odd"){
        //     strength = 0;
        // }
        // else{
        //     strength = exp(OH_strength);
        // }
    }
    return strength;
}

double b_plus(){
    return 1;
}

double b_min(double G_bb){
    return exp(-G_bb);
}

double c_plus(int base_template, int base_copy, double repulsion_same_base){
    double strength;
    string templated;
    string copy;
    if (base_template % 2 == 0){
        templated ="even";
    }
    else{
        templated = "odd";
    }
    if (base_copy % 2 == 0){
        copy = "even";
    }
    else{
        copy = "odd";
    }
    std::vector<int> base_i = base_OH_groups(base_template);
    std::vector<int> base_j = base_OH_groups(base_copy);
    double OH_strength = OH_matches_strength(base_i, base_j);
    if (repulsion_same_base >=0){
        if(copy =="even" && templated=="even"){
            OH_strength += repulsion_same_base;
            strength = exp(OH_strength);
        }
        else if(copy =="odd" && templated=="odd"){
            OH_strength += repulsion_same_base;
            strength = exp(OH_strength);
        }
        else{
            strength = exp(OH_strength);
        }
    }
    else{
        strength = exp(OH_strength);
        // if(copy =="even" && templated=="even"){
        // strength = 0;
        // }
        // else if(copy =="odd" && templated=="odd"){
        //     strength = 0;
        // }
        // else{
        //     strength = exp(OH_strength);
        // }
    }
    return strength;
}

double c_min(){
    return 1;
}


string declare_filename(bool all_steps, string DIR, string x_def, double delta_g_pol, int step, string purine_1, string pyrimidines_1, string purine_2, string pyrimidines_2, vector<double> ratio, double repulsion_same_base){
    string date_now = date(time(0));
    int rounded_delta_g_pol = round(delta_g_pol);
    string str_rounded_delta_g_pol = std::to_string(rounded_delta_g_pol);
    if(delta_g_pol == -log(2)){
        str_rounded_delta_g_pol = "-ln2";
    }
    if(delta_g_pol == -log(4)){
        str_rounded_delta_g_pol = "-ln4";
    }
    int rounded_ratio_0 = round(ratio[0]);
    string rounded_ratio_0_str = std::to_string(rounded_ratio_0);
    int rounded_ratio_1 = round(ratio[1]);
    string rounded_ratio_1_str = std::to_string(rounded_ratio_1);
    int rounded_ratio_2 = round(ratio[2]);
    string rounded_ratio_2_str = std::to_string(rounded_ratio_2);
    int rounded_ratio_3 = round(ratio[3]);
    string rounded_ratio_3_str = std::to_string(rounded_ratio_3);
    int rounded_repulsion = round(repulsion_same_base);
    string rounded_repulsion_str = std::to_string(rounded_repulsion);

    string filename_output = "";
    if(all_steps ==true){
        string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/"+ DIR + "/";
        filename_output = path + date_now + "model_muriel_big_abs_BF_x" + x_def + purine_1 + pyrimidines_1 + purine_2 + pyrimidines_2 + rounded_ratio_0_str + rounded_ratio_1_str + rounded_ratio_2_str + rounded_ratio_3_str + "_delta_g_pol_" + str_rounded_delta_g_pol + "_repulsion" + rounded_repulsion_str +"_allsteps.csv";
    }
    else if(all_steps == false){
        int step_rounded = round(step);
        string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/" + DIR + "/";
        filename_output = path + date_now + "model_muriel_big_abs_step" + std::to_string(step_rounded) + "_x" + x_def + purine_1 + pyrimidines_1 + purine_2 + pyrimidines_2 + rounded_ratio_0_str + rounded_ratio_1_str + rounded_ratio_2_str + rounded_ratio_3_str  + "_delta_g_pol_" + str_rounded_delta_g_pol + "_repulsion" + rounded_repulsion_str + ".csv";            
    }
    return filename_output;
}

string declare_filename_parameter(bool all_steps, string DIR, string x_def, string purine_1, string pyrimidines_1, string purine_2, string pyrimidines_2, vector<double> ratio, double repulsion_same_base){
    string date_now = date(time(0));
    int rounded_repulsion = round(repulsion_same_base);
    string rounded_repulsion_str = std::to_string(rounded_repulsion);
    string filename_output = "";
    string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/" + DIR + "/";
    if(all_steps ==true){
        // filename_output = path + date_now + "parameter_file_model_big_abs_BF_x" + x_def + purine_1 + pyrimidines_1 + purine_2 + pyrimidines_2 + ratio[0] + ratio[1] + ratio[2] + ratio[3] +"_allsteps.csv";
        filename_output = path + date_now + "parameter_file_model_big_abs_BF_x" + x_def + purine_1 + pyrimidines_1 + purine_2 + pyrimidines_2 + "_repulsion" + rounded_repulsion_str +"_allsteps.csv";}
    else{
        // filename_output = path + date_now + "parameter_file_model_big_abs_BF_x" + x_def + purine_1 + pyrimidines_1 + purine_2 + pyrimidines_2 + ratio[0] + ratio[1] + ratio[2] + ratio[3] + "_seperate_steps.csv";
        filename_output = path + date_now + "parameter_file_model_big_abs_BF_x" + x_def + purine_1 + pyrimidines_1 + purine_2 + pyrimidines_2  + "_repulsion" + rounded_repulsion_str + "_seperate_steps.csv";}
     
    return filename_output;
}

std::vector<int> template_string(int string_length, vector<double> ratio_bases_template, int purine_1, int purine_2, int pyrimidines_1, int pyrimidines_2){
    // create a (uniformly) randomized template string 
    std::vector<int> string;
    std::default_random_engine generator;
    double ratios = ratio_bases_template[0] +ratio_bases_template[1]+ratio_bases_template[2]+ratio_bases_template[3];
    std::uniform_real_distribution<double> distribution(0,ratios);
    for(int j=0; j<string_length; j++){
        double random_number = distribution(generator);
        if(random_number < ratio_bases_template[0]){
            string.push_back(purine_1);
        }
        else if(random_number < ratio_bases_template[0]+ratio_bases_template[1]){
            string.push_back(pyrimidines_1);
        }
        else if(random_number < ratio_bases_template[0]+ratio_bases_template[1]+ratio_bases_template[2]){
            string.push_back(purine_2);
        }
        else{
            string.push_back(pyrimidines_2);
        }
    }
    return string;
}

string declare_polymers_string(bool all_steps, string DIR, string x_def, string purine_1, string pyrimidines_1, string purine_2, string pyrimidines_2, vector<double> ratio, double repulsion_same_base, double delta_g_pol){
    string date_now = date(time(0));

    int rounded_ratio_0 = round(ratio[0]);
    string rounded_ratio_0_str = std::to_string(rounded_ratio_0);
    int rounded_ratio_1 = round(ratio[1]);
    string rounded_ratio_1_str = std::to_string(rounded_ratio_1);
    int rounded_ratio_2 = round(ratio[2]);
    string rounded_ratio_2_str = std::to_string(rounded_ratio_2);
    int rounded_ratio_3 = round(ratio[3]);
    string rounded_ratio_3_str = std::to_string(rounded_ratio_3);
    int rounded_repulsion = round(repulsion_same_base);
    string rounded_repulsion_str = std::to_string(rounded_repulsion);
    int rounded_delta_g_pol = round(delta_g_pol);
    string str_rounded_delta_g_pol = std::to_string(rounded_delta_g_pol);
    if(delta_g_pol == -log(2)){
        str_rounded_delta_g_pol = "-ln2";
    }
    if(delta_g_pol == -log(4)){
        str_rounded_delta_g_pol = "-ln4";
    }
    string filename_output = "";
    string path = "/home/ipausers/louman/Documents/programming/DNA_replication_muriel/outs/" + DIR + "/";
    filename_output = path + date_now + "polymer_strings_model_muriel_big_abs_BF_x" + x_def + purine_1 + pyrimidines_1 + purine_2 + pyrimidines_2 + rounded_ratio_0_str + rounded_ratio_1_str + rounded_ratio_2_str + rounded_ratio_3_str  + "_delta_g_pol_" + str_rounded_delta_g_pol + "_repulsion" + rounded_repulsion_str +"_allsteps.csv";
    return filename_output;
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
        else if(i == 5 && j==6){
            E = 5;
            }
        else if(i == 6 && j==5){
            E = 5;}
        else if(i == 7 && j==8){
            E = 5;}
        else if(i == 8 && j==7){
            E = 5;
            }
        else if(i == 9 && j==10){
            E = 5;}
        else if(i == 10 && j==9){
            E = 5;}
        else if(i == 11 && j==12){
            E = 5;}
        else if(i == 12 && j==11){
            E = 5;}
        else if(i == 13 && j==14){
            E = 5;}
        else if(i == 14 && j==13){
            E = 5;}
        else{
            E = -5;}
    return E;
    }

std::vector<int> base_probs(std::vector<int> template_str, std::vector<int> copy, int purine_1_int, int purine_2_int, int pyrimidine_1_int, int pyrimidine_2_int){
    int length_string = template_str.size() -1;
    int pur1pur1 = 0;
    int pur1pur2 = 0;
    int pur2pur2 = 0;
    int py1py1 = 0;
    int py1py2 = 0;
    int py2py2 = 0;
    int py1pur1 = 0;
    int py1pur2 = 0;
    int py2pur1 = 0;
    int py2pur2 = 0;
    for (int i = 2; i < length_string; i++){
        int base_template = template_str[i];
        int base_copy = copy[i];
        if(base_template == purine_1_int && base_copy ==purine_1_int){
            pur1pur1 += 1;
        }
        else if(base_template == purine_1_int && base_copy ==purine_2_int){
            pur1pur2 += 1;
        }
        else if(base_template == purine_2_int && base_copy ==purine_1_int){
            pur1pur2 += 1;
        }
        else if(base_template == purine_2_int && base_copy ==purine_2_int){
            pur2pur2 += 1;
        }
        else if(base_template == pyrimidine_1_int && base_copy ==pyrimidine_1_int){
            py1py1 += 1;
        }
        else if(base_template == pyrimidine_1_int && base_copy ==pyrimidine_2_int){
            py1py2 += 1;
        }
        else if(base_template == pyrimidine_2_int && base_copy ==pyrimidine_1_int){
            py1py2 += 1;
        }
        else if(base_template == pyrimidine_2_int && base_copy ==pyrimidine_2_int){
            py2py2 += 1;
        }
        else if(base_template == purine_1_int && base_copy ==pyrimidine_1_int){
            py1pur1 += 1;
        }
        else if(base_template == pyrimidine_1_int && base_copy ==purine_1_int){
            py1pur1 += 1;
        }
        else if(base_template == purine_2_int && base_copy ==pyrimidine_1_int){
            py1pur2 += 1;
        }
        else if(base_template == pyrimidine_1_int && base_copy ==purine_2_int){
            py1pur2 += 1;
        }
        else if(base_template == purine_1_int && base_copy ==pyrimidine_2_int){
            py2pur1 += 1;
        }
        else if(base_template == pyrimidine_2_int && base_copy ==purine_1_int){
            py2pur1 += 1;
        }
        else if(base_template == purine_2_int && base_copy ==pyrimidine_2_int){
            py2pur2 += 1;
        }
        else if(base_template == pyrimidine_2_int && base_copy ==purine_2_int){
            py2pur2 += 1;
        }
    }
    return {pur1pur1, pur1pur2, pur2pur2, py1py1, py1py2, py2py2, py1pur1, py1pur2, py2pur1, py2pur2};
}

int main(int ac, char* av[]){
    //in main
    // string DIR = "230417output"; //"230331"
    // bool multiple_DNA = true;
    // int number_steps = 100;
    // string x_def = "1";
    // int length_mRNA = 300;
    // vector<double> L_delta_g_pol = {0,1,2,3,4,5,6,7,8,9,10,20};
    // vector<double> ratio_bases_template = {1,1,1,1};
    // string purine_1 = "iG";
    // string purine_2 = "X";
    // string pyrimidines_1 = "iC";
    // string pyrimidines_2 =  "K";
    // double repulsion_00_11 = 1;


    bool variables_in_main =true;
    po::variables_map vm;
    string DIR = "";
    bool multiple_DNA;
    int number_steps;
    string x_def = "";
    int length_mRNA;
    vector<double> L_delta_g_pol;
    vector<double> ratio_bases_template;
    string purine_1 = "";
    string purine_2 = "";
    string pyrimidines_1 = "";
    string pyrimidines_2 = "";
    int repulsion_00_11;

    if(variables_in_main == true){
        po::options_description desc("Allowed options");
        try {

            po::options_description desc("Allowed options");
            desc.add_options()
                ("help", "produce help message")
                ("multiple_DNA", po::value<bool>(), "compute multiple DNA string") //set declaration of values
                ("DIR", po::value<string>(), "output directory name") 
                ("number_DNA", po::value<int>(), "number of DNA computed")
                ("x_def", po::value<string>(), "strength x, backbone") 
                ("mRNA_length", po::value<int>(), "length of the mRNA string formed")
                ("L_g_pol", po::value< vector<double> >(), "delta G_pol definitions list") 
                ("ratio_template", po::value< vector<double> >(), "the ratio of the four bases we are using")
                ("purine_1", po::value<string>(), "purine 1 base definition")
                ("purine_2", po::value<string>(), "purine 2 base definition")
                ("pyrimidines_1", po::value<string>(), "pyrimidines 1 base definition")
                ("pyrimidines_2", po::value<string>(), "pyrimidines 2 base definition")
                ("repulsion_11_00", po::value<int>(), "strength repulsion purine-pruine or pyrimidine-pyrimidine binding")


            ;
            po::store(po::parse_command_line(ac, av, desc), vm);
            po::notify(vm);

            std::cout << "--multiple_DNA=" << vm["multiple_DNA"].as<bool>() << ".\n" ;
            std::cout << "--DIR=" << vm["DIR"].as<string>() << ".\n";
            std::cout << "--number_DNA=" << vm["number_DNA"].as<int>() << ".\n";
            std::cout << "--x_def=" << vm["x_def"].as<string>() << ".\n";
            std::cout << "--mRNA_length=" << vm["mRNA_length"].as<int>() << ".\n";
            std::cout << "--L_g_pol=";
            print_matrix_double(vm["L_g_pol"].as< vector<double> >()) ;
            std::cout <<".\n"  ;
            std::cout << "--ratio_template=";
            print_matrix_double(vm["ratio_template"].as< vector<double> >());
            std::cout << ".\n";
            std::cout << "--purine_1=" << vm["purine_1"].as<string>() << ".\n";
            std::cout << "--purine_2=" << vm["purine_2"].as<string>() << ".\n";
            std::cout << "--pyrimidines_1=" << vm["pyrimidines_1"].as<string>() << ".\n";
            std::cout << "--pyrimidines_2=" << vm["pyrimidines_2"].as<string>() << ".\n";
            std::cout << "--repulsion_11_00=" << vm["repulsion_11_00"].as<int>() << ".\n";
            std::cout << "\n";
            std::cout << "model used: Big_abs.\n";

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
            if (vm.count("x_def")) {
                std::cout << "x value: "
                    << vm["x_def"].as<string>() << ".\n"; //print value definition if set
                x_def = vm["x_def"].as<string>();
                }
            if (vm.count("mRNA_length")) {
                std::cout << "length of the mRNA formed: "
                    << vm["mRNA_length"].as<int>() << ".\n"; //print value definition if set
                length_mRNA = vm["mRNA_length"].as<int>();
                }
            if (vm.count("L_g_pol")){
                std::cout << "delta G_pol list definition: " ;
                print_matrix_double(vm["L_g_pol"].as< vector<double> >());
                std::cout << ".\n";
                L_delta_g_pol = vm["L_g_pol"].as< vector<double> >();
            }
            if (vm.count("ratio_template")){
                std::cout << "the ratio of the four bases we are using: " ;
                print_matrix_double(vm["ratio_template"].as< vector<double> >());
                std::cout << ".\n";
                ratio_bases_template = vm["ratio_template"].as< vector<double> >();
            }
            if (vm.count("purine_1")) {
                std::cout << "purine 1 base definition: "
                    << vm["purine_1"].as<string>() << ".\n"; //print value definition if set
                purine_1 = vm["purine_1"].as<string>();
                } 
            if (vm.count("purine_2")) {
                std::cout << "purine 2 base definition: "
                    << vm["purine_2"].as<string>() << ".\n"; //print value definition if set
                purine_2 = vm["purine_2"].as<string>();
                } 
            if (vm.count("pyrimidines_1")) {
                std::cout << "pyrimidines 1 base definition: "
                    << vm["pyrimidines_1"].as<string>() << ".\n"; //print value definition if set
                pyrimidines_1 = vm["pyrimidines_1"].as<string>();
                } 
            if (vm.count("pyrimidines_2")) {
                std::cout << "pyrimidines 2 base definition: "
                    << vm["pyrimidines_2"].as<string>() << ".\n"; //print value definition if set
                pyrimidines_2 = vm["pyrimidines_2"].as<string>();
                } 
            if (vm.count("repulsion_11_00")) {
                std::cout << "strength repulsion purine-purine or pyrimidine-pyrimi9dine binding: "
                    << vm["repulsion_11_00"].as<int>() << ".\n"; //print value definition if set
                repulsion_00_11 = vm["repulsion_11_00"].as<int>();
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

     // define variables
    string model_1 = "";
    double x = 0;
    int step = 0;
    string x_string = "";
    int purine_1_int = change_base_string_to_int(purine_1);
    int pyrimidines_1_int = change_base_string_to_int(pyrimidines_1);
    int purine_2_int = change_base_string_to_int(purine_2);
    int pyrimidines_2_int = change_base_string_to_int(pyrimidines_2);

    // random number generator for uniform distr between 0 and 1
    // std::default_random_engine generator;
    std::random_device seed; 
    // uint32_t seed = device();
    std::mt19937 generator(seed());
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    std::vector<int> mRNA_template = template_string(length_mRNA, ratio_bases_template, purine_1_int, purine_2_int, pyrimidines_1_int, pyrimidines_2_int);
      

    fstream fss;
    string parameter_output_file = declare_filename_parameter(multiple_DNA, DIR, x_def, purine_1, pyrimidines_1, purine_2, pyrimidines_2, ratio_bases_template, repulsion_00_11);
    fss.open(parameter_output_file.c_str(), ios::out | ios::app);
    fss << "--multiple_DNA=" << multiple_DNA << ".\n" ;
    fss << "--DIR=" << DIR << ".\n";
    fss << "--number_DNA=" << number_steps << ".\n";
    fss << "--x_def=" << x_def << ".\n";
    fss << "--mRNA_length=" << length_mRNA << ".\n";
    fss << "--L_g_pol=";
    std::vector<double> in = L_delta_g_pol;
    // loop through the array elements
    for(int j=0; j< in.size(); ++j){
        fss << in.at(j);
        fss << " ";}
    fss <<".\n";
    fss << "--ratio_template=";
    std::vector<double> in3 = ratio_bases_template;
    // loop through the array elements
    for(int j=0; j< in3.size(); ++j){
        fss << in3.at(j);
        fss << " ";}       
    fss << ".\n";
    fss << "--purine_1=" << purine_1 << ".\n";
    fss << "--purine_2=" << purine_2 << ".\n";
    fss << "--pyrimidines_1=" << pyrimidines_1 << ".\n";
    fss << "--pyrimidines_2=" << pyrimidines_2 << ".\n";
    fss << "--repulsion_00_11=" << repulsion_00_11 << ".\n";

    fss << "\n";
    fss << "model used: Big_abs.\n";   
    fss << "template string:" ;
    std::vector<int> in4 = mRNA_template;
    // loop through the array elements
    for(int j=0; j< in4.size(); ++j){
        fss << in4.at(j);
        fss << " ";}       
    fss << ".\n";
    fss << "the seed:" << seed()<< ". \n";
    fss.close();

   
    
    //loop over different g_tt and g_bb
    for (double delta_g_pol : L_delta_g_pol){
        if(delta_g_pol == 3000){
            // delta_g_pol = -log(2);
            if(repulsion_00_11==-1){
                delta_g_pol = -log(2);}
            else if(repulsion_00_11 > -1){
                delta_g_pol = -log(4);
            }
        }
        
        if(x_def == "1"){
            x = 1;
            int rounded_x = round(x);
            x_string = std::to_string(rounded_x);}          
        
        std::vector<int> mRNA_string_out;
        double A1m, A2p,A2m, A3p,A3m, A4p,A4m, A5p,A5m, B1p,B1m, B2p,B2m, B3p,B3m, B4p,B4m, B5p,B5m, G1p, G1m, G2p, G3p, G4p, G5p;
        A1m=A2p=A2m=A3p=A3m=A4p=A4m=A5p=A5m=B1p=B1m=B2p=B2m=B3p=B3m=B4p=B4m=B5p=B5m=G1p=G1m=G2p=G3p=G4p=G5p = 0;
        
        fstream fss;
        
        
        // for saving multiple steps end error 
        if(multiple_DNA == true){
            string filename = declare_filename(multiple_DNA, DIR, x_def, delta_g_pol, step, purine_1, pyrimidines_1, purine_2, pyrimidines_2, ratio_bases_template, repulsion_00_11);
            fss.open(filename.c_str(), ios::out | ios::app);
            // write the file headers
            fss << "mRNA polymer" << "," << "error probability"  << "\n";
                }

        for (step = 0; step <= number_steps; step += 1) {
            // for saving all moves to make the polymer 
            if(multiple_DNA == false){
                string filename = declare_filename(multiple_DNA, DIR, x_def, delta_g_pol, step, purine_1, pyrimidines_1, purine_2, pyrimidines_2, ratio_bases_template, repulsion_00_11);
                fss.open(filename.c_str(), ios::out | ios::app);
                // write the file headers
                fss << "monomer added removed" << "," << "length polymer" << "," << "transition state used" <<"," << "error probability"  << "\n";
            }
            
            std::vector<int> mRNA_string;
            double error_prob = 0;
            std::vector<int> error_string;
            mRNA_string.push_back(purine_1_int); 
            error_string.push_back(error(mRNA_string[0], mRNA_template[0]));
            mRNA_string.push_back(purine_1_int); //we begin with two right monomers
            error_string.push_back(error(mRNA_string[1], mRNA_template[1]));
            
            
            int transition_state = 0;
            structure_mRNA_move mRNA_moves;
            int iteration = 0;
            string transition_state_string = "t";
            string monomer_status = "0";   
            int monomer_adding = 0;


            while (mRNA_string.size()<length_mRNA){
                // Go to transition state 1r or 1w or 2, first MC step
                // create a random number betwen 0 and 1
                
                int Q = mRNA_string.size() -1;
                

                // define the variables used in the absorption matrix

                
                
                A1m= a_min(mRNA_template[Q], mRNA_string[Q], repulsion_00_11)/(a_min(mRNA_template[Q], mRNA_string[Q],repulsion_00_11) + b_plus());
                B1p= b_plus() /(b_plus() + a_min(mRNA_template[Q], mRNA_string[Q], repulsion_00_11));
                B1m= b_min(delta_g_pol)/(b_min(delta_g_pol) + c_plus(mRNA_template[Q-1], mRNA_string[Q-1],repulsion_00_11));
                G1p= c_plus(mRNA_template[Q-1], mRNA_string[Q-1],repulsion_00_11)/(c_plus(mRNA_template[Q-1], mRNA_string[Q-1],repulsion_00_11) + b_min(delta_g_pol));
                G1m= c_min() /(c_min()+ a_plus(mRNA_template[Q+1], purine_1_int, repulsion_00_11)+ a_plus(mRNA_template[Q+1], pyrimidines_1_int, repulsion_00_11)+ a_plus(mRNA_template[Q+1], purine_2_int, repulsion_00_11)+a_plus(mRNA_template[Q+1], pyrimidines_2_int, repulsion_00_11));
                
                A2p= a_plus(mRNA_template[Q+1], purine_1_int,repulsion_00_11) /(c_min()+ a_plus(mRNA_template[Q+1], purine_1_int, repulsion_00_11)+ a_plus(mRNA_template[Q+1], pyrimidines_1_int, repulsion_00_11)+ a_plus(mRNA_template[Q+1], purine_2_int, repulsion_00_11)+a_plus(mRNA_template[Q+1], pyrimidines_2_int, repulsion_00_11));
                A2m= a_min(mRNA_template[Q+1], purine_1_int,repulsion_00_11)/(a_min(mRNA_template[Q+1], purine_1_int,repulsion_00_11)+ b_plus());
                B2p= b_plus() /(b_plus() + a_min(mRNA_template[Q+1], purine_1_int,repulsion_00_11));
                B2m= b_min(delta_g_pol)/(b_min(delta_g_pol) + c_plus(mRNA_template[Q], mRNA_string[Q],repulsion_00_11));
                G2p= c_plus(mRNA_template[Q], mRNA_string[Q],repulsion_00_11)/(c_plus(mRNA_template[Q], mRNA_string[Q],repulsion_00_11) + b_min(delta_g_pol));
                
                A3p= a_plus(mRNA_template[Q+1], pyrimidines_1_int,repulsion_00_11) /(c_min()+ a_plus(mRNA_template[Q+1], purine_1_int, repulsion_00_11)+ a_plus(mRNA_template[Q+1], pyrimidines_1_int, repulsion_00_11)+ a_plus(mRNA_template[Q+1], purine_2_int, repulsion_00_11)+a_plus(mRNA_template[Q+1], pyrimidines_2_int, repulsion_00_11));
                A3m= a_min(mRNA_template[Q+1], pyrimidines_1_int,repulsion_00_11)/(a_min(mRNA_template[Q+1], pyrimidines_1_int,repulsion_00_11)+ b_plus());
                B3p= b_plus() /(b_plus() + a_min(mRNA_template[Q+1], pyrimidines_1_int,repulsion_00_11));
                B3m= b_min(delta_g_pol)/(b_min(delta_g_pol) + c_plus(mRNA_template[Q], mRNA_string[Q],repulsion_00_11));
                G3p= c_plus(mRNA_template[Q], mRNA_string[Q],repulsion_00_11)/(c_plus(mRNA_template[Q], mRNA_string[Q],repulsion_00_11) + b_min(delta_g_pol));
                
                A4p= a_plus(mRNA_template[Q+1], purine_2_int,repulsion_00_11) /(c_min()+ a_plus(mRNA_template[Q+1], purine_1_int, repulsion_00_11)+ a_plus(mRNA_template[Q+1], pyrimidines_1_int, repulsion_00_11)+ a_plus(mRNA_template[Q+1], purine_2_int, repulsion_00_11)+a_plus(mRNA_template[Q+1], pyrimidines_2_int, repulsion_00_11));
                A4m= a_min(mRNA_template[Q+1], purine_2_int,repulsion_00_11)/(a_min(mRNA_template[Q+1], purine_2_int,repulsion_00_11)+ b_plus());
                B4p= b_plus() /(b_plus() + a_min(mRNA_template[Q+1], purine_2_int,repulsion_00_11));
                B4m= b_min(delta_g_pol)/(b_min(delta_g_pol) + c_plus(mRNA_template[Q], mRNA_string[Q],repulsion_00_11));
                G4p= c_plus(mRNA_template[Q], mRNA_string[Q],repulsion_00_11)/(c_plus(mRNA_template[Q], mRNA_string[Q],repulsion_00_11) + b_min(delta_g_pol));

                A5p= a_plus(mRNA_template[Q+1], pyrimidines_2_int,repulsion_00_11) /(c_min()+ a_plus(mRNA_template[Q+1], purine_1_int, repulsion_00_11)+ a_plus(mRNA_template[Q+1], pyrimidines_1_int, repulsion_00_11)+ a_plus(mRNA_template[Q+1], purine_2_int, repulsion_00_11)+a_plus(mRNA_template[Q+1], pyrimidines_2_int, repulsion_00_11));
                A5m= a_min(mRNA_template[Q+1], pyrimidines_2_int,repulsion_00_11)/(a_min(mRNA_template[Q+1], pyrimidines_2_int,repulsion_00_11) + b_plus());
                B5p= b_plus() /(b_plus() + a_min(mRNA_template[Q+1], pyrimidines_2_int,repulsion_00_11));
                B5m= b_min(delta_g_pol)/(b_min(delta_g_pol) + c_plus(mRNA_template[Q], mRNA_string[Q],repulsion_00_11));
                G5p = c_plus(mRNA_template[Q], mRNA_string[Q],repulsion_00_11)/(c_plus(mRNA_template[Q], mRNA_string[Q],repulsion_00_11) + b_min(delta_g_pol));
                

                
                __float128 Z = (1 - A2m*A2p - A3m*A3p - A4m*A4p - A5m*A5p - B1m*B1p + A2m*A2p*B1m*B1p + A3m*A3p*B1m*B1p + A4m*A4p*B1m*B1p + A5m*A5p*B1m*B1p 
                - B2m*B2p + A3m*A3p*B2m*B2p + A4m*A4p*B2m*B2p + A5m*A5p*B2m*B2p + B1m*B1p*B2m*B2p - A3m*A3p*B1m*B1p*B2m*B2p - A4m*A4p*B1m*B1p*B2m*B2p 
                - A5m*A5p*B1m*B1p*B2m*B2p - B3m*B3p + A2m*A2p*B3m*B3p + A4m*A4p*B3m*B3p + A5m*A5p*B3m*B3p + B1m*B1p*B3m*B3p - A2m*A2p*B1m*B1p*B3m*B3p 
                - A4m*A4p*B1m*B1p*B3m*B3p -  A5m*A5p*B1m*B1p*B3m*B3p + B2m*B2p*B3m*B3p - A4m*A4p*B2m*B2p*B3m*B3p - A5m*A5p*B2m*B2p*B3m*B3p 
                - B1m*B1p*B2m*B2p*B3m*B3p + A4m*A4p*B1m*B1p*B2m*B2p*B3m*B3p + A5m*A5p*B1m*B1p*B2m*B2p*B3m*B3p - B4m*B4p + A2m*A2p*B4m*B4p + A3m*A3p*B4m*B4p 
                + A5m*A5p*B4m*B4p + B1m*B1p*B4m*B4p - A2m*A2p*B1m*B1p*B4m*B4p - A3m*A3p*B1m*B1p*B4m*B4p - A5m*A5p*B1m*B1p*B4m*B4p + B2m*B2p*B4m*B4p 
                - A3m*A3p*B2m*B2p*B4m*B4p - A5m*A5p*B2m*B2p*B4m*B4p - B1m*B1p*B2m*B2p*B4m*B4p + A3m*A3p*B1m*B1p*B2m*B2p*B4m*B4p 
                + A5m*A5p*B1m*B1p*B2m*B2p*B4m*B4p + B3m*B3p*B4m*B4p - A2m*A2p*B3m*B3p*B4m*B4p - A5m*A5p*B3m*B3p*B4m*B4p - B1m*B1p*B3m*B3p*B4m*B4p 
                + A2m*A2p*B1m*B1p*B3m*B3p*B4m*B4p + A5m*A5p*B1m*B1p*B3m*B3p*B4m*B4p - B2m*B2p*B3m*B3p*B4m*B4p + A5m*A5p*B2m*B2p*B3m*B3p*B4m*B4p 
                + B1m*B1p*B2m*B2p*B3m*B3p*B4m*B4p - A5m*A5p*B1m*B1p*B2m*B2p*B3m*B3p*B4m*B4p - B5m*B5p + A2m*A2p*B5m*B5p + A3m*A3p*B5m*B5p 
                + A4m*A4p*B5m*B5p + B1m*B1p*B5m*B5p - A2m*A2p*B1m*B1p*B5m*B5p - A3m*A3p*B1m*B1p*B5m*B5p - A4m*A4p*B1m*B1p*B5m*B5p + B2m*B2p*B5m*B5p 
                - A3m*A3p*B2m*B2p*B5m*B5p - A4m*A4p*B2m*B2p*B5m*B5p - B1m*B1p*B2m*B2p*B5m*B5p + A3m*A3p*B1m*B1p*B2m*B2p*B5m*B5p 
                + A4m*A4p*B1m*B1p*B2m*B2p*B5m*B5p + B3m*B3p*B5m*B5p - A2m*A2p*B3m*B3p*B5m*B5p - A4m*A4p*B3m*B3p*B5m*B5p -  B1m*B1p*B3m*B3p*B5m*B5p 
                + A2m*A2p*B1m*B1p*B3m*B3p*B5m*B5p + A4m*A4p*B1m*B1p*B3m*B3p*B5m*B5p - B2m*B2p*B3m*B3p*B5m*B5p + A4m*A4p*B2m*B2p*B3m*B3p*B5m*B5p 
                + B1m*B1p*B2m*B2p*B3m*B3p*B5m*B5p - A4m*A4p*B1m*B1p*B2m*B2p*B3m*B3p*B5m*B5p + B4m*B4p*B5m*B5p - A2m*A2p*B4m*B4p*B5m*B5p 
                - A3m*A3p*B4m*B4p*B5m*B5p - B1m*B1p*B4m*B4p*B5m*B5p + A2m*A2p*B1m*B1p*B4m*B4p*B5m*B5p + A3m*A3p*B1m*B1p*B4m*B4p*B5m*B5p 
                - B2m*B2p*B4m*B4p*B5m*B5p + A3m*A3p*B2m*B2p*B4m*B4p*B5m*B5p + B1m*B1p*B2m*B2p*B4m*B4p*B5m*B5p - A3m*A3p*B1m*B1p*B2m*B2p*B4m*B4p*B5m*B5p 
                - B3m*B3p*B4m*B4p*B5m*B5p + A2m*A2p*B3m*B3p*B4m*B4p*B5m*B5p + B1m*B1p*B3m*B3p*B4m*B4p*B5m*B5p - A2m*A2p*B1m*B1p*B3m*B3p*B4m*B4p*B5m*B5p 
                + B2m*B2p*B3m*B3p*B4m*B4p*B5m*B5p - B1m*B1p*B2m*B2p*B3m*B3p*B4m*B4p*B5m*B5p - G1m*G1p+B2m*B2p*G1m*G1p + B3m*B3p*G1m*G1p 
                - B2m*B2p*B3m*B3p*G1m*G1p + B4m*B4p*G1m*G1p - B2m*B2p*B4m*B4p*G1m*G1p - B3m*B3p*B4m*B4p*G1m*G1p + B2m*B2p*B3m*B3p*B4m*B4p*G1m*G1p 
                + B5m*B5p*G1m*G1p - B2m*B2p*B5m*B5p*G1m*G1p - B3m*B3p*B5m*B5p*G1m*G1p + B2m*B2p*B3m*B3p*B5m*B5p*G1m*G1p - B4m*B4p*B5m*B5p*G1m*G1p 
                + B2m*B2p*B4m*B4p*B5m*B5p*G1m*G1p + B3m*B3p*B4m*B4p*B5m*B5p*G1m*G1p - B2m*B2p*B3m*B3p*B4m*B4p*B5m*B5p*G1m*G1p);
                __float128 backward = (A1m*(B1m*G1m - B1m*B2m*B2p*G1m - B1m*B3m*B3p*G1m + B1m*B2m*B2p*B3m*B3p*G1m - B1m*B4m*B4p*G1m 
                + B1m*B2m*B2p*B4m*B4p*G1m + B1m*B3m*B3p*B4m*B4p*G1m - B1m*B2m*B2p*B3m*B3p*B4m*B4p*G1m - B1m*B5m*B5p*G1m+ B1m*B2m*B2p*B5m*B5p*G1m 
                + B1m*B3m*B3p*B5m*B5p*G1m - B1m*B2m*B2p*B3m*B3p*B5m*B5p*G1m + B1m*B4m*B4p*B5m*B5p*G1m - B1m*B2m*B2p*B4m*B4p*B5m*B5p*G1m 
                - B1m*B3m*B3p*B4m*B4p*B5m*B5p*G1m + B1m*B2m*B2p*B3m*B3p*B4m*B4p*B5m*B5p*G1m))/Z; 
                __float128 forward_2 = ((A2p*B2p - A2p*B1m*B1p*B2p - A2p*B2p*B3m*B3p + A2p*B1m*B1p*B2p*B3m*B3p - A2p*B2p*B4m*B4p + A2p*B1m*B1p*B2p*B4m*B4p 
                + A2p*B2p*B3m*B3p*B4m*B4p - A2p*B1m*B1p*B2p*B3m*B3p*B4m*B4p - A2p*B2p*B5m*B5p+ A2p*B1m*B1p*B2p*B5m*B5p + A2p*B2p*B3m*B3p*B5m*B5p 
                - A2p*B1m*B1p*B2p*B3m*B3p*B5m*B5p + A2p*B2p*B4m*B4p*B5m*B5p - A2p*B1m*B1p*B2p*B4m*B4p*B5m*B5p - A2p*B2p*B3m*B3p*B4m*B4p*B5m*B5p 
                + A2p*B1m*B1p*B2p*B3m*B3p*B4m*B4p*B5m*B5p)*G2p)/Z; 
                __float128 forward_3 = ((A3p*B3p - A3p*B1m*B1p*B3p - A3p*B2m*B2p*B3p + A3p*B1m*B1p*B2m*B2p*B3p - A3p*B3p*B4m*B4p + A3p*B1m*B1p*B3p*B4m*B4p 
                + A3p*B2m*B2p*B3p*B4m*B4p - A3p*B1m*B1p*B2m*B2p*B3p*B4m*B4p - A3p*B3p*B5m*B5p+A3p*B1m*B1p*B3p*B5m*B5p + A3p*B2m*B2p*B3p*B5m*B5p 
                - A3p*B1m*B1p*B2m*B2p*B3p*B5m*B5p + A3p*B3p*B4m*B4p*B5m*B5p - A3p*B1m*B1p*B3p*B4m*B4p*B5m*B5p - A3p*B2m*B2p*B3p*B4m*B4p*B5m*B5p 
                + A3p*B1m*B1p*B2m*B2p*B3p*B4m*B4p*B5m*B5p)*G3p)/Z; 
                __float128 forward_4 =((A4p*B4p - A4p*B1m*B1p*B4p - A4p*B2m*B2p*B4p + A4p*B1m*B1p*B2m*B2p*B4p - A4p*B3m*B3p*B4p + A4p*B1m*B1p*B3m*B3p*B4p 
                + A4p*B2m*B2p*B3m*B3p*B4p - A4p*B1m*B1p*B2m*B2p*B3m*B3p*B4p - A4p*B4p*B5m*B5p+A4p*B1m*B1p*B4p*B5m*B5p + A4p*B2m*B2p*B4p*B5m*B5p 
                - A4p*B1m*B1p*B2m*B2p*B4p*B5m*B5p + A4p*B3m*B3p*B4p*B5m*B5p - A4p*B1m*B1p*B3m*B3p*B4p*B5m*B5p - A4p*B2m*B2p*B3m*B3p*B4p*B5m*B5p 
                + A4p*B1m*B1p*B2m*B2p*B3m*B3p*B4p*B5m*B5p)*G4p)/Z; 
                __float128 forward_5 = ((A5p*B5p - A5p*B1m*B1p*B5p - A5p*B2m*B2p*B5p + A5p*B1m*B1p*B2m*B2p*B5p - A5p*B3m*B3p*B5p + A5p*B1m*B1p*B3m*B3p*B5p 
                + A5p*B2m*B2p*B3m*B3p*B5p - A5p*B1m*B1p*B2m*B2p*B3m*B3p*B5p - A5p*B4m*B4p*B5p+A5p*B1m*B1p*B4m*B4p*B5p + A5p*B2m*B2p*B4m*B4p*B5p 
                - A5p*B1m*B1p*B2m*B2p*B4m*B4p*B5p + A5p*B3m*B3p*B4m*B4p*B5p - A5p*B1m*B1p*B3m*B3p*B4m*B4p*B5p - A5p*B2m*B2p*B3m*B3p*B4m*B4p*B5p 
                + A5p*B1m*B1p*B2m*B2p*B3m*B3p*B4m*B4p*B5p)*G5p)/Z; 

                
                double random_number1 = distribution(generator);   
                if (mRNA_string.size() ==2){
                    double P_forward_2 = forward_2 /(forward_2 + forward_3 + forward_4 + forward_5 );
                    double P_forward_3 = forward_3 /(forward_2 + forward_3 + forward_4 + forward_5 );
                    double P_forward_4 = forward_4 /(forward_2 + forward_3 + forward_4 + forward_5 );
                    double P_forward_5 = forward_5 /(forward_2 + forward_3 + forward_4 + forward_5 );
                    if (random_number1 < P_forward_2){
                            mRNA_string.push_back(purine_1_int);
                            monomer_status = purine_1;
                            error_string.push_back(error(mRNA_string[Q+1], mRNA_template[Q+1]));
                        }
                    else if (random_number1 < (P_forward_2 + P_forward_3)){
                            mRNA_string.push_back(pyrimidines_1_int);
                            monomer_status = pyrimidines_1;
                            error_string.push_back(error(mRNA_string[Q+1], mRNA_template[Q+1]));
                        }
                    else if (random_number1 < (P_forward_2 + P_forward_3 + P_forward_4)){
                        mRNA_string.push_back(purine_2_int);
                        monomer_status = purine_2;
                        error_string.push_back(error(mRNA_string[Q+1], mRNA_template[Q+1]));
                        }
                    else{ //if (random_number1 < (P_forward_2 + P_forward_3 + P_forward_4 + P_forward_5)){
                        mRNA_string.push_back(pyrimidines_2_int);
                        monomer_status = pyrimidines_2;
                        error_string.push_back(error(mRNA_string[Q+1], mRNA_template[Q+1]));
                        }
                }
                else{
                    //hier ben ik
                    double P_backward = backward/(forward_2 + forward_3 + forward_4 + forward_5 + backward);
                    double P_forward_2 = forward_2 /(forward_2 + forward_3 + forward_4 + forward_5 + backward);
                    double P_forward_3 = forward_3 /(forward_2 + forward_3 + forward_4 + forward_5 + backward);
                    double P_forward_4 = forward_4 /(forward_2 + forward_3 + forward_4 + forward_5 + backward);
                    double P_forward_5 = forward_5 /(forward_2 + forward_3 + forward_4 + forward_5 + backward);
                    if (random_number1 < P_forward_2){
                            mRNA_string.push_back(purine_1_int);
                            monomer_status = purine_1;
                            error_string.push_back(error(mRNA_string[Q+1], mRNA_template[Q+1]));
                        }
                    else if (random_number1 < (P_forward_2 + P_forward_3)){
                            mRNA_string.push_back(pyrimidines_1_int);
                            monomer_status = pyrimidines_1;
                            error_string.push_back(error(mRNA_string[Q+1], mRNA_template[Q+1]));
                        }
                    else if (random_number1 < (P_forward_2 + P_forward_3 + P_forward_4)){
                        mRNA_string.push_back(purine_2_int);
                        monomer_status = purine_2;
                        error_string.push_back(error(mRNA_string[Q+1], mRNA_template[Q+1]));
                        }
                    else if (random_number1 < (P_forward_2 + P_forward_3 + P_forward_4 + P_forward_5)){
                        mRNA_string.push_back(pyrimidines_2_int);
                        monomer_status = pyrimidines_2;
                        error_string.push_back(error(mRNA_string[Q+1], mRNA_template[Q+1]));
                        }
                    else{
                        mRNA_string.pop_back();
                        monomer_status = "-";
                        error_string.pop_back();
                    }
                }
            
            if(monomer_status =="-"){

            }
            
            error_prob = count(error_string.begin(), error_string.end(), -5)/(double)error_string.size();
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
            error_prob = count(error_string.begin()+2, error_string.end()-15, -5)/((double)error_string.size()-17);
            std::vector<int> in = mRNA_string;
            // loop through the array elements
            for(int j=0; j< in.size(); ++j){
                fss << in.at(j) ;
                }
            fss << "," << error_prob << "\n";
        }
        mRNA_string_out = mRNA_string;
        std::vector<int> probabilities_bases = base_probs(mRNA_template, mRNA_string_out, purine_1_int, purine_2_int, pyrimidines_1_int, pyrimidines_2_int);
    if(multiple_DNA == true){
    fstream fss;
    string out_file_strings = declare_polymers_string(multiple_DNA, DIR, x_def, purine_1, pyrimidines_1, purine_2, pyrimidines_2, ratio_bases_template, repulsion_00_11, delta_g_pol);
    fss.open(out_file_strings.c_str(), ios::out | ios::app);
    fss << "pur1pur1" << "," << "pur1pur2" << "," << "pur2pur2" <<","<< "py1py1" <<","<< "py1py2"<<","<< "py2py2"<<"," <<"py1pur1" <<","<< "py1pur2" <<"," <<"py2pur1" <<","<< "py2pur2" << "," << "nothing" <<"\n";
    fss << "template string:" ;
    std::vector<int> in4 = mRNA_template;
    // loop through the array elements
    for(int j=0; j< in4.size(); ++j){
        fss << in4.at(j);
        fss << " ";}       
    fss << ".\n";
    fss << "copy polymer string:" ;
    std::vector<int> in7 = mRNA_string_out;
    // loop through the array elements
    for(int j=0; j< in7.size(); ++j){
        fss << in7.at(j);
        fss << " ";}       
    fss << ".\n";

    std::vector<int> in8 = probabilities_bases;
    // loop through the array elements
    for(int j=0; j< in8.size(); ++j){
        fss << in8.at(j);
        fss << ",";}       
    fss << ".\n";
    fss.close();}
        }
    //save information for all steps, error info
    if(multiple_DNA == true){
        fss << std ::endl;
        fss.close();
    }
    
    }
}
