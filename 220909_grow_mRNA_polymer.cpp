/*
09-09-22, Muriel Louman, AMOLF

Grow a polymer string, without keeping track what you are adding. 
There is no wrong or right when adding a monomer.
*/


#include <iostream>
#include <cmath> // math functions
#include <vector> // vector container
#include <string> // text strings
#include <fstream> // output to file
#include <sstream> 
#include <random> // random number generator

// define rates
double a_1_func(double k_on, double delta_g_r){
    return k_on * exp(-1 * delta_g_r);}

double a_2_func(double k_on, double r_con){
    return k_on * r_con;}

double b_1_func(double k_bb, double delta_g_bb){
    return k_bb * exp(-1 * delta_g_bb);}

double b_2_func(double k_bb){
    return k_bb;}

double c_1_func(double k_on, double r_con_star){
    return k_on * r_con_star;}

double c_2_func(double k_on, double delta_g_r){
    return k_on * exp(-1* delta_g_r);}


// define probabilities between states instead of rates
double prob_alpha_1(double a_1, double b_2){
    return a_1/(a_1 + b_2);}

double prob_beta_2(double a_1, double b_2){
    return b_2/(a_1 + b_2);}

double prob_alpha_2(double a_2, double c_1){
    return a_2/(a_2 + c_1);}

double prob_gamma_1(double a_2, double c_1){
    return c_1/(a_2 + c_1);}

double prob_beta_1(double b_1, double c_2){
    return b_1/(c_2 + b_1);}

double prob_gamma_2(double b_1, double c_2){
    return c_2/(c_2 + b_1);}


void print_matrix(std::vector<int> input){
    int n = sizeof(input)/sizeof(input[0]);
 
    // loop through the array elements
    for(int i=0; i< input.size(); ++i){
        std::cout << input.at(i) << ' ';
        }
}

int main(){
    // define variables
    double k_on = 1.0;
    double k_bb = 1.0;
    double r_con = 1.0;
    double r_con_star = 1.5;
    double delta_g_r = 0;
    double delta_g_bb = -4;
    double a_1 = a_1_func(k_on, delta_g_r);
    double a_2 = a_2_func(k_on, r_con);
    double b_1 = b_1_func(k_bb, delta_g_bb);
    double b_2 = b_2_func(k_bb);
    double c_1 = c_1_func(k_on, r_con_star);
    double c_2 = c_2_func(k_on, delta_g_r);
    double alpha_1 = prob_alpha_1(a_1, b_2);
    double alpha_2 = prob_alpha_2(a_2, c_1);
    double beta_1 = prob_beta_1(b_1, c_2);
    double beta_2 = prob_beta_2(a_1, b_2);
    double gamma_1 = prob_gamma_1(a_2, c_1);
    double gamma_2 = prob_gamma_2(b_1, c_2);
    // double alpha_1 = 0.2;
    // double alpha_2 = 0.8;
    // double beta_1 = 0.2;
    // double beta_2 = 0.8;
    // double gamma_1 = 0.2;
    // double gamma_2= 0.8;
    int length_mRNA = 10;
    std::vector<int> mRNA_string;
    mRNA_string.push_back(1);
    int i = 1;
    // print_matrix(mRNA_string);

    // random number generator for uniform distr between 0 and 1
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    
    while(mRNA_string[i-1]< length_mRNA){
        // Go to transition state 1 or 2, first MC step
        // create a random number betwen 0 and 1
        double random_number1 = distribution(generator);
        // std::cout << random_number << std::endl;

        if(random_number1< gamma_1 || mRNA_string.at(i-1)==1){
            //go to transition state 1, when the state from where you come is 1 then you always go to transition state 1
            mRNA_string.push_back(-1);}
        else{
            //go to transition state 2
            mRNA_string.push_back(-2);}
        ++i;

        // if in transition state 1
        if(mRNA_string[i-1] == -1){
            double random_number2 = distribution(generator);
            double prob1 = (beta_1 * gamma_2)/(1 - beta_1*beta_2);
            if(random_number2 < prob1){
                //add nucleoide and keep track of how much nucleoide you added
                mRNA_string.push_back(mRNA_string.at(i-2) + 1);
                std::cout << "moved to a new state"  << std::endl;
                std::cout << mRNA_string.at(i-2) + 1  << std::endl;
            }
            else{
                //do not add a nucleoide (go back to the number of nucleoide you had before going to transition state 1)
                mRNA_string.push_back(mRNA_string.at(i-2));
            }
            
        }

        // if in transition state 2
        if(mRNA_string[i-1] == -2){
            double random_number3 = distribution(generator);
            double prob2 = (alpha_1 * beta_1)/(1 - beta_1*beta_2);
            if(random_number3 < prob2){
                //remove nucleoide and keep track of how much nucleoide you added
                mRNA_string.push_back(mRNA_string.at(i-2) - 1);
            }
            else{
                //do not remove a nucleoide (go back to the number of nucleoide you had before going to transition state 2)
                mRNA_string.push_back(mRNA_string.at(i-2));  
            }
        }
        ++i;
     
        // print_matrix(mRNA_string);
        // std::cout << "new computation"  << std::endl;
    }
print_matrix(mRNA_string);
}

