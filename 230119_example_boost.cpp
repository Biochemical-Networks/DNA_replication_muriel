// #include <boost/lambda/lambda.hpp>
// #include <iostream>
// #include <iterator>
// #include <algorithm>


// int main(){
//     // using namespace boost::lambda;
//     // typedef std::istream_iterator<int> in;

//     // std::for_each(
//     //     in(std::cin), in(), std::cout << (_1 * 3) << " " );}

#include <iostream>
#include <iterator>
#include <algorithm>

#include <boost/program_options.hpp>
namespace po = boost::program_options;


using namespace std;


int main(int ac, char* av[])
{po::options_description desc("Allowed options");
    try {

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("compression", po::value<int>(), "set compression level") //set declaration of values
        ;
        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 0;}
        if (vm.count("compression")) {
            std::cout << "Compression level was set to "
                 << vm["compression"].as<int>() << ".\n"; //print value definition if set
            // double x = vm["compression"].as<double>();
            // double y = x*10;
            // std::cout << y <<std::endl;
            } 
        else {
            std::cout << "Compression level was not set.\n";}
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;}
    catch(...) {
        cerr << "Exception of unknown type!\n";}

    return 0;
}
