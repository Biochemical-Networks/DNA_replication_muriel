
#include <iostream>
#include <cmath> // math functions
#include <vector> // vector container
#include <string> // text strings
#include <fstream> // output to file
#include <sstream> 

// helloworld example
// using namespace std;

// int main()
// {
//     vector<string> msg {"Hello", "C++", "World", "from", "VS Code", "and the C++ extension!"};

//     for (const string& word : msg)
//     {
//         cout << word << " ";
//     }
//     cout << endl;
// }



// #include <iostream>
// // function declaration
// int three(){
//     return 3;
//     }

// int main(){
//     int a=5;
//     a = three();
//     std::cout << a << std::endl;
//     return 0;
//     }   

// #include <iostream>
// // function to convert Fahrenheit to Celcius
// double FtoC(double F){
//     return (F-32.0)*5.0/9.0;
//     }
// int main(){
//     double degreesF = 57.0;
//     double degreesC = FtoC(degreesF);
//     std::cout << degreesF << " F is " << degreesC << " C" << std::endl;
//     degreesF = 86.0;
//     degreesC = FtoC(degreesF);
//     std::cout << degreesF << " F is " << degreesC << " C" << std::endl;
//     return 0;
//     }

//arrays
// #include <iostream>
// int main(){
//     int array[5]; // array with storage for five int variables
//     double array_b[3]={0.0,1.0,2.0}; // array of three doubles
//     // initialised to 0.0,1.0,2.0
//     // je kan geen arrays printen in C++ alleen de losse elementen 
//     return 0;
// }

// for loop and if statements
// int main(){
// for(int i=0; i<10; ++i){
//     if(i == 3 || i == 5 || i == 8){  
//         }
//     else{
//         // std::cout << i << std::endl;
//         }
//     if(i%2 ==0){
//         std::cout << i << std::endl;
//     }
//     }
// }

// making functions
// int add_numbers(int a, int b, int c){
//         int d;
//         d = a + b + c;
//         return d;
// }

// int main(){
//     int added_numbers;
//     int a = 2;
//     int b = 3; 
//     int c = 4; 

//     added_numbers = add_numbers(a,b,c);
//     std::cout << added_numbers << std::endl;
// }

// using namespace
// namespace car{
//     // namespace variables
//     int num_passengers;
//     double position;
//     double speed;
//     // namespace functions
//     double move_forward(){ return car::speed*10.0;}
// }     

// int main(){
//     // set namespace variables
//     car::position = 10.0;
//     car::speed = 30.0;
//     // use namespace function
//     car::position += car::move_forward();
//     std::cout << car::position << std::endl;
//     return 0;
// }

//using structures
// #include <iostream>
// using namespace std;
// struct bookStruct
// {
//     string bookTitle;
//     int bookPageN;
//     int bookReview;
//     float bookPrice;
//     // in the for loop you can now set your values
//     void set(string t,int pn,int r,float pr){
//         bookTitle = t;
//         bookPageN = pn;
//         bookReview = r;
//         bookPrice = pr;
//     }
// };

// ostream& operator << (ostream& os, const bookStruct& book)
// {
//    return os << "Title: " << book.bookTitle <<
//               " Pages: " << book.bookPageN << 
//                " Review: " << book.bookReview 
//               << " Price:" << book.bookPrice << endl;
// }

// int main() { 
//      bookStruct books[10];
//      for(int i =0;i<10;i++){
//          books[i].set("Book " + to_string(i),i*20,i,(float) i*54);
//      }
//      for(int i=0;i<10;i++){
//          cout << books[i];
//      }
//     return 0;
// }

// Create a structure variable called myStructure
// struct {
//   int myNum;
//   string myString;
// } myStructure;

// // Assign values to members of myStructure
// myStructure.myNum = 1;
// myStructure.myString = "Hello World!";

// // Print members of myStructure
// cout << myStructure.myNum << "\n";
// cout << myStructure.myString << "\n";


// //meaning using & after variable: when using & outside the main the variable will keep the updated value (changed by the function) 
// //instead of changing back to the old value it had in main. 
// using namespace std;

// // without using &
// void af(int g)
// {
//     g++;
//     cout<<g;
// }

// int main()
// {
//     int g = 123;
//     cout << g;
//     af(g);
//     cout << g;
//     return 0;
// }

// //using &
// void af(int& g)
// {
//     g++;
//     cout<<g;
// }

// int main()
// {
//     int g = 123;
//     cout << g;
//     af(g);
//     cout << g;
//     return 0;
// }