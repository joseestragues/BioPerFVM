#ifndef __utils_h__
#define __utils_h__
#include<vector>
#include<string>
#include<iostream>
#include<chrono>
using namespace std;

//Vector operations
class timer {
    public:
    chrono::high_resolution_clock timer;

    //Methods
    void init_register(std::string path);
    void star_clock();
    void stop_clock();
    void end_register(std::string path);

};


#endif