/**
    @author Laurens van Dam
    @date 19/02/2018
    @copyright 2018 All rights reserved.
**/

#ifndef PAIRHMM_UTILS_HPP
#define PAIRHMM_UTILS_HPP

#include <iostream>
#include <iomanip>
#include <boost/range/combine.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "config.hpp"
#include "pairhmm.hpp"
#include "pairhmm_posit.hpp"

using namespace std;
using boost::multiprecision::cpp_dec_float_50;

void writeBenchmarkText(const char *format, ...) {
    char buf[1024];

    va_list arglist;
    va_start(arglist, format);
    vsprintf(buf, format, arglist);
    va_end(arglist);

    ofstream outfile("pairhmm_values.txt", ios::out|ios::app);
    outfile << buf << endl;
    outfile.close();
}

void writeBenchmark(PairHMM<auto, auto>& pairhmm_dec50, PairHMM<auto, auto>& pairhmm_float, PairHMMPosit<auto, auto>& pairhmm_posit, std::string filename = "pairhmm_values.txt", bool printDate = true, bool overwrite = false) {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    ofstream outfile;
    if(overwrite) {
        outfile = ofstream(filename, ios::out);
    }
    else {
        outfile = ofstream(filename, ios::out|ios::app);
    }

    if(printDate)
        outfile << endl << put_time(&tm, "%d-%m-%Y %H:%M:%S") << endl << "===================" << endl;

    auto names = pairhmm_dec50.debug_values.getNames();
    auto dec_values = pairhmm_dec50.debug_values.getValues();
    auto float_values = pairhmm_float.debug_values.getValues();
    auto posit_values = pairhmm_posit.debug_values.getValues();

    outfile << "name,dE_f,dE_p,log(abs(dE_f)),log(abs(dE_p))" << endl;
    for(auto tup : boost::combine(names, dec_values, float_values, posit_values)) {
        string name;
        cpp_dec_float_50 E, dE_f, E_f, dE_p, E_p;

        boost::tie(name, E, E_f, E_p) = tup;

        if(E == 0) {
            dE_f = 0; dE_p = 0;
        } else {
            dE_f = (E_f - E) / E;
            dE_p = (E_p - E) / E;
        }

        // Relative error values
        outfile << setprecision(50) << fixed << name <<","<< dE_f <<","<< dE_p <<","<< log10(abs(dE_f)) <<","<< log10(abs(dE_p)) << endl;
    }
    outfile.close();
}

#endif //PAIRHMM_UTILS_HPP
