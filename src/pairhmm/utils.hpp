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

#include "config.hpp"
#include "pairhmm.hpp"

using namespace std;

void printDebug(const char* format, ...) {
#ifdef DEBUG_VERBOSE
    char buf[1024];

    va_list arglist;
    va_start(arglist, format);
    vsprintf(buf, format, arglist);
    va_end(arglist);

    cout << fixed << setprecision(DEBUG_PRECISION) << buf << endl;
#endif
}

void writeBenchmark(PairHMM<long double>& pairhmm_ld, PairHMM<float>& pairhmm_float, PairHMM<posit<32,2>>& pairhmm_posit) {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    ofstream outfile("pairhmm_values.txt", ios::out|ios::app);
    outfile << endl << put_time(&tm, "%d-%m-%Y %H:%M:%S") << endl << "===================" << endl;

    auto names = pairhmm_ld.debug_values.getNames();
    auto ld_values = pairhmm_ld.debug_values.getValues();
    auto float_values = pairhmm_float.debug_values.getValues();
    auto posit_values = pairhmm_posit.debug_values.getValues();

    outfile << "name,dE_f,dE_p,log(abs(dE_f)),log(abs(dE_p))" << endl;
    for(auto tup : boost::combine(names, ld_values, float_values, posit_values)) {
        string name;
        long double E;
        float E_f;
        posit<32,2> E_p;

        boost::tie(name, E, E_f, E_p) = tup;

        long double dE_f = (E_f - E) / E;
        long double dE_p = (static_cast<long double>(E_p) - E) / E;

        // Relative error values
        outfile << setprecision(50) << fixed << name <<","<< dE_f <<","<< dE_p <<","<< log10l(abs(dE_f)) <<","<< log10l(abs(dE_p)) << endl;
    }
    outfile.close();
}

#endif //PAIRHMM_UTILS_HPP
