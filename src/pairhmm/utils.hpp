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
#include "pairhmm_decimal.hpp"

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

cpp_dec_float_50 decimal_accuracy(cpp_dec_float_50 exact, cpp_dec_float_50 computed) {
    if(isnan(exact) || isnan(computed) || (sign(exact) != sign(computed))) {
        return std::numeric_limits<cpp_dec_float_50>::quiet_NaN();
    } else if(exact == computed) {
        return std::numeric_limits<cpp_dec_float_50>::infinity();
    } else if((exact == std::numeric_limits<cpp_dec_float_50>::infinity() && computed != std::numeric_limits<cpp_dec_float_50>::infinity()) || (exact != std::numeric_limits<cpp_dec_float_50>::infinity() && computed == std::numeric_limits<cpp_dec_float_50>::infinity()) || (exact == 0 && computed != 0) || (exact != 0 && computed == 0)) {
        return -std::numeric_limits<cpp_dec_float_50>::infinity();
    } else {
        return -log10(abs(log10(computed/exact)));
    }
}

void writeBenchmark(PairHMMDecimal<auto>& pairhmm_dec50, PairHMM<auto, auto>& pairhmm_float, PairHMMPosit<auto, auto>& pairhmm_posit, std::string filename = "pairhmm_values.txt", bool printDate = true, bool overwrite = false) {
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

    auto dec_values = pairhmm_dec50.debug_values.items;
    auto float_values = pairhmm_float.debug_values.items;
    auto posit_values = pairhmm_posit.debug_values.items;

    outfile << "name,dE_f,dE_p,log(abs(dE_f)),log(abs(dE_p)),E,E_f,E_p,da_F,da_P" << endl;
    for(int i = 0; i < dec_values.size(); i++) {
        cpp_dec_float_50 E, E_f, E_p, dE_f, dE_p;
        cpp_dec_float_50 da_F, da_P; // decimal accuracies

        string name = dec_values[i].name;
        E = dec_values[i].value;

        auto E_f_entry = std::find_if(float_values.begin(), float_values.end(), find_entry(name));
        E_f = E_f_entry->value;

        auto E_p_entry = std::find_if(posit_values.begin(), posit_values.end(), find_entry(name));
        E_p = E_p_entry->value;

        if(name != E_f_entry->name || name != E_p_entry->name) {
            cout << "Error: mismatching names! Could not find name '" << E_f_entry->name << "' in PairHMMPosit" << endl;
        }

        da_F = decimal_accuracy(E, E_f);
        da_P = decimal_accuracy(E, E_p);

        if(E == 0) {
            dE_f = 0; dE_p = 0;
        } else {
            dE_f = (E_f - E) / E;
            dE_p = (E_p - E) / E;
        }

        // Relative error values
        outfile << setprecision(50) << fixed << name <<","<< dE_f <<","<< dE_p <<","<< log10(abs(dE_f)) <<","<< log10(abs(dE_p)) <<","<< E <<","<< E_f <<","<< E_p <<","<< da_F <<","<< da_P << endl;
    }
    outfile.close();
}

#endif //PAIRHMM_UTILS_HPP
