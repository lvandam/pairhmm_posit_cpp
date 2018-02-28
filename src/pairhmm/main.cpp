/**
    @author Laurens van Dam
    @date 02/02/2018
    @copyright 2018 All rights reserved.
**/

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <boost/range/combine.hpp>

#include "inputreader.hpp"
#include "pairhmm.hpp"

using namespace std;

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

int main(int argc, char *argv[]) {
    cout.precision(50);
    cout.flags(cout.fixed);

    InputReader reader {};


    std::vector<float> results_ld, results_float, results_posit;

    PairHMM<long double> pairhmm_ld;
    PairHMM<float> pairhmm_float;
    PairHMM<posit<32,2>> pairhmm_posit;

    // Perform calculation per testcase for various number representations
    std::vector<Testcase> testcases = reader.from_file(argv[1]);
    for(Testcase testcase : testcases)
    {
        float result;

        result = pairhmm_ld.compute_full_prob(&testcase);
        results_ld.push_back(result);
        cout << "-- LONG DOUBLE -- " << result << " -- log10 = " << log10(result) << endl;
//        pairhmm_ld.debug_values.printDebugValues();
        pairhmm_ld.debug_values.exportDebugValues("pairhmm_ld.txt");

        result = pairhmm_float.compute_full_prob(&testcase);
        results_float.push_back(result);
        cout << "-- FLOAT -- = " << result << " -- log10 = " << log10(result) << endl;
//        pairhmm_float.debug_values.printDebugValues();
        pairhmm_float.debug_values.exportDebugValues("pairhmm_float.txt");

        result = pairhmm_posit.compute_full_prob(&testcase);
        results_posit.push_back(result);
        cout << "-- POSIT<32,2> -- = " << result << " -- log10 = " << log10(result) << endl;
//        pairhmm_posit.debug_values.printDebugValues();
        pairhmm_posit.debug_values.exportDebugValues("pairhmm_posit.txt");

        // For benchmarking (print without labels)
        writeBenchmark(pairhmm_ld, pairhmm_float, pairhmm_posit);
    }

    return 0;
}

