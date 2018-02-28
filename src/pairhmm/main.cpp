/**
    @author Laurens van Dam
    @date 02/02/2018
    @copyright 2018 All rights reserved.
**/

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

#include "inputreader.hpp"
#include "pairhmm.hpp"
#include "utils.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    cout.precision(50);
    cout.flags(cout.fixed);

    InputReader reader {};

    std::vector<float> results_ld, results_float, results_posit;

    const long double initial_constant = ldexpf(1.0f, 5);

    PairHMM<long double> pairhmm_ld(initial_constant);
    PairHMM<float> pairhmm_float(initial_constant);
    PairHMM<posit<32,2>> pairhmm_posit(initial_constant);

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

