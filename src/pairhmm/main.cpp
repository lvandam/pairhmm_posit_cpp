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

using namespace std;

int main(int argc, char *argv[]) {
    cout.precision(50);
    cout.flags(cout.fixed);

    InputReader reader {};

    float result;

    std::vector<float> results_ld;
    PairHMM<long double> pairhmm_ld;

    std::vector<float> results_float;
    PairHMM<float> pairhmm_float;

    std::vector<float> results_posit;
    PairHMM<posit<32,2>> pairhmm_posit;

    // Perform calculation per testcase for various number representations
    std::vector<Testcase> testcases = reader.from_file(argv[1]);
    for(Testcase testcase : testcases)
    {
        result = pairhmm_ld.compute_full_prob(&testcase);
        results_ld.push_back(result);
        cout << "-- LONG DOUBLE -- " << result << " -- log10 = " << log10(result) << endl;
        pairhmm_ld.debug_values.printDebugValues(); pairhmm_ld.debug_values.exportDebugValues("pairhmm_ld.txt");

        result = pairhmm_float.compute_full_prob(&testcase);
        results_float.push_back(result);
        cout << "-- FLOAT -- = " << result << " -- log10 = " << log10(result) << endl;
        pairhmm_float.debug_values.printDebugValues(); pairhmm_float.debug_values.exportDebugValues("pairhmm_float.txt");

        result = pairhmm_posit.compute_full_prob(&testcase);
        results_posit.push_back(result);
        cout << "-- POSIT<32,2> -- = " << result << " -- log10 = " << log10(result) << endl;
        pairhmm_posit.debug_values.printDebugValues(); pairhmm_posit.debug_values.exportDebugValues("pairhmm_posit.txt");
    }

    return 0;
}

