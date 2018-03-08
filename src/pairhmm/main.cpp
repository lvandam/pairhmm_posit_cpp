/**
    @author Laurens van Dam
    @date 02/02/2018
    @copyright 2018 All rights reserved.
**/

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <posit/posit>
#include <float/quire.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "inputreader.hpp"
#include "pairhmm.hpp"
#include "pairhmm_posit.hpp"
#include "utils.hpp"

#define POSIT_NBITS 32
#define POSIT_EBITS 2

using namespace std;
using boost::multiprecision::cpp_dec_float_50;

typedef sw::unum::posit<POSIT_NBITS, POSIT_EBITS> Posit;
typedef sw::unum::quire<32, 8> QuireFloat;
typedef sw::unum::quire<80, 15> QuireDecimal50;
typedef sw::unum::quire<POSIT_NBITS, POSIT_EBITS> QuirePosit;

int main(int argc, char *argv[]) {
    cout.precision(50);
    cout.flags(cout.fixed);

    InputReader reader {};

    std::vector<cpp_dec_float_50> results_decimal, results_float, results_posit;

    const long double initial_constant = ldexpf(1.0f, 100);

    PairHMM<cpp_dec_float_50, QuireDecimal50> pairhmm_dec50(initial_constant);
    PairHMM<float, QuireFloat> pairhmm_float(initial_constant);
    PairHMMPosit<Posit, QuirePosit> pairhmm_posit(initial_constant);

    // Perform calculation per testcase for various number representations
    std::vector<Testcase> testcases = reader.from_file(argv[1]);
    for(Testcase testcase : testcases)
    {
        cpp_dec_float_50 result;

        result = pairhmm_dec50.compute_full_prob(&testcase);
        results_decimal.push_back(result);
        cout << "-- DECIMAL50 -- " << result << " -- log10 = " << log10(result) << endl;
        pairhmm_dec50.debug_values.exportDebugValues("pairhmm_decimal50.txt");

        result = pairhmm_float.compute_full_prob(&testcase);
        results_float.push_back(result);
        cout << "-- FLOAT -- = " << result << " -- log10 = " << log10(result) << endl;
        pairhmm_float.debug_values.exportDebugValues("pairhmm_float.txt");

        result = pairhmm_posit.compute_full_prob(&testcase);
        results_posit.push_back(result);
        cout << "-- POSIT<32,2> -- = " << result << " -- log10 = " << log10(result) << endl;
        pairhmm_posit.debug_values.exportDebugValues("pairhmm_posit.txt");

        // For benchmarking (print without labels)
        writeBenchmark(pairhmm_dec50, pairhmm_float, pairhmm_posit);

//        pairhmm_dec50.debug_values.printDebugValues();
//        pairhmm_float.debug_values.printDebugValues();
//        pairhmm_posit.debug_values.printDebugValues();
    }

    return 0;
}

