/**
    @author Laurens van Dam
    @date 02/02/2018
    @copyright 2018 All rights reserved.
**/

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

#include "inputreader.hpp"
#include "pairhmm_ld.hpp"
#include "pairhmm_float.hpp"
#include "pairhmm_posit.hpp"

using namespace std;

int main(int argc, char *argv[]) {

    InputReader reader {};
    PairHMMLongDouble pairhmm_ld {};
    PairHMMFloat pairhmm_float {};
    PairHMMPosit pairhmm_posit {};

    std::vector<Testcase> testcases = reader.from_file(argv[1]);

    // Perform calculation per testcase for various number representations
    std::vector<float> results_ld;
    std::vector<float> results_float;
    std::vector<float> results_posit;
    for(Testcase testcase : testcases)
    {
        printDebug(">>> LONG DOUBLE");
        results_ld.push_back(pairhmm_ld.compute_full_prob(&testcase));
        printDebug(">>> FLOAT");
        results_float.push_back(pairhmm_float.compute_full_prob(&testcase));
        printDebug("\n>>> POSIT");
        results_posit.push_back(pairhmm_posit.compute_full_prob(&testcase));
    }
    cout << endl;

    // Results
    cout << "RESULTS" << endl;
    for(int i = 0; i < testcases.size(); i++) {
        printDebug("Read #%d", i);

        cout << "-- LONG DOUBLE = " << fixed << setprecision(50) << results_ld[i] << " -- log10 = " << log10(results_ld[i]) << endl;
        cout << "-- FLOAT = " << fixed << setprecision(50) << results_float[i] << " -- log10 = " << log10(results_float[i]) << endl;
        cout << "-- POSIT = " << fixed << setprecision(50) << results_posit[i] << " -- log10 = " << log10(results_posit[i]) << endl;
    }
    cout << endl;

    // Print intermediate debug values
    cout << "DEBUG VALUES" << endl;
    cout << ">> LONG DOUBLE" << endl;
    pairhmm_ld.debug_values.printDebugValues(); pairhmm_ld.debug_values.exportDebugValues("pairhmm_ld.txt");
    cout << ">> FLOAT" << endl;
    pairhmm_float.debug_values.printDebugValues(); pairhmm_float.debug_values.exportDebugValues("pairhmm_float.txt");
    cout << ">> POSIT" << endl;
    pairhmm_posit.debug_values.printDebugValues(); pairhmm_posit.debug_values.exportDebugValues("pairhmm_posit.txt");

    return 0;
}

