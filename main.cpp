#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

#include "inputreader.hpp"
#include "pairhmm_float.hpp"
#include "pairhmm_posit.hpp"

using namespace std;

int main(int argc, char *argv[]) {

    InputReader reader {};
    PairHMMFloat pairhmm_float {};
    PairHMMPosit pairhmm_posit {};

    std::vector<Testcase> testcases = reader.from_file(argv[1]);

    // Perform calculation per testcase (for double and posit)
    std::vector<double> results_double;
    std::vector<double> results_posit;
    for(Testcase testcase : testcases)
    {
        cout << "FLOAT" << endl;
        results_double.push_back(pairhmm_float.compute_full_prob(&testcase));
        cout << endl << "POSIT" << endl;
        results_posit.push_back(pairhmm_posit.compute_full_prob(&testcase));
    }

    cout << endl;
    for(int i = 0; i < testcases.size(); i++) {
        cout << "Read #" << i << ":" << endl;
        cout << "-- Result double = " << fixed << setprecision(30) << results_double[i] << " -- log10 = " << log10(results_double[i]) << endl;
        cout << "-- Result posit  = " << fixed << setprecision(30) << results_posit[i] << " -- log10 = " << log10(results_posit[i]) << endl;
    }

    return 0;
}

