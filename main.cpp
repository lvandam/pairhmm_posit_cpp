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

    // Perform calculation per testcase
    std::vector<double> results_double;
    std::vector<double> results_posit;
    for(Testcase testcase : testcases)
    {
        results_double.push_back(pairhmm_float.compute_full_prob(&testcase));
        results_posit.push_back(pairhmm_posit.compute_full_prob(&testcase));
    }

    for(int i = 0; i < testcases.size(); i++) {
        cout << "Read #" << i << " -- Result double = " << setprecision(20) << log10(results_double[i]) << ", result posit = " << log10(results_posit[i]) << endl;
    }

    return 0;
}

