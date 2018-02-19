#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

#include "inputreader.hpp"
#include "pairhmm_float.hpp"
#include "pairhmm_posit.hpp"

#define DEBUG_PRECISION 40

using namespace std;

int main(int argc, char *argv[]) {

    InputReader reader {};
    PairHMMFloat pairhmm_float {};
    PairHMMPosit pairhmm_posit {};

    std::vector<Testcase> testcases = reader.from_file(argv[1]);

    // Perform calculation per testcase for various number representations
    std::vector<float> results_float;
    std::vector<float> results_posit;
    for(Testcase testcase : testcases)
    {
        printDebug(">>> FLOAT");
        results_float.push_back(pairhmm_float.compute_full_prob(&testcase));
        printDebug("\n>>> POSIT");
        results_posit.push_back(pairhmm_posit.compute_full_prob(&testcase));
    }

    cout << endl;
    for(int i = 0; i < testcases.size(); i++) {
        printDebug("Read #%d", i);

        double dec_error_float = abs(log10(results_float[i] / 1140165690094602607336957012851.948));
        cout << "-- Result float = " << fixed << setprecision(50) << results_float[i] << " -- log10 = " << log10(results_float[i]) << " -- decimal error: " << dec_error_float << endl;

        double dec_error_posit = abs(log10(results_posit[i] / 1140165690094602607336957012851.948));
        cout << "-- Result posit = " << fixed << setprecision(50) << results_posit[i] << " -- log10 = " << log10(results_posit[i]) << " -- decimal error: " << dec_error_posit << endl;
    }

    return 0;
}

