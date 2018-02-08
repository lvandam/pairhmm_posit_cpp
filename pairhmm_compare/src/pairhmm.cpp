// Copyright 2018 Delft University of Technology
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <cstdint>
#include <memory>
#include <vector>
#include <string>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

#include "pairhmm.hpp"
#include "debug_values.hpp"
#include "utils.hpp"
#include "batch.hpp"

using namespace std;

/**
 * Main function for pair HMM accelerator
 */
int main(int argc, char **argv) {
    srand(0);
    flush(cout);

//    posit<32,2> in1 = 0;
//    posit<32,2> in2 = -1;
//    posit<32,2> res;
//    cout << pretty_print(in1) << endl;
//    cout << pretty_print(in2) << endl;
//
//    res = in1 + in2;
//    cout << pretty_print(res) << endl;
//
//    exit(1);

    t_workload *workload;
    std::vector <t_batch> batches;

    unsigned long pairs, x, y = 0;
    int initial_constant_power = 1;
    bool calculate_sw = true;
    bool show_results = false;
    bool show_table = false;

    DEBUG_PRINT("Parsing input arguments...\n");
    if (argc > 4) {
        pairs = strtoul(argv[1], NULL, 0);
        x = strtoul(argv[2], NULL, 0);
        y = strtoul(argv[3], NULL, 0);
        initial_constant_power = strtoul(argv[4], NULL, 0);

        workload = gen_workload(pairs, x, y);

        BENCH_PRINT("M, ");
        BENCH_PRINT("%8d, %8d, %8d, ", workload->pairs, x, y);
    } else {
        fprintf(stderr,
                "ERROR: Correct usage is: %s <pairs> <X> <Y> <initial constant power>\n",
                "pairhmm");
        return (EXIT_FAILURE);
    }

    batches = std::vector<t_batch>(workload->batches);

    // Generate random basepair strings for reads and haplotypes
    std::string x_string = randomBasepairs(workload->batches * (px(x, y) + x - 1));
    std::string y_string = randomBasepairs(workload->batches * (py(y) + y - 1));

    for (int q = 0; q < workload->batches; q++) {
        fill_batch(batches[q], x_string, y_string, q, workload->bx[q], workload->by[q], powf(2.0, initial_constant_power));
        print_batch_info(batches[q]);
    }

    PairHMMPosit pairhmm_posit(workload, show_results, show_table);
    PairHMMFloat<float> pairhmm_float(workload, show_results, show_table);
    PairHMMFloat<cpp_dec_float_100> pairhmm_dec50(workload, show_results, show_table);

    pairhmm_posit.calculate(batches);
    pairhmm_float.calculate(batches);
    pairhmm_dec50.calculate(batches);

    writeBenchmark(pairhmm_dec50, pairhmm_float, pairhmm_posit,
                   "pairhmm_es" + std::to_string(ES) + "_" + std::to_string(pairs) +
                   "_" + std::to_string(x) + "_" + std::to_string(y) + "_" + std::to_string(initial_constant_power) +
                   ".txt",
                   false, true);

    return 0;
}
