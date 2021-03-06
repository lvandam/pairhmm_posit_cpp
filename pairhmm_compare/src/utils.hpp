#ifndef __UTILS_H
#define __UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <posit/posit>

#include "debug_values.hpp"
#include "pairhmm_float.hpp"
#include "pairhmm_posit.hpp"
#include "defines.hpp"
#include "batch.hpp"

using namespace std;
using namespace sw::unum;

cpp_dec_float_100 decimal_accuracy(cpp_dec_float_100 exact, cpp_dec_float_100 computed);

void writeBenchmark(PairHMMFloat<cpp_dec_float_100> &pairhmm_dec50, PairHMMFloat<float> &pairhmm_float,
                    PairHMMPosit &pairhmm_posit, std::string filename = "pairhmm_values.txt", bool printDate = true, bool overwrite = false);

void print_batch_info(t_batch& batch);

int px(int x, int y);

int pbp(int x);

int py(int y);

t_workload *gen_workload(unsigned long pairs, unsigned long fixedX, unsigned long fixedY);

void copyProbBytes(t_probs& probs, uint8_t bytesArray[]);

int roundToMultiple(int toRound, int multiple);

std::string randomBasepairs(int len);

#endif //__UTILS_H
