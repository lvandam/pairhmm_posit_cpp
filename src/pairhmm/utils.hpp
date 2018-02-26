//
// Created by Laurens van Dam on 19/02/2018.
//

#ifndef PAIRHMM_SIMPLE_UTILS_HPP
#define PAIRHMM_SIMPLE_UTILS_HPP

#include <iostream>
#include <iomanip>
#include "config.hpp"

using namespace std;

template<typename T>
T score_to_probability(int Q) {
    return powl(10, -static_cast<T>(Q) / 10);
}

void printDebug(const char* format, ...) {
#ifdef DEBUG_VERBOSE
    char buf[1024];

    va_list arglist;
    va_start(arglist, format);
    vsprintf(buf, format, arglist);
    va_end(arglist);

    cout << fixed << setprecision(DEBUG_PRECISION) << buf << endl;
#endif
}

#endif //PAIRHMM_SIMPLE_UTILS_HPP
