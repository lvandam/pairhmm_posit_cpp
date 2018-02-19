//
// Created by Laurens van Dam on 19/02/2018.
//

#ifndef PAIRHMM_SIMPLE_UTILS_HPP
#define PAIRHMM_SIMPLE_UTILS_HPP

#define PRECISION 40

#include <iostream>
#include <iomanip>

using namespace std;

void printDebug(const char* format, ...) {
#ifdef PRINTVALUES
    char buf[1024];

    va_list arglist;
    va_start(arglist, format);
    vsprintf(buf, format ,arglist);
    va_end(arglist);

    cout << buf << endl;
#endif
}

template<class T>
void printValue(T &value) {
#ifdef PRINTVALUES
    cout << fixed << setprecision(PRECISION) << value << endl;
#endif
}

template<class T>
void printDebugValue(T &value, string name) {
#ifdef PRINTVALUES
    cout << name << " = " << fixed << setprecision(PRECISION) << value << endl;
#endif
}

#endif //PAIRHMM_SIMPLE_UTILS_HPP
