//
// Created by Laurens van Dam on 19/02/2018.
//

#ifndef PAIRHMM_SIMPLE_INTERMEDIATE_HPP
#define PAIRHMM_SIMPLE_INTERMEDIATE_HPP

#include <iostream>
#include <iomanip>
#include <vector>
#include "config.hpp"

template<class T>
class Intermediate {
private:
    std::vector<T> items;

public:
    Intermediate() = default;

    void debugValue(T value, const char* format, ...) {
#ifdef DEBUG_VALUES
        char buf[1024];

        va_list arglist;
        va_start(arglist, format);
        vsprintf(buf, format ,arglist);
        va_end(arglist);

        items.push_back(value);

        cout << buf << " = " << fixed << setprecision(DEBUG_PRECISION) << value << endl;
#endif
    }
};

#endif //PAIRHMM_SIMPLE_INTERMEDIATE_HPP
