/**
    @author Laurens van Dam
    @date 19/02/2018
    @copyright 2018 All rights reserved.
**/

#ifndef PAIRHMM_DEBUG_VALUES_HPP
#define PAIRHMM_DEBUG_VALUES_HPP

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <ctime>
#include "config.hpp"

using namespace std;

template<class T>
class DebugValues {
private:
    struct Entry {
        char name[1024];
        T value;
    };

    std::vector<Entry> items;

public:
    DebugValues() = default;

    template<class Q>
    void debugValue(Q val, const char* format, ...) {
        char buf[1024];

        T value = static_cast<T>(val);

        va_list arglist;
        va_start(arglist, format);
        vsprintf(buf, format ,arglist);
        va_end(arglist);

        Entry entry;
        strcpy(entry.name, buf);
        entry.value = value;
        items.push_back(entry);

#ifdef DEBUG_VALUES
        cout << buf << " = " << fixed << setprecision(DEBUG_PRECISION) << value << endl;
#endif
    }

    void printDebugValues() {
        for(auto el : items) {
            cout << setw(20) << el.name << " = " << fixed << setprecision(DEBUG_PRECISION) << el.value << endl;
        }
    }

    void exportDebugValues(string filename) {
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        ofstream outfile(filename, ios::out|ios::app);
        outfile << endl << put_time(&tm, "%d-%m-%Y %H:%M:%S") << endl << "===================" << endl;

        for(auto el : items) {
            outfile << el.name << "," << fixed << setprecision(DEBUG_PRECISION) << el.value << endl;
        }

        outfile.close();
    }

    std::vector<std::string> getNames() {
        std::vector<std::string> result;
        for(auto el : items) {
            result.push_back(el.name);
        }
        return result;
    }

    std::vector<T> getValues() {
        std::vector<T> result;
        for(auto el : items) {
            result.push_back(el.value);
        }
        return result;
    }
};

#endif //PAIRHMM_DEBUG_VALUES_HPP
