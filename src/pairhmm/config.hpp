/**
    @author Laurens van Dam
    @date 19/02/2018
    @copyright 2018 All rights reserved.
**/

#ifndef PAIRHMM_CONFIG_HPP
#define PAIRHMM_CONFIG_HPP

#include <posit/posit>
#include <float/quire.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#ifndef DEBUG_PRECISION
#define DEBUG_PRECISION 40
#endif // DEBUG_PRECISION

#define POSIT_NBITS 32
#define POSIT_EBITS 2

using namespace std;
using boost::multiprecision::cpp_dec_float_50;

typedef sw::unum::posit<POSIT_NBITS, POSIT_EBITS> Posit;
typedef sw::unum::quire<32, 8> QuireFloat;
typedef sw::unum::quire<80, 15> QuireDecimal50;
typedef sw::unum::quire<POSIT_NBITS, POSIT_EBITS> QuirePosit;

struct Entry {
    string name;
    cpp_dec_float_50 value;
};

struct find_entry {
    string name;
    find_entry(string name) : name(name) {}
    bool operator () (const Entry& m) const
    {
        return m.name == name;
    }
};

#endif //PAIRHMM_CONFIG_HPP
