/**
    @author Laurens van Dam
    @date 21/02/2018
    @copyright 2018 All rights reserved.
**/

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <posit/posit>
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;
using namespace sw::unum;
using boost::multiprecision::cpp_dec_float_50;

cpp_dec_float_50 decimal_accuracy(cpp_dec_float_50 exact, cpp_dec_float_50 computed) {
    if(isnan(exact) || isnan(computed) || (sign(exact) != sign(computed))) {
        return std::numeric_limits<cpp_dec_float_50>::quiet_NaN();
    } else if(exact == computed) {
        return std::numeric_limits<cpp_dec_float_50>::infinity();
    } else if((exact == std::numeric_limits<cpp_dec_float_50>::infinity() && computed != std::numeric_limits<cpp_dec_float_50>::infinity()) || (exact != std::numeric_limits<cpp_dec_float_50>::infinity() && computed == std::numeric_limits<cpp_dec_float_50>::infinity()) || (exact == 0 && computed != 0) || (exact != 0 && computed == 0)) {
        return -std::numeric_limits<cpp_dec_float_50>::infinity();
    } else {
        return -log10(abs(log10(computed/exact)));
    }
}

cpp_dec_float_50 score_to_probability(int Q) {
    return boost::multiprecision::pow(static_cast<cpp_dec_float_50>(10), -Q / 10);
}

int main() {
    ofstream outfile("phred_float_posit.txt", ios::out);
    outfile << "Q,da_f,da_p2,da_p3" << endl;

    for(int Q = 1; Q <= 100; Q++) {
        cpp_dec_float_50 E = score_to_probability(Q);
        float E_f = static_cast<float>(E);
        posit<32,2> E_p2(E);
        posit<32,3> E_p3(E);

        cpp_dec_float_50 da_f = decimal_accuracy(E, E_f);
        cpp_dec_float_50 da_p2 = decimal_accuracy(E, static_cast<cpp_dec_float_50>(E_p2));
        cpp_dec_float_50 da_p3 = decimal_accuracy(E, static_cast<cpp_dec_float_50>(E_p3));

        outfile << setprecision(50) << fixed << Q <<","<< da_f <<","<< da_p2 <<","<< da_p3 << endl;
    }

    outfile.close();

    return 0;
}