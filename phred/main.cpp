#include <iostream>

#include <posit/posit>
#include <boost/range/combine.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;
using namespace sw::unum;
using boost::multiprecision::cpp_dec_float_100;

cpp_dec_float_100 decimal_accuracy(cpp_dec_float_100 exact, cpp_dec_float_100 computed) {
    if (boost::math::isnan(exact) || boost::math::isnan(computed) ||
        (boost::math::sign(exact) != boost::math::sign(computed))) {
        return std::numeric_limits<cpp_dec_float_100>::quiet_NaN();
    } else if (exact == computed) {
        return std::numeric_limits<cpp_dec_float_100>::infinity();
    } else if ((exact == std::numeric_limits<cpp_dec_float_100>::infinity() &&
                computed != std::numeric_limits<cpp_dec_float_100>::infinity()) ||
               (exact != std::numeric_limits<cpp_dec_float_100>::infinity() &&
                computed == std::numeric_limits<cpp_dec_float_100>::infinity()) || (exact == 0 && computed != 0) ||
               (exact != 0 && computed == 0)) {
        return -std::numeric_limits<cpp_dec_float_100>::infinity();
    } else {
        return -log10(abs(log10(computed / exact)));
    }
}

int main() {

    ofstream outfile("da.txt", ios::out);
    outfile << "Q,da_f,da_p2,da_p3" << endl;

    for(int Q = 1; Q <= 100; Q++) {
        cpp_dec_float_100 dec = pow(10.0, -(cpp_dec_float_100)Q/10);

        float f = powl(10.0, -(long double)Q/10);
        posit<32,2> p2 = powl(10.0, -(long double)Q/10);
        posit<32,3> p3 = powl(10.0, -(long double)Q/10);

        cpp_dec_float_100 da_f, da_p2, da_p3;
        da_f = decimal_accuracy(dec, static_cast<cpp_dec_float_100>(f));
        da_p2 = decimal_accuracy(dec, static_cast<cpp_dec_float_100>(p2));
        da_p3 = decimal_accuracy(dec, static_cast<cpp_dec_float_100>(p3));

        outfile << Q << ",";
        outfile << setprecision(100) << fixed << da_f << "," << da_p2 << "," << da_p3 << endl << flush;
    }
    outfile.close();

}