/**
    @author Laurens van Dam
    @date 21/02/2018
    @copyright 2018 All rights reserved.
**/

#include <iostream>
#include <posit/posit>
#include <vector>
#include <cmath>

#include <gnuplot-iostream.h>

#define NBITS 32
#define ES 2

using namespace std;
using namespace sw::unum;

long double score_to_probability(int Q) {
    return powl(10, -static_cast<long double>(Q) / 10);
}

int main() {
    Gnuplot gp;

    int Q_min = 1;
    int Q_max = 100;

    vector<tuple<int, long double, long double>> results;
    for(int Q = Q_min; Q <= Q_max; Q++) {
        long double E = score_to_probability(Q);
        auto E_f = static_cast<float>(E);
        posit<32,2> E_p(E);

        results.push_back(make_tuple(Q, ((E_f - E) / E), ((static_cast<long double>(E_p) - E) / E)));
    }

    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    ofstream outfile("phred_float_posit.txt", ios::out|ios::app);
    outfile << endl << put_time(&tm, "%d-%m-%Y %H:%M:%S") << endl << "===================" << endl;

    cout << "Q,dE_f,log(abs(dE_f)),dE_p,log(abs(dE_p))" << endl;
    outfile << "Q,dE_f,log(abs(dE_f)),dE_p,log(abs(dE_p))" << endl;
    for(auto el : results) {
        cout << setprecision(40) << fixed << get<0>(el) <<","<< get<1>(el) <<","<< log10l(abs(get<1>(el))) <<","<< get<2>(el) <<","<< log10l(abs(get<2>(el))) << endl;
        outfile << setprecision(40) << fixed << get<0>(el) <<","<< get<1>(el) <<","<< log10l(abs(get<1>(el))) <<","<< get<2>(el) <<","<< log10l(abs(get<2>(el))) << endl;
    }
    outfile.close();

    return 0;

//    gp << "set grid mxtics xtics\n";
//    gp << "set xtics 1\nset mxtics 10\n";
//    gp << "set xlabel 'Phred Quality Score'\nset ylabel 'Relative error'\n";
//    gp << "set style fill solid\nset style circle radius 0.4\nset xrange [1:100]\nset yrange [-1e-7:1e-7]\n";
//    gp << "set arrow 1 from 0,0 to 100,0 nohead\n";
//    gp << "plot '-' with circles title 'dE_f', '-' with circles title 'dE_p'\n";
//    gp.send1d(boost::make_tuple(vec_Q, dE_f));
//    gp.send1d(boost::make_tuple(vec_Q, dE_p));
//    gp << "replot\n";
}