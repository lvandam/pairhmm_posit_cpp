//
// Created by Laurens van Dam on 21/02/2018.
//

#include <iostream>
#include <posit>

#define NBITS 32
#define ES 2

using namespace std;
using namespace sw::unum;

long double score_to_probability(int Q) {
    return powl(10, -static_cast<long double>(Q) / 10);
}

int main() {
    int Q_min = 1;
    int Q_max = 100;

    cout << "Q,E_ld,E_f,E_p" << endl;
    for(int Q = Q_min; Q <= Q_max; Q++) {
        long double E = score_to_probability(Q);
        posit<32,2> E_posit(E);

        cout << fixed << setprecision(50) << Q <<","<< E <<","<< static_cast<float>(E) <<","<< E_posit << endl;
    }

    return 0;
}