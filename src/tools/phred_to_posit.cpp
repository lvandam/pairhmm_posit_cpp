//
// Created by Laurens van Dam on 21/02/2018.
//

#include <iostream>
#include <posit>
#include <quadmath.h>

#define NBITS 32
#define ES 2

using namespace sw::unum;
using namespace std;

double score_to_probability(int Q) {
    return powf(10.f, -((double) Q) / 10.f);
}

int main() {
    int Q_min = 1;
    int Q_max = 100;

    cout << "Q,E_quad,E_float,E_posit" << endl;
    for(int Q = Q_min; Q <= Q_max; Q++) {
        __float128 E_quad = score_to_probability(Q);

        printf("%d,%.60f,%.60f,%.60f\n", Q, E_quad, float(E_quad), float(posit<NBITS,ES>(E_quad)));
    }

    return 0;
}