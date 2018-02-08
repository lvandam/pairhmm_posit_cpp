//
// Created by Laurens van Dam on 08/02/2018.
//

#ifndef PAIRHMM_SIMPLE_PAIRHMM_HPP
#define PAIRHMM_SIMPLE_PAIRHMM_HPP

#include <cmath>
#include "testcase.hpp"

class PairHMMFloat {
private:
    static float ph2pr(int i) {
        return powf(10.f, -((float) i) / 10.f);
    }

public:
    double compute_full_prob(Testcase *testcase) {
        int r, c;
        int ROWS = testcase->read_size;
        int COLS = testcase->haplotype_size;

        double M[350][350];
        double X[350][350];
        double Y[350][350];
        double p[350][6];

        int MM = 0, GapM = 1, MX = 2, XX = 3, MY = 4, YY = 5;

        p[0][MM] = 0;
        p[0][GapM] = 0;
        p[0][MX] = 0;
        p[0][XX] = 0;
        p[0][MY] = 0;
        p[0][YY] = 0;

        for (r = 1; r <= ROWS; r++) {
            int _i = testcase->ins_quals[r - 1] & 127;
            int _d = testcase->del_quals[r - 1] & 127;
            int _c = testcase->gcp_quals[r - 1] & 127;

            p[r][MM] = 1.0f - ph2pr((_i + _d) & 127);

            p[r][GapM] = 0.9;
            p[r][MX] = ph2pr(_i);
            p[r][XX] = 0.1;
            p[r][MY] = ph2pr(_i);
            p[r][YY] = 0.1;
        }

        for (c = 0; c <= COLS; c++) {
            M[0][c] = 0;
            X[0][c] = 0;
            Y[0][c] = ldexpf(1.f, 100) / (float) (testcase->haplotype_size);
        }

        for (r = 1; r <= ROWS; r++) {
            M[r][0] = 0;
            X[r][0] = X[r - 1][0] * p[r - 1][XX];
            Y[r][0] = 0;
        }

        for (r = 1; r <= ROWS; r++)
            for (c = 1; c <= COLS; c++) {

                char _rs = testcase->read_base[r - 1];
                char _hap = testcase->haplotype_base[c - 1];

                int _q = testcase->base_quals[r - 1] & 127;
                float distm = ph2pr(_q);

                if (_rs == _hap || _rs == 'N' || _hap == 'N')
                    distm = 1 - distm;
                else
                    distm = distm / 3;

                M[r][c] = distm * (M[r - 1][c - 1] * p[r][MM] + X[r - 1][c - 1] * p[r][GapM] + Y[r - 1][c - 1] * p[r][GapM]);
                X[r][c] = M[r - 1][c] * p[r][MX] + X[r - 1][c] * p[r][XX];
                Y[r][c] = M[r][c - 1] * p[r][MY] + Y[r][c - 1] * p[r][YY];
            }

        double result = 0;
        for (c = 1; c <= COLS; c++) {
            result += M[ROWS][c] + X[ROWS][c];
        }

        return result;
    }
};


#endif //PAIRHMM_SIMPLE_PAIRHMM_HPP
