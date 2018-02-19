//
// Created by Laurens van Dam on 08/02/2018.
//

#ifndef PAIRHMM_SIMPLE_PAIRHMM_HPP
#define PAIRHMM_SIMPLE_PAIRHMM_HPP

#include <cmath>
#include "testcase.hpp"
#include "utils.hpp"

class PairHMMFloat {
private:
    static float score_to_probability(int i) {
        return powf(10.f, -((float) i) / 10.f);
    }

public:
    float compute_full_prob(Testcase *testcase) {
        int r, c;
        int ROWS = testcase->read_size;
        int COLS = testcase->haplotype_size;

        float M[350][350];
        float X[350][350];
        float Y[350][350];
        float p[350][6];
        float distm[350][350];

        int MM = 0, GapM = 1, MX = 2, XX = 3, MY = 4, YY = 5;

        p[0][MM] = 0;
        p[0][GapM] = 0;
        p[0][MX] = 0;
        p[0][XX] = 0;
        p[0][MY] = 0;
        p[0][YY] = 0;

        for (r = 1; r <= ROWS; r++) {
            int score_ins = testcase->ins_quals[r - 1] & 127;
            int score_del = testcase->del_quals[r - 1] & 127;
            int score_con = testcase->gcp_quals[r - 1] & 127;

            p[r][MM] = 1.0f - score_to_probability((score_ins + score_del) & 127);
            printDebug("p[%d][MM]=%.40f", r, p[r][MM]);

            p[r][GapM] = 0.9;
            printDebug("p[%d][GapM]=%.40f", r, p[r][GapM]);

            p[r][MX] = score_to_probability(score_ins);
            printDebug("p[%d][MX]=%.40f", r, p[r][MX]);

            p[r][XX] = 0.1;
            printDebug("p[%d][XX]=%.40f", r, p[r][XX]);

            p[r][MY] = score_to_probability(score_ins);
            printDebug("p[%d][MY]=%.40f", r, p[r][MY]);

            p[r][YY] = 0.1;
            printDebug("p[%d][YY]=%.40f", r, p[r][YY]);
        }

        // Initialize first row of every column
        for (c = 0; c <= COLS; c++) {
            M[0][c] = 0;
            X[0][c] = 0;
            Y[0][c] = ldexpf(1.f, 100) / (float) (testcase->haplotype_size);
            printDebug("Y[0][%d]=%.40f", c, Y[0][c]);
        }

        for (r = 1; r <= ROWS; r++) {
            M[r][0] = 0;
            X[r][0] = X[r - 1][0] * p[r - 1][XX];
            printDebug("X[%d][0]=%.40f", r, X[r][0]);
            Y[r][0] = 0;
        }

        for (r = 1; r <= ROWS; r++) {
            for (c = 1; c <= COLS; c++) {
                printDebug(">> [r,c] = %d,%d", r, c);

                char _rs = testcase->read_base[r - 1];
                printDebug("_rs=%.40f", _rs);

                char _hap = testcase->haplotype_base[c - 1];
                printDebug("_hap=%.40f", _hap);

                int score_base = testcase->base_quals[r - 1] & 127;
                printDebug("score_base=%.40f", score_base);

                distm[r][c] = score_to_probability(score_base);
                printDebug("distm[%d][%d]=%.40f", r, c, distm[r][c]);

                if (_rs == _hap || _rs == 'N' || _hap == 'N') {
                    distm[r][c] = 1 - distm[r][c];
                } else {
                    distm[r][c] = distm[r][c] / 3;
                }
                printDebug("distm_after[%d][%d]=%.40f", r, c, distm[r][c]);

                printDebug("");
            }
        }

        for (r = 1; r <= ROWS; r++) {
            for (c = 1; c <= COLS; c++) {
                M[r][c] = distm[r][c] * (M[r - 1][c - 1] * p[r][MM] + X[r - 1][c - 1] * p[r][GapM] + Y[r - 1][c - 1] * p[r][GapM]);
                printDebugValue(M[r][c], string("M[")+to_string(r)+"]["+to_string(c)+"]");

                X[r][c] = M[r - 1][c] * p[r][MX] + X[r - 1][c] * p[r][XX];
                printDebugValue(X[r][c], string("X[")+to_string(r)+"]["+to_string(c)+"]");

                Y[r][c] = M[r][c - 1] * p[r][MY] + Y[r][c - 1] * p[r][YY];
                printDebugValue(Y[r][c], string("Y[")+to_string(r)+"]["+to_string(c)+"]");
            }
        }

        printDebug("RESULT ACCUMULATION");
        float result = 0;
        for (c = 1; c <= COLS; c++) {
            result += M[ROWS][c] + X[ROWS][c];
            printDebugValue(result, string("result"));
        }

        return result;
    }
};


#endif //PAIRHMM_SIMPLE_PAIRHMM_HPP
