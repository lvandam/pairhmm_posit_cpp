/**
    @author Laurens van Dam
    @date 27/02/2018
    @copyright 2018 All rights reserved.
**/

#ifndef PAIRHMM_PAIRHMM_HPP
#define PAIRHMM_PAIRHMM_HPP

#include <cmath>
#include <posit/posit>
#include "utils.hpp"
#include "testcase.hpp"
#include "debug_values.hpp"

using namespace std;
using namespace sw::unum;

template<class T, class QUIRE>
class PairHMM {
private:
    static T score_to_probability(int Q) {
        return powl(10.f, -((long double) Q) / 10.f);
    }

public:
    PairHMM() : INITIAL_CONSTANT(ldexpf(1.f, 100)) {
    }

    explicit PairHMM(long double initial_constant) : INITIAL_CONSTANT(initial_constant) {
    }

    T INITIAL_CONSTANT;

    DebugValues<cpp_dec_float_50> debug_values;

    cpp_dec_float_50 compute_full_prob(Testcase *testcase) {
        int r, c;
        int ROWS = testcase->read_size;
        int COLS = testcase->haplotype_size;

        std::vector<std::vector<T>> M;
        std::vector<std::vector<T>> X;
        std::vector<std::vector<T>> Y;
        std::vector<std::vector<T>> p;
        std::vector<std::vector<T>> distm;

        // Initialize matrices
        std::vector<T> row_m_x_y(350), row_p(6);
        for (int i = 0; i < 350; i++) {
            M.push_back(row_m_x_y);
            X.push_back(row_m_x_y);
            Y.push_back(row_m_x_y);
            distm.push_back(row_m_x_y);
            p.push_back(row_p);
        }

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
            p[r][GapM] = 0.9;
            p[r][MX] = score_to_probability(score_ins);
            p[r][XX] = 0.1;
            p[r][MY] = score_to_probability(score_ins);
            p[r][YY] = 0.1;

            debug_values.debugValue(p[r][MM], "p[%d][MM]", r);
            debug_values.debugValue(p[r][GapM], "p[%d][GapM]", r);
            debug_values.debugValue(p[r][MX], "p[%d][MX]", r);
            debug_values.debugValue(p[r][XX], "p[%d][XX]", r);
            debug_values.debugValue(p[r][MY], "p[%d][MY]", r);
            debug_values.debugValue(p[r][YY], "p[%d][YY]", r);
        }

        // Initialize first row of every column
        for (c = 0; c <= COLS; c++) {
            M[0][c] = 0;
            X[0][c] = 0;
            Y[0][c] = static_cast<T>(INITIAL_CONSTANT) / testcase->haplotype_size;

            debug_values.debugValue(Y[0][c], "Y[0][%d]", c);
        }

        // Initialize first column of every row
        for (r = 1; r <= ROWS; r++) {
            M[r][0] = 0;
            X[r][0] = X[r - 1][0] * p[r - 1][XX];
            Y[r][0] = 0;

            debug_values.debugValue(X[r][0], "X[%d][0]", r);
        }

        for (r = 1; r <= ROWS; r++) {
            for (c = 1; c <= COLS; c++) {
                char _rs = testcase->read_base[r - 1];
                char _hap = testcase->haplotype_base[c - 1];
                int score_base = testcase->base_quals[r - 1] & 127;

                distm[r][c] = score_to_probability(score_base);

                if (_rs == _hap || _rs == 'N' || _hap == 'N') {
                    distm[r][c] = 1 - distm[r][c];
                } else {
                    distm[r][c] = distm[r][c] / 3;
                }

                debug_values.printDebug(">> [r,c] = %d,%d", r, c);
                debug_values.debugValue(_rs, "_rs");
                debug_values.debugValue(_hap, "_hap");
                debug_values.debugValue(score_base, "score_base");
                debug_values.debugValue(distm[r][c], "distm[%d][%d]", r, c);
                debug_values.debugValue(distm[r][c], "distm_after[%d][%d]", r, c);
                debug_values.printDebug("");
            }
        }

        for (r = 1; r <= ROWS; r++) {
            for (c = 1; c <= COLS; c++) {
                M[r][c] = distm[r][c] * (M[r - 1][c - 1] * p[r][MM] + X[r - 1][c - 1] * p[r][GapM] + Y[r - 1][c - 1] * p[r][GapM]);
                X[r][c] = M[r - 1][c] * p[r][MX] + X[r - 1][c] * p[r][XX];
                Y[r][c] = M[r][c - 1] * p[r][MY] + Y[r][c - 1] * p[r][YY];

                debug_values.debugValue(M[r][c], "M[%d][%d]", r, c);
                debug_values.debugValue(X[r][c], "X[%d][%d]", r, c);
                debug_values.debugValue(Y[r][c], "Y[%d][%d]", r, c);
            }
        }

        // Result accumulation
        T result = 0;
        for(c = 1; c <= COLS; c++) {
            result += M[ROWS][c] + X[ROWS][c];
            debug_values.debugValue((cpp_dec_float_50)result, "result[%d]", c);
        }

        return cpp_dec_float_50(result);
    }
};

#endif //PAIRHMM_PAIRHMM_HPP
