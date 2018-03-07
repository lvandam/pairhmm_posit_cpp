/**
    @author Laurens van Dam
    @date 01/03/2018
    @copyright 2018 All rights reserved.
**/

#ifndef PAIRHMM_PAIRHMM_POSIT_HPP
#define PAIRHMM_PAIRHMM_POSIT_HPP

#include <cmath>
#include <posit/posit>
#include "utils.hpp"
#include "testcase.hpp"
#include "debug_values.hpp"

using namespace std;
using namespace sw::unum;

void printDebug(const char* format, ...);

template<class T, class QUIRE>
class PairHMMPosit {
private:
    static T score_to_probability(int Q) {
        return powl(10.f, -((long double) Q) / 10.f);
    }

public:
    PairHMMPosit() : INITIAL_CONSTANT(ldexpf(1.f, 100)) {
    }

    explicit PairHMMPosit(long double initial_constant) : INITIAL_CONSTANT(initial_constant) {
    }

    long double INITIAL_CONSTANT;

    DebugValues<T> debug_values;

    float compute_full_prob(Testcase *testcase) {
        int r, c;
        int ROWS = testcase->read_size;
        int COLS = testcase->haplotype_size;

        std::vector<std::vector<T>> M;
        std::vector<std::vector<T>> X;
        std::vector<std::vector<T>> Y;
        std::vector<std::vector<T>> p;
        std::vector<std::vector<T>> distm;

        // Initialize matrices
        for (int i = 0; i < 350; i++) {
            std::vector<T> row_m_x_y(350);
            M.push_back(row_m_x_y);
            X.push_back(row_m_x_y);
            Y.push_back(row_m_x_y);
            distm.push_back(row_m_x_y);
            std::vector<T> row_p(6);
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
            debug_values.debugValue(p[r][MM], "p[%d][MM]", r);

            p[r][GapM] = 0.9;
            debug_values.debugValue(p[r][GapM], "p[%d][GapM]", r);

            p[r][MX] = score_to_probability(score_ins);
            debug_values.debugValue(p[r][MX], "p[%d][MX]", r);

            p[r][XX] = 0.1;
            debug_values.debugValue(p[r][XX], "p[%d][XX]", r);

            p[r][MY] = score_to_probability(score_ins);
            debug_values.debugValue(p[r][MY], "p[%d][MY]", r);

            p[r][YY] = 0.1;
            debug_values.debugValue(p[r][YY], "p[%d][YY]", r);
        }

        // Initialize first row of every column
        for (c = 0; c <= COLS; c++) {
            M[0][c] = 0;
            X[0][c] = 0;
            Y[0][c] = INITIAL_CONSTANT / (float) (testcase->haplotype_size);
            debug_values.debugValue(Y[0][c], "Y[0][%d]", c);
        }

        // Initialize first column of every row
        for (r = 1; r <= ROWS; r++) {
            M[r][0] = 0;
            X[r][0] = X[r - 1][0] * p[r - 1][XX];
            debug_values.debugValue(X[r][0], "X[%d][0]", r);
            Y[r][0] = 0;
        }

        for (r = 1; r <= ROWS; r++) {
            for (c = 1; c <= COLS; c++) {
                printDebug(">> [r,c] = %d,%d", r, c);

                char _rs = testcase->read_base[r - 1];
                debug_values.debugValue(_rs, "_rs");

                char _hap = testcase->haplotype_base[c - 1];
                debug_values.debugValue(_hap, "_hap");

                int score_base = testcase->base_quals[r - 1] & 127;
                debug_values.debugValue(score_base, "score_base");

                distm[r][c] = score_to_probability(score_base);
                debug_values.debugValue(distm[r][c], "distm[%d][%d]", r, c);

                if (_rs == _hap || _rs == 'N' || _hap == 'N') {
                    distm[r][c] = 1 - distm[r][c];
                } else {
                    distm[r][c] = distm[r][c] / 3;
                }
                debug_values.debugValue(distm[r][c], "distm_after[%d][%d]", r, c);

                printDebug("");
            }
        }

        for (r = 1; r <= ROWS; r++) {
            for (c = 1; c <= COLS; c++) {
                QUIRE Mq(0), Xq(0), Yq(0);
                posit<32,2> Mposit, Xposit, Yposit;
                Mq += quire_mul(distm[r][c] * M[r - 1][c - 1], p[r][MM]);
                Mq += quire_mul(distm[r][c] * X[r - 1][c - 1], p[r][GapM]);
                Mq += quire_mul(distm[r][c] * Y[r - 1][c - 1], p[r][GapM]);
                Mposit.convert(Mq.to_value());
                M[r][c] = Mposit;

                debug_values.debugValue(M[r][c], "M[%d][%d]", r, c);

                Xq += quire_mul(M[r - 1][c], p[r][MX]);
                Xq += quire_mul(X[r - 1][c], p[r][XX]);
                Xposit.convert(Xq.to_value());
                X[r][c] = Xposit;

                debug_values.debugValue(X[r][c], "X[%d][%d]", r, c);

                Yq += quire_mul(M[r][c - 1], p[r][MY]);
                Yq += quire_mul(Y[r][c - 1], p[r][YY]);
                Yposit.convert(Yq.to_value());
                Y[r][c] = Yposit;

                debug_values.debugValue(Y[r][c], "Y[%d][%d]", r, c);
            }
        }

        printDebug("RESULT ACCUMULATION");
        QUIRE quire(0);
        for (c = 1; c <= COLS; c++) {
            constexpr int bits = std::numeric_limits<long double>::digits - 1;
            quire += M[ROWS][c].convert_to_scientific_notation();
            quire += X[ROWS][c].convert_to_scientific_notation();
            debug_values.debugValue(float(quire.to_value()), "result[%d]", c);
        }

        // Convert back to float
        return float(quire.to_value());
    }
};

#endif //PAIRHMM_PAIRHMM_POSIT_HPP
