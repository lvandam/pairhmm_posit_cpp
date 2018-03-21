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
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;
using namespace sw::unum;
using boost::multiprecision::cpp_dec_float_50;

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

    long double INITIAL_CONSTANT;

    DebugValues<T> debug_values;

    cpp_dec_float_50 compute_full_prob(Testcase *testcase) {
        int r, c;
        int ROWS = testcase->read_size;
        int COLS = testcase->haplotype_size;

        std::vector<std::vector<T>> M, X, Y, p, distm;

        // Initialize matrices
        for (int i = 0; i < 350; i++) {
            std::vector<T> row_m_x_y(350), row_p(6);

            M.push_back(row_m_x_y);
            X.push_back(row_m_x_y);
            Y.push_back(row_m_x_y);
            distm.push_back(row_m_x_y);
            p.push_back(row_p);
        }

        const int MM = 0, GapM = 1, MX = 2, XX = 3, MY = 4, YY = 5;

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
            }
        }

        constexpr int bits = std::numeric_limits<T>::digits - 1;

        QUIRE result_quire(0);
        for (r = 1; r <= ROWS; r++) {
            for (c = 1; c <= COLS; c++) {
                QUIRE Mq(0), Xq(0), Yq(0);

                const int bits_mult = bits * 2 + 2;
                const int bits_mult_mult = bits_mult * 2 + 2;

                // Calculation of M[r][c]
                value<bits_mult> distm_M, distm_X, distm_Y, pMM, pGapM;
                value<bits_mult_mult> resM, resX, resY;

                module_multiply<bits, bits_mult>(value<bits>(distm[r][c]), value<bits>(M[r - 1][c - 1]), distm_M);
                module_multiply<bits, bits_mult>(value<bits>(distm[r][c]), value<bits>(X[r - 1][c - 1]), distm_X);
                module_multiply<bits, bits_mult>(value<bits>(distm[r][c]), value<bits>(Y[r - 1][c - 1]), distm_Y);

                module_multiply<bits_mult, bits_mult_mult>(distm_M, value<bits_mult>(p[r][MM]), resM);
                module_multiply<bits_mult, bits_mult_mult>(distm_X, value<bits_mult>(p[r][GapM]), resX);
                module_multiply<bits_mult, bits_mult_mult>(distm_Y, value<bits_mult>(p[r][GapM]), resY);

                Mq += resM;
                Mq += resX;
                Mq += resY;

                M[r][c] = static_cast<T>(Mq.to_value());

                // Calculation of X[r][c]
                value<bits_mult> M_p_MX, X_p_XX;
                module_multiply<bits, bits_mult>(value<bits>(M[r - 1][c]), value<bits>(p[r][MX]), M_p_MX);
                module_multiply<bits, bits_mult>(value<bits>(X[r - 1][c]), value<bits>(p[r][XX]), X_p_XX);

                Xq += M_p_MX;
                Xq += X_p_XX;

                X[r][c] = static_cast<T>(Xq.to_value());

                // Calculation of Y[r][c]
                value<bits_mult> M_p_MY, Y_p_YY;
                module_multiply<bits, bits_mult>(value<bits>(M[r][c - 1]), value<bits>(p[r][MY]), M_p_MY);
                module_multiply<bits, bits_mult>(value<bits>(Y[r][c - 1]), value<bits>(p[r][YY]), Y_p_YY);

                Yq += M_p_MY;
                Yq += Y_p_YY;

                Y[r][c] = static_cast<T>(Yq.to_value());

                debug_values.debugValue(M[r][c], "M[%d][%d]", r, c);
                debug_values.debugValue(X[r][c], "X[%d][%d]", r, c);
                debug_values.debugValue(Y[r][c], "Y[%d][%d]", r, c);

                // Result accumulation
                if(r == ROWS) {
                    result_quire += Mq.to_value();
                    result_quire += Xq.to_value();

                    debug_values.debugValue((T)result_quire.to_value(), "result[%d]", c);
                }
            }
        }

        // Convert back to decimal50
        return cpp_dec_float_50(result_quire.to_value());
    }
};

#endif //PAIRHMM_PAIRHMM_HPP
