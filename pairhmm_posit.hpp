//
// Created by Laurens van Dam on 08/02/2018.
//

#ifndef PAIRHMM_SIMPLE_PAIRHMM_POSIT_HPP
#define PAIRHMM_SIMPLE_PAIRHMM_POSIT_HPP

#include <cmath>
#include <posit>
#include "testcase.hpp"

#define NBITS 32
#define ES 2

using namespace std;
using namespace sw::unum;

class PairHMMPosit {
private:
    static constexpr size_t quire_capacity = 6;

    static posit<NBITS, ES> score_to_probability(int i) {
        cout << "score to probability " << i << endl;
        return powf(10.f, -((float) i) / 10.f);
    }

public:
    float compute_full_prob(Testcase *testcase) {
        int r, c;
        int ROWS = testcase->read_size;
        int COLS = testcase->haplotype_size;

        std::vector<std::vector<posit<NBITS, ES>>> M;
        std::vector<std::vector<posit<NBITS, ES>>> X;
        std::vector<std::vector<posit<NBITS, ES>>> Y;
        std::vector<std::vector<posit<NBITS, ES>>> p;
        std::vector<std::vector<posit<NBITS, ES>>> distm;

        // Initialize matrices
        for (int i = 0; i < 350; i++) {
            std::vector<posit<NBITS, ES>> row_m_x_y(350);
            M.push_back(row_m_x_y);
            X.push_back(row_m_x_y);
            Y.push_back(row_m_x_y);
            distm.push_back(row_m_x_y);
            std::vector<posit<NBITS, ES>> row_p(6);
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
        }

        // Initialize first row of every column
        for (c = 0; c <= COLS; c++) {
            M[0][c] = 0;
            X[0][c] = 0;
            Y[0][c] = ldexpf(1.f, 100) / (float) (testcase->haplotype_size);
        }

        // Initialize first column of every row
        for (r = 1; r <= ROWS; r++) {
            M[r][0] = 0;
            X[r][0] = X[r - 1][0] * p[r - 1][XX];
            Y[r][0] = 0;
        }

        for (r = 1; r <= ROWS; r++) {
            for (c = 1; c <= COLS; c++) {

                char _rs = testcase->read_base[r - 1];
                char _hap = testcase->haplotype_base[c - 1];
                int score_base = testcase->base_quals[r - 1] & 127;

                distm[r][c] = score_to_probability(score_base);

                if (_rs == _hap || _rs == 'N' || _hap == 'N') {
                    distm[r][c] = 1 - distm[r][c];
                }
                else {
                    distm[r][c] = distm[r][c] / 3;
                }
            }
        }

        // TODO vectorize this
        /* The vector to keep the anti-diagonals of M, X, Y*/
        /* Current: M_t, X_t, Y_t */
        /* Previous: M_t_1, X_t_1, Y_t_1 */
        /* Previous to previous: M_t_2, X_t_2, Y_t_2 */
        for (r = 1; r <= ROWS; r++) {
            for (c = 1; c <= COLS; c++) {
                M[r][c] = distm[r][c] * (M[r - 1][c - 1]*p[r][MM] + X[r - 1][c - 1]*p[r][GapM] + Y[r - 1][c - 1]*p[r][GapM]);
                X[r][c] = M[r - 1][c]*p[r][MX] + X[r - 1][c]*p[r][XX];
                Y[r][c] = M[r][c - 1]*p[r][MY] + Y[r][c - 1]*p[r][YY];
            }
        }

        // Accumulator
        quire<NBITS, ES, quire_capacity> result;
        for (c = 1; c <= COLS; c++) {
            result += (M[ROWS][c] + X[ROWS][c]).convert_to_scientific_notation();
        }

        // Convert back to float
        return float(result.to_value());
    }
};

#endif //PAIRHMM_SIMPLE_PAIRHMM_POSIT_HPP
