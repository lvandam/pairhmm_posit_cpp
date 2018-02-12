//
// Created by Laurens van Dam on 08/02/2018.
//

#ifndef PAIRHMM_SIMPLE_PAIRHMM_POSIT_HPP
#define PAIRHMM_SIMPLE_PAIRHMM_POSIT_HPP

#include <cmath>
#include <posit>
#include "testcase.hpp"

using namespace std;
using namespace sw::unum;

class PairHMMPosit {
private:
    static constexpr size_t quire_capacity = 6;

    static posit<32,2> score_to_probability(int i) {
        return powf(10.f, -((float) i) / 10.f);
    }

public:
    double compute_full_prob(Testcase *testcase) {
        int r, c;
        int ROWS = testcase->read_size;
        int COLS = testcase->haplotype_size;

        std::vector<std::vector<posit<32, 2>>> M;
        std::vector<std::vector<posit<32, 2>>> X;
        std::vector<std::vector<posit<32, 2>>> Y;
        std::vector<std::vector<posit<32, 2>>> p;

        // Initialize matrices
        for (int i = 0; i < 350; i++) {
            std::vector<posit<32, 2>> row_m_x_y(350);
            M.push_back(row_m_x_y);
            X.push_back(row_m_x_y);
            Y.push_back(row_m_x_y);
            std::vector<posit<32, 2>> row_p(6);
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
            cout << "p["<<r<<"][MM] " << fixed << setprecision(30) << p[r][MM] << endl;

            p[r][GapM] = 0.9;
            cout << "p["<<r<<"][GapM] " << fixed << setprecision(30) << p[r][GapM] << endl;
            p[r][MX] = score_to_probability(score_ins);
            cout << "p["<<r<<"][MX] " << fixed << setprecision(30) << p[r][MX] << endl;
            p[r][XX] = 0.1;
            cout << "p["<<r<<"][XX] " << fixed << setprecision(30) << p[r][XX] << endl;
            p[r][MY] = score_to_probability(score_ins);
            cout << "p["<<r<<"][MY] " << fixed << setprecision(30) << p[r][MY] << endl;
            p[r][YY] = 0.1;
            cout << "p["<<r<<"][YY] " << fixed << setprecision(30) << p[r][YY] << endl;
        }

        // Initialize first row of every column
        for (c = 0; c <= COLS; c++) {
            M[0][c] = 0;
            X[0][c] = 0;
            Y[0][c] = ldexpf(1.f, 100) / (float) (testcase->haplotype_size);
            cout << "Y[0]["<<c<<"] " << fixed << setprecision(30) << Y[0][c] << endl;
        }

        // Initialize first column of every row
        for (r = 1; r <= ROWS; r++) {
            M[r][0] = 0;
            X[r][0] = X[r - 1][0] * p[r - 1][XX];
            cout << "X["<<r<<"][0] " << fixed << setprecision(30) << X[r][0] << endl;
            Y[r][0] = 0;
        }

        for (r = 1; r <= ROWS; r++) {
            for (c = 1; c <= COLS; c++) {

                char _rs = testcase->read_base[r - 1];
                char _hap = testcase->haplotype_base[c - 1];
                int score_base = testcase->base_quals[r - 1] & 127;

                posit<32,2> distm = score_to_probability(score_base);
                cout << "["<<r<<"]["<<c<<"] score_base="<<score_base<<" -- distm=" << fixed << setprecision(30) << distm << endl;

                if (_rs == _hap || _rs == 'N' || _hap == 'N') {
                    distm = 1 - distm;
                }
                else {
                    distm = distm / 3;
                }
                cout << "["<<r<<"]["<<c<<"] score_base="<<score_base<<" -- distm_after=" << fixed << setprecision(30) << distm << endl;

                M[r][c] = distm * (M[r - 1][c - 1] * p[r][MM] + X[r - 1][c - 1] * p[r][GapM] + Y[r - 1][c - 1] * p[r][GapM]);
                cout << "M["<<r<<"]["<<c<<"]=" << fixed << setprecision(30) << M[r][c] << endl;
                X[r][c] = M[r - 1][c] * p[r][MX] + X[r - 1][c] * p[r][XX];
                cout << "X["<<r<<"]["<<c<<"]=" << fixed << setprecision(30) << X[r][c] << endl;
                Y[r][c] = M[r][c - 1] * p[r][MY] + Y[r][c - 1] * p[r][YY];
                cout << "Y["<<r<<"]["<<c<<"]=" << fixed << setprecision(30) << Y[r][c] << endl;
            }
        }

        // Accumulator
        quire<32, 2, quire_capacity> result;
        for (c = 1; c <= COLS; c++) {
            result += (M[ROWS][c] + X[ROWS][c]).convert_to_scientific_notation();
        }

        // Convert back to 64-bit double
        return double(result.to_value());
    }
};

#endif //PAIRHMM_SIMPLE_PAIRHMM_POSIT_HPP
