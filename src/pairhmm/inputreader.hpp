/**
    @author Laurens van Dam
    @date 08/02/2018
    @copyright 2018 All rights reserved.
**/

#ifndef PAIRHMM_SIMPLE_INPUTREADER_HPP
#define PAIRHMM_SIMPLE_INPUTREADER_HPP

#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include "testcase.hpp"
#include "utils.hpp"

using namespace std;

class InputReader {

public:
    std::vector<Testcase> from_file(const std::string& filename) {
        std::string line;
        std::ifstream input(filename);

        int read_size;
        input >> read_size;

        printDebug("Number of reads: %d", read_size);

        std::vector<Testcase> testcases;
        testcases.reserve(read_size);

        for(int i = 0; i < read_size; i++)
        {
            Testcase testcase;

            input >> testcase.read_size;
            input >> testcase.read_base;

            std::string read;
            std::getline(input, read);

            std::vector<std::string> tokens;
            boost::split(tokens, read, boost::is_any_of(" "));
            tokens.erase(tokens.begin());

            int token_idx = 0;

            try {
                // Base qualities
                for (int j = 0; j < testcase.read_size; j++)
                    testcase.base_quals[j] = (char)stoi(tokens[token_idx++]);

                // Insertion qualities
                for (int j = 0; j < testcase.read_size; j++)
                    testcase.ins_quals[j] = (char)stoi(tokens[token_idx++]);

                // Deletion qualities
                for (int j = 0; j < testcase.read_size; j++)
                    testcase.del_quals[j] = (char)stoi(tokens[token_idx++]);

                // Gap continuity qualities
                for (int j = 0; j < testcase.read_size; j++)
                    testcase.gcp_quals[j] = (char)stoi(tokens[token_idx++]);
            } catch(...) {
                throw std::runtime_error("Invalid input reads file");
            }

            input >> testcase.haplotype_size;
            input >> testcase.haplotype_base;

            testcases.push_back(testcase);
        }

        return testcases;
    }
};

#endif //PAIRHMM_SIMPLE_INPUTREADER_HPP
