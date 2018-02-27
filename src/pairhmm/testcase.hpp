/**
    @author Laurens van Dam
    @date 08/02/2018
    @copyright 2018 All rights reserved.
**/

#ifndef PAIRHMM_SIMPLE_TESTCASE_HPP
#define PAIRHMM_SIMPLE_TESTCASE_HPP

struct Testcase {
    int read_size;
    char read_base[501];
    char base_quals[501];
    char ins_quals[501];
    char del_quals[501];
    char gcp_quals[501];
    int haplotype_size;
    char haplotype_base[501];
};

#endif //PAIRHMM_SIMPLE_TESTCASE_HPP
