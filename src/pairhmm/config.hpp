/**
    @author Laurens van Dam
    @date 19/02/2018
    @copyright 2018 All rights reserved.
**/

#ifndef PAIRHMM_SIMPLE_CONFIG_HPP
#define PAIRHMM_SIMPLE_CONFIG_HPP

#ifndef DEBUG_PRECISION
#define DEBUG_PRECISION 40
#endif // DEBUG_PRECISION

#ifndef INITIAL_CONSTANT
#define INITIAL_CONSTANT (ldexpf(1.f, 100)) // 2^100
//#define INITIAL_CONSTANT (ldexpf(1.f, 10)) // 2^10
#endif // INITIAL_CONSTANT

#endif //PAIRHMM_SIMPLE_CONFIG_HPP
