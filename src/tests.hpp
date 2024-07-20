#ifndef __TEST_HPP
#define __TEST_HPP

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <numeric>
#include <complex>
#include <string>
#include <cstring>
#include <chrono>
#include <algorithm>
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <time.h>

#include "Signal.hpp"
#include "Spectrum.hpp"
#include "Filter.hpp"
#include "Demodulator.hpp"
#include "PID.hpp"
#include "CSVFile.hpp"
#include "Noise.hpp"
#include "Window.hpp"

/**
 * @brief Test algorithm
 * @param[in] args Arguments
 * @note Write help message if the argument "help" is provided
 */
int test_simple(const std::vector<std::string> &args);

/**
 * @brief Test the spectrum algorithm
 * @param[in] args Arguments
 * @note Write help message if the argument "help" is provided
 */
int test_spectrum(const std::vector<std::string> &args);

/**
 * @brief Test the demodulation algorithm
 * @param[in] args Arguments
 * @note Write help message if the argument "help" is provided
 */
int test_demodulation(const std::vector<std::string> &args);


/**
 * @brief Test the demodulation algorithm n°2
 * @param[in] args Arguments
 * @note Write help message if the argument "help" is provided
 */
int test_demodulation2(const std::vector<std::string> &args);


#endif // __TEST_HPP