#ifndef __UTILS_HPP
#define __UTILS_HPP

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

bool stringToBool(const std::string &str);

std::string waveformTypeToString(WaveformType waveform);

WaveformType stringToWaveformType(const std::string &str);

// for string delimiter
std::vector<std::string> split(std::string s, std::string delimiter);

int calculateDecimation(double f, int N_per_period);

uint32_t getTimeDelay(int decimation);

int convertToInteger(const std::string& str);

// Fonction pour calculer le PGCD de deux nombres
int gcd(int a, int b);

// Fonction pour calculer le PGCD de plusieurs nombres
int gcd_multiple(const std::vector<int>& numbers);

// Fonction pour générer une période complète du signal
std::vector<double> generate_waveform_with_n_sinus(size_t buffer_size, std::vector<int> frequencies, std::vector<double> amplitudes, int &fundamental_frequency);

#endif // __UTILS_HPP