#ifndef __NOISE_HPP
#define __NOISE_HPP

#include <vector>
#include <random>
#include <cmath>
#include <complex>
#include "Signal.hpp"

enum class NoiseType {
    WHITE,
    PINK,
    BROWN
};

class Noise {
public:
    void set(NoiseType type, double samplingFreq, double amplitude);
    Signal process(const Signal &input);

private:
    NoiseType currentType;
    double fs; // Sampling frequency
    double amp; // Amplitude of the noise

    std::vector<double> whiteNoise(const Signal &input);
    std::vector<double> pinkNoise(const Signal &input);
    std::vector<double> brownNoise(const Signal &input);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> uniformDistribution;
    std::normal_distribution<double> normalDistribution;
    std::vector<double> pinkFilterCoefficients;
    std::vector<double> brownFilterState;
};

#endif // __NOISE_HPP