#include "Noise.hpp"

void Noise::setParams(NoiseType type, double samplingFreq, double amplitude) {
    currentType = type;
    fs = samplingFreq;
    amp = amplitude;

    // Initialize random number generators
    generator = std::default_random_engine(std::random_device{}());
    uniformDistribution = std::uniform_real_distribution<double>(-1.0, 1.0);
    normalDistribution = std::normal_distribution<double>(0.0, 1.0);

    // Initialize filter coefficients and state if needed
    if (type == NoiseType::PINK) {
        pinkFilterCoefficients = {0.02109238, -0.02613793, 0.03399677, -0.04781565, 0.07111543, -0.11625247, 0.22685173, -0.62218065};
    } else if (type == NoiseType::BROWN) {
        brownFilterState = {0.0};
    }
}

Signal Noise::process(const Signal &input) {
    Signal result(input.size(), input.getSamplingFrequency());

    // Generate noise samples based on the selected noise type
    std::vector<double> noiseSamples;
    switch (currentType) {
        case NoiseType::WHITE:
            noiseSamples = whiteNoise(input);
            break;
        case NoiseType::PINK:
            noiseSamples = pinkNoise(input);
            break;
        case NoiseType::BROWN:
            noiseSamples = brownNoise(input);
            break;
    }

    // Apply the noise to the input signal
    for (size_t i = 0; i < input.size(); ++i) {
        result[i] = input[i] + amp * noiseSamples[i];
    }

    return result;
}

std::vector<double> Noise::whiteNoise(const Signal &input) {
    size_t numSamples = input.size();
    std::vector<double> noise(numSamples);

    for (size_t i = 0; i < numSamples; ++i) {
        noise[i] = uniformDistribution(generator);
    }

    return noise;
}

std::vector<double> Noise::pinkNoise(const Signal &input) {
    size_t numSamples = input.size();
    std::vector<double> noise(numSamples);

    double pinkedSample = 0.0;
    for (size_t i = 0; i < numSamples; ++i) {
        pinkedSample += normalDistribution(generator);

        // Apply the pink noise filter coefficients
        for (size_t j = 0; j < pinkFilterCoefficients.size(); ++j) {
            if (i >= j) {
                pinkedSample -= pinkFilterCoefficients[j] * noise[i - j];
            }
        }

        noise[i] = pinkedSample;
    }

    return noise;
}

std::vector<double> Noise::brownNoise(const Signal &input) {
    size_t numSamples = input.size();
    std::vector<double> noise(numSamples);

    double brownerSample = brownFilterState[0];
    for (size_t i = 0; i < numSamples; ++i) {
        brownerSample += normalDistribution(generator);

        // Update the brown noise filter state
        brownFilterState[0] = brownerSample;

        noise[i] = brownerSample;
    }

    return noise;
}