#include "Noise.hpp"
#include "Signal.hpp"


WhiteNoise::WhiteNoise() : uniformDistribution(-1.0, 1.0) {}

void WhiteNoise::compute() {
    // Pas de coefficients spécifiques à calculer pour le bruit blanc
}

Signal WhiteNoise::process(const Signal &input) {
    Signal output = input;
    for (auto &sample : output) {
        sample += gain * uniformDistribution(generator);
    }
    return output;
}

double WhiteNoise::process(double input) {
    return input + gain * uniformDistribution(generator);
}




PinkNoise::PinkNoise() : gain(1.0), b(7, 0.0), a(7, 0.0), x_hist(7, 0.0), y_hist(7, 0.0) {}

void PinkNoise::compute() {
    // Coefficients based on Voss-McCartney algorithm for pink noise
    // You can adjust these coefficients as needed
    b = {0.99886, 0.99332, 0.96900, 0.86650, 0.55000, -0.7616, 0.115926};
    a = {1.0, -1.6915, 0.8928, -0.3092, 0.1990, -0.1606, 0.1065};
}

Signal PinkNoise::process(const Signal &input) {
    Signal output = input;
    for (size_t i = 0; i < input.size(); ++i) {
        output[i] = process(input[i]);
    }
    return output;
}

double PinkNoise::process(double input) {
    x_hist.insert(x_hist.begin(), input); // Add current input to history
    x_hist.pop_back(); // Remove oldest input

    double yn = 0.0;
    for (size_t i = 0; i < b.size(); ++i) {
        yn += b[i] * x_hist[i];
    }
    for (size_t i = 1; i < a.size(); ++i) {
        yn -= a[i] * y_hist[i - 1];
    }
    yn /= a[0];

    y_hist.insert(y_hist.begin(), yn); // Add current output to history
    y_hist.pop_back(); // Remove oldest output

    return yn * gain;
}




BrownNoise::BrownNoise() : previous(0.0), distribution(-1.0, 1.0) {}

void BrownNoise::compute() {
    // Pas de coefficients spécifiques à calculer pour le bruit marron
}

Signal BrownNoise::process(const Signal &input) {
    Signal output = input;
    for (auto &sample : output) {
        double white = distribution(generator);
        previous = (previous + (gain * white)) / 1.02;
        sample += previous * 0.1; // Ajuster le facteur si nécessaire
    }
    return output;
}

double BrownNoise::process(double input) {
    double white = distribution(generator);
    previous = (previous + (gain * white)) / 1.02;
    return input + previous * 0.1; // Ajuster le facteur si nécessaire
}
