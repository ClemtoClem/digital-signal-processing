#ifndef __FILTER_HPP
#define __FILTER_HPP

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <stdexcept>
#include "Signal.hpp"


class IIRFilter {
public:
    IIRFilter() : m_numerator(), m_denominator(), m_isSetup(false) {}

    void butterworthLowPass(int order, double fc, int fs);;

    void butterworthHighPass(int order, double fc, int fs);

    void butterworthBandPass(int order, double fc1, double fc2, int fs);

    void butterworthBandStop(int order, double fc1, double fc2, int fs);

    Signal filter(const Signal& input);

    // Méthode pour afficher les coefficients du filtre (à des fins de débogage)
    void printCoefficients();

private:

    // Méthode de convolution de deux vecteurs
    std::vector<double> convolve(const std::vector<double>& x, const std::vector<double>& y) {
        int m = x.size();
        int n = y.size();
        int p = m + n - 1;
        std::vector<double> result(p, 0.0);
        for (int i = 0; i < p; ++i) {
            for (int j = std::max(0, i - n + 1); j <= std::min(i, m - 1); ++j) {
                result[i] += x[j] * y[i - j];
            }
        }
        return result;
    }

    std::vector<double> m_numerator; // Coefficients du numérateur
    std::vector<double> m_denominator; // Coefficients du dénominateur
    bool m_isSetup; // Indicateur de configuration
};

#endif // __FILTER_HPP
