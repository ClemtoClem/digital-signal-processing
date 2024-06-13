#ifndef __MATH_HPP
#define __MATH_HPP

#include <complex>
#include <cmath>
#include <vector>

// calul les racines d'un racine n ième
inline std::vector<std::complex<double>> sqrtn(const std::complex<double> &value, int N) {
    std::vector<std::complex<double>> roots;
    if (N <= 0) {
        throw invalid_argument("N must be greater than 0");
    }

    double r = std::abs(value); // Module of the complex number
    double theta = std::arg(value); // Argument of the complex number

    double root_r = pow(r, 1.0 / N); // Nth root of the module

    for (int k = 0; k < N; ++k) {
        double angle = (theta + 2 * M_PI * k) / N;
        std::complex<double> root = std::polar(root_r, angle);
        roots.push_back(root);
    }

    return roots;
}

// calul de la racine n ième inverse
inline std::complex<double> sqrtin(const std::complex<double> &value, int N) {
    if (N == 0) {
        throw invalid_argument("N must be non-zero");
    }

    double r = std::abs(value); // Module of the complex number
    double theta = std::arg(value); // Argument of the complex number

    // Calculate the nth root of the module
    double root_r = pow(r, 1.0 / N);

    // Calculate the argument for the principal root
    double root_theta = theta / N;

    // Return the principal root in polar form
    return std::polar(root_r, root_theta);
}

#endif // __MATH_HPP