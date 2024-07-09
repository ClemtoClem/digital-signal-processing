#ifndef __FILTER_HPP
#define __FILTER_HPP

#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include "Signal.hpp"
#include "Spectrum.hpp"

enum class FilterGabarit {
    LOW_PASS,
    HIGH_PASS,
    BAND_PASS,
    BAND_STOP,
};

enum class AnalogFilter {
    BESSEL,      // IIR Bessel /Thomsonfilter design
    BUTTERWORTH, // IIR Butterworth filter design
    CHEBYSHEV1,  // IIR Chebyshev Type I filter design
    CHEBYSHEV2,  // IIR Chebyshev Type II filter design
    ELLIPTIC,    // IIR Cauer/Elliptic filter design
};

// classe abstraire d'un filtre
class BaseFilter {
public:
    virtual bool set(int order, double fc1, double fc2, FilterGabarit gabarit, AnalogFilter analogFilter, double rp = 0, double rs = 0) = 0;
    virtual bool isSetup() const = 0;
    virtual void setup() = 0;
    virtual void reset() = 0;
    virtual Signal apply(const Signal &input) = 0;
    virtual double apply(double x) = 0;
    virtual void printCoefficients() = 0;
};

class IIRFilter: public BaseFilter {
public:
    IIRFilter();
    ~IIRFilter();

    // paramétrage du filtre
    bool set(int order, double fc1, double fc2, FilterGabarit gabarit, AnalogFilter analogFilter, double rp = 0, double rs = 0);
    
    /**
     * @brief Check if the filter is setup
     * @return True if the filter is setup
     */
    bool isSetup() const { return _isSetup; }

    void reset();

    // Méthode pour afficher les coefficients du filtre (à des fins de débogage)
    void printCoefficients();

    /* ----- Méthode 1 ----- */

    // calcul des coefficients
    void setup();
    
    Signal apply(const Signal &input);

    // calcul de y(n) par application de l‘équation aux différences
    double apply(double y);

    Spectrum frequency_response(size_t num_points);

    /* ----- Méthode 2 ----- */

    // calcul des coefficients
    void setup2();

private:
    void ButterworthCoefficients(); /* < Méthode 1 */

    bool _isSetup;

    int _order;
    double _fc1, _fc2;
    double _rp; // passband ripple
    double _rs; // stopband attenuation
    FilterGabarit _gabarit;
    AnalogFilter _analogFilter;

    std::vector<double> _a, _b;
    std::vector<std::complex<double>> _z, _p;
    double _k; // gain
    std::vector<double> _xMem, _yMem;
};


/* ------------------------------------------------------------------------ */
/* Méthode 2 */

// Return (z,p,k) for analog prototype of Nth-order Butterworth filter
void ButterworthAnalogPrototype(int N, std::vector<std::complex<double>> &z, std::vector<std::complex<double>> &p, double &k);

void TransformLowpassToLowpass(std::vector<std::complex<double>> &z, std::vector<std::complex<double>> &p, double &k, double warped);

void TransformLowpassToHighpass(std::vector<std::complex<double>> &z, std::vector<std::complex<double>> &p, double &k, double warped);

void TransformLowpassToBandpass(std::vector<std::complex<double>> &z, std::vector<std::complex<double>> &p, double &k, double warped1, double warped2);

void TransformLowpassToBandStop(std::vector<std::complex<double>> &z, std::vector<std::complex<double>> &p, double &k, double warped1, double warped2);

void IIRBilinearTransformation(std::vector<std::complex<double>> &z, std::vector<std::complex<double>> &p, double &k, double fs);

void PolynomialTransfer(const std::vector<std::complex<double>> &z, const std::vector<std::complex<double>> &p, double k, std::vector<double> &a, std::vector<double> &b);

#endif // __FILTER_HPP