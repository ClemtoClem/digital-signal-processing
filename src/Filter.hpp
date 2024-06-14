#ifndef __FILTER_HPP
#define __FILTER_HPP

#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <iostream>

class Signal;

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
    virtual bool set(int order, int fc1, int fc2, int fs, FilterGabarit gabarit, AnalogFilter analogFilter, int rp = 0, int rs = 0) = 0;
    virtual void setup() = 0;
    virtual void reset() = 0;
    virtual Signal process(const Signal &input) = 0;
    virtual double eqdiff(double x) = 0;
    virtual void printCoefficients() = 0;
};

class IIRFilter: public BaseFilter {
public:
    IIRFilter();
    ~IIRFilter();

    // paramétrage du filtre
    bool set(int order, int fc1, int fc2, int fs, FilterGabarit gabarit, AnalogFilter analogFilter, int rp = 0, int rs = 0);
    

    void reset();

    // Méthode pour afficher les coefficients du filtre (à des fins de débogage)
    void printCoefficients();

    /* ----- Méthode 1 ----- */

    // calcul des coefficients
    void setup();
    
    Signal process(const Signal &input);

    /* ----- Méthode 2 ----- */

    // calcul des coefficients
    void setup2();

    // calcul de y(n) par application de l‘équation aux différences
    double eqdiff(double x);

private:
    void ButterworthCoefficients(); /* < Méthode 1 */

    bool mIsSetup;

    int mOrder;
    int mFc1, mFc2;
    int mFs; // Sample frequency
    int mRp; // passband ripple
    int mRs; // stopband attenuation
    FilterGabarit mGabarit;
    AnalogFilter mAnalogFilter;

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