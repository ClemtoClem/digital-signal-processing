#include "Filter.hpp"
#include <algorithm>
#include <numeric>
#include "Signal.hpp"

IIRFilter::IIRFilter() : mIsSetup(false) {}

IIRFilter::~IIRFilter() {
    // Destructor logic if needed (automatic cleanup of vectors)
}

bool IIRFilter::set(int order, int fc1, int fc2, int fs, FilterGabarit gabarit, AnalogFilter analogFilter, int rp, int rs)
{
    if (order < 1) {
        std::cerr << "Warning: L'ordre du filtre doit être supérieur ou égal à 1.\n";
        return false;
    }
    if (fc1 < 0 || fc2 < 0) {
        std::cerr << "Warning: Les fréquences de coupures doivent être positives.\n";
        return false;
    }
    if (fc1 >= fc2 && (gabarit == FilterGabarit::BAND_PASS || gabarit == FilterGabarit::BAND_STOP)) {
        std::cerr << "Warning: la fréquence de coupure fc1 doit être inférieur à la fréquence de coupre fc2.\n";
        return false;
    }

    mOrder = order;
    mFc1 = fc1;
    mFc2 = fc2;
    mFs = fs;
    mGabarit = gabarit;
    mAnalogFilter = analogFilter;
    mRp = rp; // passband ripple      (Chebyshev I, Chebyshev II & elliptic filter)
    mRs = rs; // stopband attenuation (Chebyshev II & elliptic filter)

    mIsSetup = false;
    return true;
}

void IIRFilter::reset() {
    // Reset memory x[n] and y[n] buffers
    _xMem.resize(mOrder+1, 0.0);
    _yMem.resize(mOrder+1, 0.0);
}

void IIRFilter::printCoefficients() {
    std::cout << "a coefficients: ";
    for (const auto& coeff : _a) {
        std::cout << coeff << " ";
    }
    std::cout << std::endl;

    std::cout << "b coefficients: ";
    for (const auto& coeff : _b) {
        std::cout << coeff << " ";
    }
    std::cout << std::endl;
}


/* -------------------------------------------------------------------------- */
/* Méthode 1 */

void IIRFilter::setup() {
    if (!mIsSetup) {

        switch (mAnalogFilter) {
            case AnalogFilter::BESSEL:
                //BesselAnalogPrototype(mOrder, _z, _p, _k);
                break;
            case AnalogFilter::BUTTERWORTH:
                ButterworthCoefficients();
                break;
            case AnalogFilter::CHEBYSHEV1:
                //Chebyshev1AnalogPrototype(mOrder, _z, _p, _k);
                break;
            case AnalogFilter::CHEBYSHEV2:
                //Chebyshev2AnalogPrototype(mOrder, _z, _p, _k);
                break;
            case AnalogFilter::ELLIPTIC:
                //EllipticAnalogPrototype(mOrder, _z, _p, _k);
                break;
            default:
                throw std::invalid_argument("Unknown analog filter type");
        }

        mIsSetup = true;
    }
}

Signal IIRFilter::process(const Signal &input) {
    if (mIsSetup == false) {
        throw std::invalid_argument("Filter is not set up");
    }

    Signal output(input.size(), mFs);
    for (size_t i = 0; i < input.size(); i++) {
        output[i] = _b[0] * input[i];
        for (size_t j = 1; j < 3*static_cast<size_t>(mOrder); j++) {
            if (i >= j) {
                output[i] += (_b[j] * input[i - j]) - (_a[j] * output[i - j]);
            }
        }
    }
    
    mIsSetup = true;
    return output;
}

void IIRFilter::ButterworthCoefficients()
{
    double wc1 = 2 * M_PI * mFc1 / mFs;
    double wc2 = 2 * M_PI * mFc2 / mFs;
    double B, W0;

    if (mGabarit == FilterGabarit::BAND_PASS || mGabarit == FilterGabarit::BAND_STOP) {
        B = wc2 - wc1;
        W0 = sqrt(wc1 * wc2);
    }

    // Variables intermédiaires pour les coefficients
    double a0, a1, a2, b0, b1, b2;
    double a0_num, a1_num, a2_num;

    _a.resize(3 * mOrder);
    _b.resize(3 * mOrder);

    bool stable = true;

    for (int i = 0; i < mOrder; i++) {
        double theta = M_PI * (2.0 * i + 1.0) / (2.0 * mOrder);
        double pole_re = -sin(theta);
        double pole_im = cos(theta);

        if (mGabarit == FilterGabarit::LOW_PASS) { // Passe-bas
            double k = tan(wc1 / 2.0);
            a0 = 1.0 + sqrt(2.0) * k + k * k;
            a1 = 2.0 * (k * k - 1.0);
            a2 = 1.0 - sqrt(2.0) * k + k * k;
            b0 = k * k;
            b1 = 2.0 * b0;
            b2 = b0;
        } else if (mGabarit == FilterGabarit::HIGH_PASS) { // Passe-haut
            double k = tan(wc1 / 2.0);
            a0 = 1.0 + sqrt(2.0) * k + k * k;
            a1 = 2.0 * (k * k - 1.0);
            a2 = 1.0 - sqrt(2.0) * k + k * k;
            b0 = 1.0;
            b1 = -2.0;
            b2 = 1.0;
        } else if (mGabarit == FilterGabarit::BAND_PASS) { // Passe-bande
            double k1 = 2 * tan(B / 2.0);
            double k2 = W0 * W0;
            a0 = k2 + k1 * pole_re + pole_re * pole_re + pole_im * pole_im;
            a1 = 2.0 * (pole_re * pole_re + pole_im * pole_im - k2);
            a2 = k2 - k1 * pole_re + pole_re * pole_re + pole_im * pole_im;
            b0 = k1;
            b1 = 0.0;
            b2 = -k1;
        } else if (mGabarit == FilterGabarit::BAND_STOP) { // Coupe-bande
            double k1 = 2 * tan(B / 2.0);
            double k2 = W0 * W0;
            a0 = 1.0 + k1 * pole_re + k2;
            a1 = 2.0 * (k2 - 1.0);
            a2 = 1.0 - k1 * pole_re + k2;
            b0 = 1.0;
            b1 = -2.0 * (1.0 - k2);
            b2 = 1.0;
        }

        // Vérifier la stabilité (pour les types de filtre autres que passe-bas et passe-haut)
        if (mGabarit == FilterGabarit::BAND_PASS || mGabarit == FilterGabarit::BAND_STOP) {
            double zero_re = -cos(theta);
            double zero_im = sin(theta);
            double zero_radius = sqrt(zero_re * zero_re + zero_im * zero_im);
            if (zero_radius >= 1.0) {
                stable = false;
                break; // Sortir de la boucle si le filtre est instable
            }
        }

        // Normalisation des coefficients
        a0_num = a0;
        a1_num = a1 / a0_num;
        a2_num = a2 / a0_num;

        _b[3 * i] = b0 / a0_num;
        _b[3 * i + 1] = b1 / a0_num;
        _b[3 * i + 2] = b2 / a0_num;

        _a[3 * i] = 1.0;
        _a[3 * i + 1] = a1_num;
        _a[3 * i + 2] = a2_num;
    }

    if (!stable) {
        throw std::invalid_argument("Le filtre est instable pour les paramètres donnés");
    }
}



/* -------------------------------------------------------------------------- */
/* Méthode 2 */

void ButterworthAnalogPrototype(int N, std::vector<std::complex<double>> &z, std::vector<std::complex<double>> &p, double &k) {
    p.clear();
    for (int i = 0; i < N; ++i) {
        double theta = M_PI * (2 * i + 1) / (2 * N);
        std::complex<double> pole = std::polar(1.0, theta);
        p.push_back(pole);
    }
    k = 1;
}

void TransformLowpassToLowpass(std::vector<std::complex<double>> &zeros, std::vector<std::complex<double>> &poles, double &k, double warped) {
    for (std::complex<double> &z : zeros) {
        z *= warped;
    }
    for (std::complex<double> &p : poles) {
        p *= warped;
    }
    int degree = poles.size() - zeros.size();
    k *= pow(warped, degree);
}

void TransformLowpassToHighpass(std::vector<std::complex<double>> &zeros, std::vector<std::complex<double>> &poles, double &k, double warped) {
    std::complex<double> mulZ = 1.0;
    for (std::complex<double> &z : zeros) {
        mulZ *= -z;
        z = warped / z;
    }
    std::complex<double> mulP = 1.0;
    for (std::complex<double> &p : poles) {
        mulP *= -p;
        p = warped / p;
    }
    //int degree = poles.size() - zeros.size();
    k *= std::real(mulZ / mulP);
}

void TransformLowpassToBandpass(std::vector<std::complex<double>> &zeros, std::vector<std::complex<double>> &poles, double &k, double warped1, double warped2) {
    double bw = warped2 - warped1;
    double wo = sqrt(warped1 * warped2);
    int degree = poles.size() - zeros.size();

    std::vector<std::complex<double>> z_bp, p_bp;
    for (const auto &z : zeros) {
        z_bp.push_back(z * bw / 2.0);
    }
    for (const auto &p : poles) {
        p_bp.push_back(p * bw / 2.0);
    }

    std::vector<std::complex<double>> z_final, p_final;
    for (const auto &z : z_bp) {
        std::complex<double> sqrt_term = std::sqrt(z * z - wo * wo);
        z_final.push_back(z + sqrt_term);
        z_final.push_back(z - sqrt_term);
    }
    for (const auto &p : p_bp) {
        std::complex<double> sqrt_term = std::sqrt(p * p - wo * wo);
        p_final.push_back(p + sqrt_term);
        p_final.push_back(p - sqrt_term);
    }

    for (int i = 0; i < degree; ++i) {
        z_final.push_back(0.0);
    }

    k *= pow(bw, degree);
    zeros = std::move(z_final);
    poles = std::move(p_final);
}


void TransformLowpassToBandStop(std::vector<std::complex<double>> &zeros, std::vector<std::complex<double>> &poles, double &k, double warped1, double warped2) {
    double bw = warped2 - warped1;
    double wo = sqrt(warped1 * warped2);
    int degree = poles.size() - zeros.size();

    std::vector<std::complex<double>> z_hp, p_hp;
    for (const auto &z : zeros) {
        z_hp.push_back((bw / 2.0) / z);
    }
    for (const auto &p : poles) {
        p_hp.push_back((bw / 2.0) / p);
    }

    std::vector<std::complex<double>> z_final, p_final;
    for (const auto &z : z_hp) {
        std::complex<double> sqrt_term = std::sqrt(z * z - wo * wo);
        z_final.push_back(z + sqrt_term);
        z_final.push_back(z - sqrt_term);
    }
    for (const auto &p : p_hp) {
        std::complex<double> sqrt_term = std::sqrt(p * p - wo * wo);
        p_final.push_back(p + sqrt_term);
        p_final.push_back(p - sqrt_term);
    }

    for (int i = 0; i < degree; ++i) {
        z_final.push_back(std::complex<double>(0.0, wo));
        z_final.push_back(std::complex<double>(0.0, -wo));
    }

    k *= std::real(std::accumulate(poles.begin(), poles.end(), std::complex<double>(1.0), 
                      [=](const std::complex<double> &accum, const std::complex<double> &pole) { return accum * (-pole); }) /
                   std::accumulate(zeros.begin(), zeros.end(), std::complex<double>(1.0), 
                      [=](const std::complex<double> &accum, const std::complex<double> &zero) { return accum * (-zero); }));

    zeros = std::move(z_final);
    poles = std::move(p_final);
}

void IIRBilinearTransformation(std::vector<std::complex<double>> &zeros, std::vector<std::complex<double>> &poles, double &k, double fs) {
    int degree = poles.size() - zeros.size();
    double fs2 = 2.0 * fs;

    std::vector<std::complex<double>> z_transformed, p_transformed;
    std::complex<double> mulZ = 1.0;
    for (const auto &z : zeros) {
        mulZ *= (std::complex<double>(fs2) - z);
        z_transformed.push_back((fs2 + z) / (fs2 - z));
    }
    std::complex<double> mulP = 1.0;
    for (const auto &p : poles) {
        mulP *= (std::complex<double>(fs2) - p);
        p_transformed.push_back((fs2 + p) / (fs2 - p));
    }
    for (int i = 0; i < degree; ++i) {
        z_transformed.push_back(std::complex<double>(-1.0, 0.0));
    }

    k *= std::real(mulZ / mulP);
    zeros = std::move(z_transformed);
    poles = std::move(p_transformed);
}


void poly(const std::vector<std::complex<double>> &roots, std::vector<double> &coeffs) {
    coeffs[0] = 1.0;
    for (const auto &root : roots) {
        std::vector<double> new_coeffs(coeffs.size() + 1);
        for (size_t i = 0; i < coeffs.size(); ++i) {
            new_coeffs[i] -= root.real() * coeffs[i];
            new_coeffs[i + 1] += coeffs[i];
        }
        coeffs = std::move(new_coeffs);
    }
}

void PolynomialTransfer(const std::vector<std::complex<double>> &zeros, const std::vector<std::complex<double>> &poles, double k, std::vector<double> &b, std::vector<double> &a) {
    poly(zeros, b);
    for (auto &coeff : b) {
        coeff *= k;
    }
    poly(poles, a);

    auto are_conjugates = [](const std::vector<std::complex<double>> &roots) {
        std::vector<std::complex<double>> pos_roots, neg_roots;
        for (const auto &root : roots) {
            if (root.imag() > 0) pos_roots.push_back(root);
            else if (root.imag() < 0) neg_roots.push_back(std::conj(root));
        }
        if (pos_roots.size() != neg_roots.size()) return false;
        std::sort(pos_roots.begin(), pos_roots.end(), [](const std::complex<double> &a, const std::complex<double> &b) { return std::arg(a) < std::arg(b); });
        std::sort(neg_roots.begin(), neg_roots.end(), [](const std::complex<double> &a, const std::complex<double> &b) { return std::arg(a) < std::arg(b); });
        return std::equal(pos_roots.begin(), pos_roots.end(), neg_roots.begin());
    };

    if (are_conjugates(zeros)) {
        for (auto &coeff : b) coeff = std::real(coeff);
    }
    if (are_conjugates(poles)) {
        for (auto &coeff : a) coeff = std::real(coeff);
    }
}

void IIRFilter::setup2() {
    if (!mIsSetup) {
        // Reset coefficient buffers
        _a.resize(mOrder+1, 0.0);
        _b.resize(mOrder+1, 0.0);
        // Reset zeros, poles and gian
        _z.clear();
        _p.clear();
        _k = 0;

        switch (mAnalogFilter) {
            case AnalogFilter::BESSEL:
                //BesselAnalogPrototype(mOrder, _z, _p, _k);
                break;
            case AnalogFilter::BUTTERWORTH:
                ButterworthAnalogPrototype(mOrder, _z, _p, _k);
                break;
            case AnalogFilter::CHEBYSHEV1:
                //Chebyshev1AnalogPrototype(mOrder, _z, _p, _k);
                break;
            case AnalogFilter::CHEBYSHEV2:
                //Chebyshev2AnalogPrototype(mOrder, _z, _p, _k);
                break;
            case AnalogFilter::ELLIPTIC:
                //EllipticAnalogPrototype(mOrder, _z, _p, _k);
                break;
            default:
                throw std::invalid_argument("Unknown analog filter type");
        }

        double fs = 2.0;
        double warped1 = 2.0 * fs* tan(M_PI * mFc1 / fs);
        double warped2 = 2.0 * fs* tan(M_PI * mFc2 / fs);

        if (mGabarit == FilterGabarit::LOW_PASS) {
            TransformLowpassToLowpass(_z, _p, _k, warped1);
        } else if (mGabarit == FilterGabarit::HIGH_PASS) {
            TransformLowpassToHighpass(_z, _p, _k, warped1);
        } else if (mGabarit == FilterGabarit::LOW_PASS) {
            TransformLowpassToBandpass(_z, _p, _k, warped1, warped2);
        } else if (mGabarit == FilterGabarit::HIGH_PASS) {
            TransformLowpassToBandStop(_z, _p, _k, warped1, warped2);
        }

        IIRBilinearTransformation(_z, _p, _k, fs);
        PolynomialTransfer(_z, _p, _k, _a, _b);

        mIsSetup = true;
    }
}

double IIRFilter::eqdiff(double x) {
    if (mIsSetup) {
        _xMem[0] = x;
        double y = 0.0;
        for (int i = 0; i <= mOrder; ++i) {
            y += _b[i] * _xMem[i];
            if (i > 0) {
                y -= _a[i] * _yMem[i];
            }
        }
        for (int i = mOrder; i > 0; --i) {
            _xMem[i] = _xMem[i - 1];
            _yMem[i] = _yMem[i - 1];
        }
        _yMem[0] = y;
        return y;
    } else {
        return x;
    }
}