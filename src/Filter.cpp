#include "Filter.hpp"
#include <algorithm>

static int relative_degree(const std::vector<complexd> &zeros, const std::vector<complexd> &poles) {
    int degree = poles.size() - zeros.size();
    if (degree < 0)
        throw std::invalid_argument("Improper transfer function. Must have at least as many poles as zeros.");
    return degree;
}

static void lowpassToLowpass(std::vector<complexd> &zeros, std::vector<complexd> &poles, double &k, double wo) {
    // Calculer le degré relatif
    int degree = relative_degree(zeros, poles);

    // Transformer les pôles et les zéros en les multipliant par wo
    for (auto &zero : zeros) zero *= wo;
    for (auto &pole : poles) pole *= wo;

    // Ajuster le gain
    k *= std::pow(wo, degree);
}

static void lowpassToHighpass(std::vector<complexd> &zeros, std::vector<complexd> &poles, double &k, double wo) {
    // Calculer le degré relatif
    int degree = relative_degree(zeros, poles);

    // Transformer les pôles et les zéros pour passe-haut
    for (auto &zero : zeros) zero = wo / zero;
    for (auto &pole : poles) pole = wo / pole;

    // Ajouter des zéros à l'origine pour les zéros à l'infini
    for (int i = 0; i < degree; ++i) zeros.push_back(0.0);

    // Compenser la variation de gain
    std::complex<double> prod_z = 1.0;
    for (const auto &zero : zeros) {
        if (zero != 0.0) {
            prod_z *= -zero;
        }
    }

    std::complex<double> prod_p = 1.0;
    for (const auto &pole : poles) {
        prod_p *= -pole;
    }

    k *= std::real(prod_z / prod_p);
}

static void lowpassToBandpass(std::vector<complexd> &zeros, std::vector<complexd> &poles, double &k, double wo, double bw) {
    // Calculer le degré relatif
    int degree = relative_degree(zeros, poles);

    // Transformer les pôles et les zéros pour passe-bande
    std::vector<complexd> z_bp, p_bp;

    for (const auto &zero : zeros) {
        complexd z_lp = zero * (bw / 2.0);
        z_bp.push_back(z_lp + std::sqrt(z_lp * z_lp - wo * wo));
        z_bp.push_back(z_lp - std::sqrt(z_lp * z_lp - wo * wo));
    }

    for (const auto &pole : poles) {
        complexd p_lp = pole * (bw / 2.0);
        p_bp.push_back(p_lp + std::sqrt(p_lp * p_lp - wo * wo));
        p_bp.push_back(p_lp - std::sqrt(p_lp * p_lp - wo * wo));
    }

    // Ajouter des zéros à l'origine pour les zéros à l'infini
    for (int i = 0; i < degree; ++i) z_bp.push_back(0.0);

    // Mettre à jour les vecteurs de zéros et de pôles
    zeros = z_bp;
    poles = p_bp;

    // Compenser la variation de gain
    k *= std::pow(bw, degree);
}

static void lowpassToBandStop(std::vector<complexd> &zeros, std::vector<complexd> &poles, double &k, double wo, double bw) {
    // Invert to a highpass filter with desired bandwidth
    std::vector<complexd> z_hp(zeros.size()), p_hp(poles.size());
    for (size_t i = 0; i < zeros.size(); ++i) {
        z_hp[i] = (bw / 2.0) / zeros[i];
    }
    for (size_t i = 0; i < poles.size(); ++i) {
        p_hp[i] = (bw / 2.0) / poles[i];
    }

    // Duplicate poles and zeros and shift from baseband to +wo and -wo
    std::vector<complexd> z_bs, p_bs;
    for (size_t i = 0; i < z_hp.size(); ++i) {
        z_bs.push_back(z_hp[i] + sqrt(z_hp[i] * z_hp[i] - wo * wo));
        z_bs.push_back(z_hp[i] - sqrt(z_hp[i] * z_hp[i] - wo * wo));
    }
    for (size_t i = 0; i < p_hp.size(); ++i) {
        p_bs.push_back(p_hp[i] + sqrt(p_hp[i] * p_hp[i] - wo * wo));
        p_bs.push_back(p_hp[i] - sqrt(p_hp[i] * p_hp[i] - wo * wo));
    }

    // Move any zeros that were at infinity to the center of the stopband
    int degree = poles.size() - zeros.size();
    for (int i = 0; i < degree; ++i) {
        z_bs.push_back(std::complex<double>(0.0, wo));
        z_bs.push_back(std::complex<double>(0.0, -wo));
    }

    // Cancel out gain change caused by inversion
    std::complex<double> prod_z = 1.0, prod_p = 1.0;
    for (const auto &zero : zeros) prod_z *= -zero;
    for (const auto &pole : poles) prod_p *= -pole;
    k *= std::real(prod_z / prod_p);
}

/* =========================================================================================== */

// return numerator and denominator
void polynomial_transfer(const std::vector<complexd> &zeros, const std::vector<complexd> &poles, double k, std::vector<double> &numerator, std::vector<double> &denominator) {
    int num_zeros = zeros.size();
    int num_poles = poles.size();

    // Calculating the numerator polynomial coefficients
    if (num_zeros > 1) {
        numerator.resize(num_zeros + 1);
        for (int i = 0; i < num_zeros; ++i) {
            std::vector<double> temp(num_zeros + 1, 0.0);
            temp[i + 1] = 1.0;
            for (int j = 0; j < num_zeros; ++j) {
                if (j != i) {
                    for (int k = num_zeros; k > 0; --k) {
                        temp[k] = temp[k] * -zeros[j].real() + temp[k - 1];
                    }
                    temp[0] *= -zeros[j].real();
                }
            }
            for (int k = 0; k <= num_zeros; ++k) {
                numerator[k] += k * k * temp[k];
            }
        }
    } else if (num_zeros == 1) {
        numerator.resize(2);
        numerator[0] = 0.0;
        numerator[1] = k * (-zeros[0].real());
    } else {
        numerator.resize(1);
        numerator[0] = k;
    }

    // Calculating the denominator polynomial coefficients
    denominator.resize(num_poles + 1);
    denominator[0] = 1.0;
    for (int i = 0; i < num_poles; ++i) {
        for (int j = num_poles; j > 0; --j) {
            denominator[j] = denominator[j] * -poles[i].real() + denominator[j - 1];
        }
        denominator[0] *= -poles[i].real();
    }
}

/* =========================================================================================== */

// Forward filtering function
Signal forward_filter(const Signal &signal, const std::vector<double> &numerator, const std::vector<double> &denominator) {
    int num_signal = signal.size();
    int num_numerator = numerator.size();
    int num_denominator = denominator.size();

    Signal result(signal.size(), signal.getSampleFrequency(), 0.0);

    for (int i = num_numerator - 1; i < num_signal; ++i) {
        double sum = 0.0;
        for (int j = 0; j < num_numerator; ++j) {
            int idx = i - j;
            if (idx >= 0) {
                sum += numerator[j] * signal[idx].real();
            }
        }
        for (int j = 1; j < num_denominator; ++j) {
            int idx = i - j;
            if (idx >= 0) {
                sum -= denominator[j] * result[idx].real();
            }
        }
        result[i].real(sum / denominator[0]);
    }

    return result;
}


// Backward filtering function
Signal backward_filter(const Signal &signal, const std::vector<double> &numerator, const std::vector<double> &denominator) {
    int num_signal = signal.size();
    int num_numerator = numerator.size();
    int num_denominator = denominator.size();

    Signal result(signal.size(), signal.getSampleFrequency(), 0.0);

    for (int i = num_signal - num_numerator; i >= 0; --i) {
        double sum = 0.0;
        for (int j = 0; j < num_numerator; ++j) {
            sum += numerator[j] * signal[i + j].real();
        }
        for (int j = 1; j < num_denominator; ++j) {
            sum -= denominator[j] * result[i + j].real();
        }
        result[i].real(sum / denominator[0]);
    }

    return result;
}

/* =========================================================================================== */


void IIRFilter::butterworthLowPass(int order, double fc, int fs) {
    // Initialisation des coefficients à zéro
    m_numerator.clear();
    m_denominator.clear();

    // Calcul des pôles du filtre Butterworth
    std::vector<std::complex<double>> poles;
    for (int i = 1; i <= order; ++i) {
        double theta = M_PI * (2.0 * i + order - 1) / (2.0 * order);
        std::complex<double> pole(-sin(theta), cos(theta));
        poles.push_back(pole);
    }

    // Calcul des coefficients du numérateur (b) en multipliant les termes (s - p) pour chaque pôle
    m_numerator.push_back(1.0);
    for (const auto& pole : poles) {
        std::vector<double> temp;
        temp.push_back(1.0);
        temp.push_back(-2.0 * pole.real());
        temp.push_back(pole.real() * pole.real() + pole.imag() * pole.imag());
        m_numerator = convolve(m_numerator, temp);
    }

    // Calcul des coefficients du dénominateur (a) à partir des pôles
    m_denominator = std::vector<double>(order + 1, 0.0);
    m_denominator[0] = 1.0;
    for (const auto& pole : poles) {
        std::vector<double> temp;
        temp.push_back(1.0);
        temp.push_back(-2.0 * pole.real());
        temp.push_back(pole.real() * pole.real() + pole.imag() * pole.imag());
        m_denominator = convolve(m_denominator, temp);
    }

    // Normalisation des coefficients du numérateur par le premier coefficient du dénominateur
    double norm_factor = m_denominator[0];
    for (auto& coef : m_numerator) {
        coef /= norm_factor;
    }
    for (auto& coef : m_denominator) {
        coef /= norm_factor;
    }

    // Configuration terminée
    m_isSetup = true;
}

void IIRFilter::butterworthHighPass(int order, double fc, int fs) {
    // Calcul de la fréquence de Nyquist
    double f_nyq = fs / 2.0;

    // Conversion de la fréquence de coupure en rapport à la fréquence de Nyquist
    double fc_ratio = fc / f_nyq;

    // Calcul des coefficients du filtre Butterworth
    double theta_c = M_PI * fc_ratio;
    double d = 1.0 / std::tan(theta_c / 2.0);
    double d_square = d * d;

    // Calcul des coefficients du polynôme de Butterworth
    m_numerator.resize(order + 1, 0.0);
    m_denominator.resize(order + 1, 0.0);

    m_numerator[0] = 1.0;

    for (int i = 1; i <= order / 2; ++i) {
        double theta = M_PI * (2 * i + order - 1) / (2.0 * order);
        double beta = 0.5 * (1.0 - std::sin(theta)) / (1.0 + std::sin(theta));
        double gamma = (0.5 + beta) * std::cos(theta);

        m_numerator[i] = gamma / (d_square + gamma);
        m_numerator[order - i] = m_numerator[i];
        m_denominator[i] = -2.0 * d * std::cos(M_PI * (2 * i + order - 1) / (2.0 * order)) / (1.0 + beta * (1.0 - d_square));
        m_denominator[order - i] = -m_denominator[i];
    }

    // Configuration terminée
    m_isSetup = true;
}

void IIRFilter::butterworthBandPass(int order, double fc1, double fc2, int fs) {
    // Calcul de la fréquence de Nyquist
    double f_nyq = fs / 2.0;

    // Conversion des fréquences de coupure en rapport à la fréquence de Nyquist
    double fc1_ratio = fc1 / f_nyq;
    double fc2_ratio = fc2 / f_nyq;

    // Calcul des coefficients du filtre Butterworth passe-bande
    double wc1 = 2.0 * M_PI * fc1_ratio;
    double wc2 = 2.0 * M_PI * fc2_ratio;
    double wc_ratio = wc2 - wc1;
    double wc_prod = wc1 * wc2;
    double sqrt_wc_prod = sqrt(wc_prod);

    // Calcul des coefficients du polynôme de Butterworth
    m_numerator.resize(order + 1, 0.0);
    m_denominator.resize(order + 1, 0.0);

    m_numerator[0] = sqrt_wc_prod;

    for (int i = 1; i <= order / 2; ++i) {
        double theta = M_PI * (2 * i + order - 1) / (2.0 * order);
        double beta = 0.5 * (1.0 - std::sin(theta)) / (1.0 + std::sin(theta));
        double gamma = (0.5 + beta) * std::cos(theta);

        m_numerator[i] = wc_prod / (gamma + wc_prod);
        m_numerator[order - i] = m_numerator[i];
        m_denominator[i] = -2.0 * wc_ratio * std::cos(M_PI * (2 * i + order - 1) / (2.0 * order)) / (1.0 + beta * (1.0 - wc_prod));
        m_denominator[order - i] = -m_denominator[i];
    }

    // Configuration terminée
    m_isSetup = true;
}

void IIRFilter::butterworthBandStop(int order, double fc1, double fc2, int fs) {
    // Normalisation des fréquences de coupure
    double wc1 = 2.0 * M_PI * fc1 / fs;
    double wc2 = 2.0 * M_PI * fc2 / fs;

    // Calcul des pôles et des zéros pour le filtre de Butterworth en passe-bande
    std::vector<complexd> zeros;
    std::vector<complexd> poles;
    double k = 1.0;

    for (int i = 1; i <= order; ++i) {
        double theta = M_PI * (2 * i - 1) / (2.0 * order);
        complexd pole1 = -sin(theta) / 2.0 + complexd(0.0, cos(theta)) * sqrt(3.0) / 2.0;
        complexd pole2 = -sin(theta) / 2.0 - complexd(0.0, cos(theta)) * sqrt(3.0) / 2.0;
        poles.push_back(exp(complexd(0.0, wc1)) * pole1);
        poles.push_back(exp(complexd(0.0, wc1)) * pole2);
        poles.push_back(exp(complexd(0.0, wc2)) * pole1);
        poles.push_back(exp(complexd(0.0, wc2)) * pole2);
    }

    // Transformations de fréquence pour obtenir la bande d'arrêt
    lowpassToBandStop(zeros, poles, k, wc1, wc2 - wc1);

    // Calcul des coefficients a_ et b_
    polynomial_transfer(zeros, poles, k, m_numerator, m_denominator);

    // Configuration terminée
    m_isSetup = true;
}

Signal IIRFilter::filter(const Signal& input) {
    if (!m_isSetup) {
        throw std::runtime_error("Filter is not set up.");
    }

    // Création d'une copie de l'entrée pour stocker le signal filtré
    Signal filtered_signal = input;

    // Application du filtrage en avant
    filtered_signal = forward_filter(filtered_signal, m_numerator, m_denominator);

    // Inversion du signal filtré
    std::reverse(filtered_signal.getBuffer().begin(), filtered_signal.getBuffer().end());

    // Application du filtrage en avant sur le signal inversé
    filtered_signal = forward_filter(filtered_signal, m_numerator, m_denominator);

    // Inversion du signal filtré à nouveau
    std::reverse(filtered_signal.getBuffer().begin(), filtered_signal.getBuffer().end());

    // Retourner le signal filtré
    return filtered_signal;
}

// Méthode pour afficher les coefficients du filtre (à des fins de débogage)
void IIRFilter::printCoefficients() {
    std::cout << "Numerator coefficients (b): ";
    for (const auto& coef : m_numerator) {
        std::cout << coef << " ";
    }
    std::cout << std::endl;

    std::cout << "Denominator coefficients (a): ";
    for (const auto& coef : m_denominator) {
        std::cout << coef << " ";
    }
    std::cout << std::endl;
}