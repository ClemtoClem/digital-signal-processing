#ifndef __SIGNAL_HPP
#define __SIGNAL_HPP

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <optional>
#include <string>

enum class WaveformType {
    SINUS,
    COSINUS,
    SQUARE,   // with duty cycle
    TRIANGLE, // with duty cycle
    RAMP_UP,
    RAMP_DOWN,
    POSITIVE_DC,
    NEGATIVE_DC,
};

class Signal : public std::vector<std::complex<double>> {
private:
    double samplingFrequency; // Fréquence d'échantillonnage en Hz
    std::optional<std::string> name; // Nom optionnel du signal

public:
    // Constructeur par défaut
    Signal();

    // Constructeur prenant une taille initiale, la fréquence d'échantillonnage et éventuellement le nom du signal
    Signal(size_t size, double samplingFrequency, const std::string& name = std::string());

    Signal(const std::vector<std::complex<double>>& values, double samplingFrequency, const std::string& name = std::string());

    // Constructeur par recopie
    Signal(const Signal& other);

    // Opérateur d'affectation (ne copie pas le nom du signal)
    Signal& operator=(const Signal& other);

    /* ------------------------------- */

    // Méthode pour obtenir le nom du signal (retourne une chaîne vide si aucun nom n'est défini)
    const std::string &getName() const;

    // Méthode pour définir le nom du signal
    void setName(const std::string &newName);

    // Retourne la fréquence d'échantillonnage
    double getSamplingFrequency() const;

    std::vector<double> getRealBuffer() const;

    /* ------------------------------- */

    // Surcharge de l'opérateur d'addition
    Signal operator+(const Signal& other) const;

    // Surcharge de l'opérateur de soustraction
    Signal operator-(const Signal& other) const;

    // Surcharge de l'opérateur de multiplication
    Signal operator*(const Signal& other) const;

    // Surcharge de l'opérateur de division
    Signal operator/(const Signal& other) const;

    // Surcharge de l'opérateur +=
    Signal& operator+=(const Signal& other);

    // Surcharge de l'opérateur -=
    Signal& operator-=(const Signal& other);

    // Surcharge de l'opérateur *=
    Signal& operator*=(const Signal& other);

    // Surcharge de l'opérateur /=
    Signal& operator/=(const Signal& other);

    /* ------------------------------- */

    // Surcharge de l'opérateur d'addition avec un complexe
    Signal operator+(const std::complex<double>& value) const;

    // Surcharge de l'opérateur de soustraction avec un complexe
    Signal operator-(const std::complex<double>& value) const;

    // Surcharge de l'opérateur de multiplication avec un complexe
    Signal operator*(const std::complex<double>& value) const;

    // Surcharge de l'opérateur de division avec un complexe
    Signal operator/(const std::complex<double>& value) const;

    // Surcharge de l'opérateur += avec un complexe
    Signal& operator+=(const std::complex<double>& value);

    // Surcharge de l'opérateur -= avec un complexe
    Signal& operator-=(const std::complex<double>& value);

    // Surcharge de l'opérateur *= avec un complexe
    Signal& operator*=(const std::complex<double>& value);

    // Surcharge de l'opérateur /= avec un complexe
    Signal& operator/=(const std::complex<double>& value);

    /* ------------------------------- */

    // Fonction amie pour l'addition d'un complexe et d'un signal
    friend Signal operator+(const std::complex<double>& value, const Signal& signal) {
        return signal + value;
    }

    // Fonction amie pour la soustraction d'un complexe et d'un signal
    friend Signal operator-(const std::complex<double>& value, const Signal& signal) {
        Signal result(signal.size(), signal.samplingFrequency);
        for (size_t i = 0; i < signal.size(); ++i) {
            result[i] = value - signal[i];
        }
        return result;
    }

    // Fonction amie pour la multiplication d'un complexe et d'un signal
    friend Signal operator*(const std::complex<double>& value, const Signal& signal) {
        return signal * value;
    }

    // Fonction amie pour la division d'un complexe et d'un signal
    friend Signal operator/(const std::complex<double>& value, const Signal& signal) {
        Signal result(signal.size(), signal.samplingFrequency);
        for (size_t i = 0; i < signal.size(); ++i) {
            result[i] = value / signal[i];
        }
        return result;
    }

    /* ------------------------------- */

    friend Signal cos(const Signal& signal) {
        Signal result(signal.size(), signal.samplingFrequency);
        for (size_t i = 0; i < signal.size(); ++i) {
            result[i] = std::cos(signal[i]);
        }
        return result;
    }

    friend Signal sin(const Signal& signal) {
        Signal result(signal.size(), signal.samplingFrequency);
        for (size_t i = 0; i < signal.size(); ++i) {
            result[i] = std::sin(signal[i]);
        }
        return result;
    }

    friend Signal tan(const Signal& signal) {
        Signal result(signal.size(), signal.samplingFrequency);
        for (size_t i = 0; i < signal.size(); ++i) {
            result[i] = std::tan(signal[i]);
        }
        return result;
    }

    friend Signal cosh(const Signal& signal) {
        Signal result(signal.size(), signal.samplingFrequency);
        for (size_t i = 0; i < signal.size(); ++i) {
            result[i] = std::cosh(signal[i]);
        }
        return result;
    }

    friend Signal sinh(const Signal& signal) {
        Signal result(signal.size(), signal.samplingFrequency);
        for (size_t i = 0; i < signal.size(); ++i) {
            result[i] = std::sinh(signal[i]);
        }
        return result;
    }

    friend Signal tanh(const Signal& signal) {
        Signal result(signal.size(), signal.samplingFrequency);
        for (size_t i = 0; i < signal.size(); ++i) {
            result[i] = std::tanh(signal[i]);
        }
        
        return result;
    }

    friend Signal acos(const Signal& signal) {
        Signal result(signal.size(), signal.samplingFrequency);
        for (size_t i = 0; i < signal.size(); ++i) {
            result[i] = std::acos(signal[i]);
        }
        return result;
    }

    friend Signal asin(const Signal& signal) {
        Signal result(signal.size(), signal.samplingFrequency);
        for (size_t i = 0; i < signal.size(); ++i) {
            result[i] = std::asin(signal[i]);
        }
        return result;
    }

    friend Signal atan(const Signal& signal) {
        Signal result(signal.size(), signal.samplingFrequency);
        for (size_t i = 0; i < signal.size(); ++i) {
            result[i] = std::atan(signal[i]);
        }
        
        return result;
    }

    friend Signal acosh(const Signal& signal) {
        Signal result(signal.size(), signal.samplingFrequency);
        for (size_t i = 0; i < signal.size(); ++i) {
            result[i] = std::acosh(signal[i]);
        }
        return result;
    }

    friend Signal asinh(const Signal& signal) {
        Signal result(signal.size(), signal.samplingFrequency);
        for (size_t i = 0; i < signal.size(); ++i) {
            result[i] = std::asinh(signal[i]);
        }
        return result;
    }

    friend Signal atanh(const Signal& signal) {
        Signal result(signal.size(), signal.samplingFrequency);
        for (size_t i = 0; i < signal.size(); ++i) {
            result[i] = std::atanh(signal[i]);
        }
        
        return result;
    }

    /* ------------------------------- */

    // Générer une forme périodique au signal
    void generateWaveform(WaveformType type, double amplitude, double frequency,
        double phase = 0.0, double offset = 0.0, double delay = 0.0, double duty_cycle = 50.0);
    
    /* ------------------------------- */

    // Fonction pour mettre au carré chaque élément du signal
    Signal square() const;

    // Fonction pour élever chaque élément du signal à une puissance donnée
    Signal pow(double exponent) const;

    // Fonction pour calculer la racine carrée de chaque élément du signal
    Signal sqrt() const;

    // Fonction pour calculer le logarithme naturel de chaque élément du signal
    Signal log() const;

    // Fonction pour calculer le logarithme décimal de chaque élément du signal
    Signal log10() const;

    // Fonction pour calculer l'exponentielle de chaque élément du signal
    Signal exp() const;

    /* ------------------------------- */

    // Fonction pour calculer le maximum du signal
    std::complex<double> max() const;

    // Fonction pour calculer le minimum du signal
    std::complex<double> min() const;

    // Fonction pour calculer la moyenne du signal
    std::complex<double> mean() const;
    
    /* ------------------------------- */

    // Fonction pour effectuer la transformée de Fourier discrète  rapide(DFT) et générer le spectre
    Signal DFT(size_t size_zero_padding = 0) const;

    // Fonction pour effectuer la transformée de Fourier inverse rapide (IDFT) et reconstruire le signal
    Signal IDFT(size_t size_zero_padding = 0) const;

    // Générer un axe fréquentiel en fonction de la fréquence d'échantillonnage
    void generateFrequencyAxis();

    /* ------------------------------- */
    // Fonction pour afficher le signal
    void display() const;
};

#endif // __SIGNAL_HPP