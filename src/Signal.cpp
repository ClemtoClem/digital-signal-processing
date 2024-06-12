#include "Signal.hpp"

Signal::Signal() : samplingFrequency(1.0) {}

Signal::Signal(size_t size, double samplingFrequency, const std::string& name) :
    std::vector<std::complex<double>>(size), samplingFrequency(samplingFrequency), name(name) {
        //std::cout << name << ":" << this->size() << std::endl;
    }

Signal::Signal(const std::vector<std::complex<double>>& values, double samplingFrequency, const std::string &name)
    : std::vector<std::complex<double>>(values), samplingFrequency(samplingFrequency), name(name) {}

Signal::Signal(const Signal& other) :
    std::vector<std::complex<double>>(other), samplingFrequency(other.samplingFrequency), name(other.name) {}

Signal& Signal::operator=(const Signal& other) {
    if (this != &other) {
        std::vector<std::complex<double>>::operator=(other); // Copie les éléments du signal

        samplingFrequency = other.samplingFrequency; // Copie la fréquence d'échantillonnage
    }
    return *this;
}

/* ------------------------------- */

const std::string &Signal::getName() const {
    if (name) {
        return *name;
    } else {
        static const std::string emptyString = ""; // Valeur par défaut si name n'est pas initialisé
        return emptyString;
    }
}

void Signal::setName(const std::string& newName) {
    name = newName;
}

double Signal::getSamplingFrequency() const {
    return samplingFrequency;
}

/* ------------------------------- */

Signal Signal::operator+(const Signal& other) const {
    Signal result(std::max(this->size(), other.size()), this->samplingFrequency);
    for (size_t i = 0; i < result.size(); ++i) {
        std::complex<double> val1 = (i < this->size() ? (*this)[i] : 0);
        std::complex<double> val2 = (i < other.size() ? other[i] : 0);
        result[i] = val1 + val2;
    }
    return result;
}

Signal Signal::operator-(const Signal& other) const {
    Signal result(std::max(this->size(), other.size()), this->samplingFrequency);
    for (size_t i = 0; i < result.size(); ++i) {
        std::complex<double> val1 = (i < this->size() ? (*this)[i] : 0);
        std::complex<double> val2 = (i < other.size() ? other[i] : 0);
        result[i] = val1 - val2;
    }
    return result;
}

Signal Signal::operator*(const Signal& other) const {
    Signal result(std::max(this->size(), other.size()), this->samplingFrequency);
    for (size_t i = 0; i < result.size(); ++i) {
        std::complex<double> val1 = (i < this->size() ? (*this)[i] : 1);
        std::complex<double> val2 = (i < other.size() ? other[i] : 1);
        result[i] = val1 * val2;
    }
    return result;
}

Signal Signal::operator/(const Signal& other) const {
    Signal result(std::max(this->size(), other.size()), this->samplingFrequency);
    for (size_t i = 0; i < result.size(); ++i) {
        std::complex<double> val1 = (i < this->size() ? (*this)[i] : 1);
        std::complex<double> val2 = (i < other.size() && other[i] != std::complex<double>(0, 0) ? other[i] : 1);
        result[i] = val1 / val2;
    }
    return result;
}

Signal &Signal::operator+=(const Signal &other){
    for (size_t i = 0; i < std::max(this->size(), other.size()); ++i) {
        if (i >= this->size()) {
            this->push_back(other[i]);
        } else {
            (*this)[i] += (i < other.size() ? other[i] : 0);
        }
    }
    return *this;
}

Signal& Signal::operator-=(const Signal& other) {
    for (size_t i = 0; i < std::max(this->size(), other.size()); ++i) {
        if (i >= this->size()) {
            this->push_back(-other[i]);
        } else {
            (*this)[i] -= (i < other.size() ? other[i] : 0);
        }
    }
    return *this;
}

Signal& Signal::operator*=(const Signal& other) {
    for (size_t i = 0; i < std::max(this->size(), other.size()); ++i) {
        if (i >= this->size()) {
            this->push_back(0);
        } else {
            (*this)[i] *= (i < other.size() ? other[i] : 1);
        }
    }
    return *this;
}

Signal& Signal::operator/=(const Signal& other) {
    for (size_t i = 0; i < std::max(this->size(), other.size()); ++i) {
        if (i >= this->size()) {
            this->push_back(0);
        } else {
            (*this)[i] /= (i < other.size() && other[i] != std::complex<double>(0, 0) ? other[i] : 1);
        }
    }
    return *this;
}

/* ------------------------------ */

Signal Signal::operator+(const std::complex<double>& value) const {
    Signal result(this->size(), this->samplingFrequency);
    for (size_t i = 0; i < this->size(); ++i) {
        result[i] = (*this)[i] + value;
    }
    return result;
}

Signal Signal::operator-(const std::complex<double>& value) const {
    Signal result(this->size(), this->samplingFrequency);
    for (size_t i = 0; i < this->size(); ++i) {
        result[i] = (*this)[i] - value;
    }
    return result;
}

Signal Signal::operator*(const std::complex<double>& value) const {
    Signal result(this->size(), this->samplingFrequency);
    for (size_t i = 0; i < this->size(); ++i) {
        result[i] = (*this)[i] * value;
    }
    return result;
}

Signal Signal::operator/(const std::complex<double>& value) const {
    Signal result(this->size(), this->samplingFrequency);
    for (size_t i = 0; i < this->size(); ++i) {
        result[i] = (*this)[i] / value;
    }
    return result;
}

Signal& Signal::operator+=(const std::complex<double>& value) {
    for (auto& element : *this) {
        element += value;
    }
    return *this;
}

Signal& Signal::operator-=(const std::complex<double>& value) {
    for (auto& element : *this) {
        element -= value;
    }
    return *this;
}

Signal& Signal::operator*=(const std::complex<double>& value) {
    for (auto& element : *this) {
        element *= value;
    }
    return *this;
}

Signal& Signal::operator/=(const std::complex<double>& value) {
    for (auto& element : *this) {
        element /= value;
    }
    return *this;
}

/* ------------------------------- */

void Signal::generateWaveform(WaveformType type, double amplitude, double frequency, double phase, double offset, double delay, double duty_cycle) {
    // Vérification de la fréquence d'échantillonnage pour éviter la division par zéro
    if (samplingFrequency <= 0) {
        std::cerr << "Error: Sampling frequency is not positive." << std::endl;
        return;
    }

    // Calcul du période d'échantillonnage
    double sample_period = 1.0 / samplingFrequency;
    // Conversion du délai en nombre d'échantillons
    double delay_samples = delay / 1000.0 * samplingFrequency;

    // Variables pour le temps et la valeur du signal
    double t, value;
    for (size_t i = 0; i < this->size(); i++) {
        // Calcul du temps pour cet échantillon en tenant compte du décalage
        t = (i - delay_samples) * sample_period;
        value = 0;

        // Sélection du type de forme d'onde
        switch (type) {
            case WaveformType::SINUS:
                value = amplitude * std::sin(2 * M_PI * frequency * t + phase) + offset;
                break;
            case WaveformType::COSINUS:
                value = amplitude * std::cos(2 * M_PI * frequency * t + phase) + offset;
                break;
            case WaveformType::SQUARE:
                if (t >= 0) {
                    double period = 1.0 / frequency;
                    double time_in_period = fmod(t, period);
                    value = (time_in_period < duty_cycle * period) ? amplitude : -amplitude;
                    value += offset;
                }
                break;
            case WaveformType::TRIANGLE:
                if (t >= 0) {
                    double period = 1.0 / frequency;
                    double time_in_period = fmod(t, period);
                    double ramp_up_duration = duty_cycle * period;
                    if (time_in_period < ramp_up_duration) {
                        value = (2 * amplitude / ramp_up_duration) * time_in_period - amplitude;
                    } else {
                        value = (-2 * amplitude / (period - ramp_up_duration)) * (time_in_period - ramp_up_duration) + amplitude;
                    }
                    value += offset;
                }
                break;
            case WaveformType::RAMP_UP:
                if (t >= 0) {
                    double period = 1.0 / frequency;
                    double time_in_period = fmod(t, period);
                    value = (2 * amplitude / period) * time_in_period - amplitude + offset;
                }
                break;
            case WaveformType::RAMP_DOWN:
                if (t >= 0) {
                    double period = 1.0 / frequency;
                    double time_in_period = fmod(t, period);
                    value = (-2 * amplitude / period) * time_in_period + amplitude + offset;
                }
                break;
            case WaveformType::POSITIVE_DC:
                value = amplitude + offset;
                break;
            case WaveformType::NEGATIVE_DC:
                value = -amplitude + offset;
                break;
            default:
                std::cerr << "Error: Unknown waveform type." << std::endl;
                return;
        }

        // Stockage de la valeur du signal à cet échantillon
        (*this)[i] = std::complex<double>(value, 0);
    }
}

/* ------------------------------- */

Signal Signal::square() const
{
    Signal result(*this);
    for (auto& element : result) {
        element *= element;
    }
    return result;
}

Signal Signal::pow(double exponent) const {
    Signal result(*this);
    for (auto& element : result) {
        element = std::pow(element, exponent);
    }
    return result;
}

Signal Signal::sqrt() const {
    Signal result(*this);
    for (auto& element : result) {
        element = std::sqrt(element);
    }
    return result;
}

Signal Signal::log() const {
    Signal result(*this);
    for (auto& element : result) {
        element = std::log(element);
    }
    return result;
}

Signal Signal::log10() const {
    Signal result(*this);
    for (auto& element : result) {
        element = std::log10(element);
    }
    return result;
}

Signal Signal::exp() const {
    Signal result(*this);
    for (auto& element : result) {
        element = std::exp(element);
    }
    return result;
}

/* ------------------------------- */

std::complex<double> Signal::max() const {
    if (this->empty()) {
        throw std::runtime_error("Signal is empty");
    }
    return *std::max_element(this->begin(), this->end(), [](const std::complex<double>& a, const std::complex<double>& b) {
        return std::abs(a) < std::abs(b);
    });
}

std::complex<double> Signal::min() const {
    if (this->empty()) {
        throw std::runtime_error("Signal is empty");
    }
    return *std::min_element(this->begin(), this->end(), [](const std::complex<double>& a, const std::complex<double>& b) {
        return std::abs(a) < std::abs(b);
    });
}

std::complex<double> Signal::mean() const {
    if (this->empty()) {
        throw std::runtime_error("Signal is empty");
    }
    std::complex<double> sum = std::accumulate(this->begin(), this->end(), std::complex<double>(0.0, 0.0));
    return sum / static_cast<double>(this->size());
}

/* ------------------------------- */

Signal Signal::DFT() const
{
    const size_t N = size();
    Signal spectrum(N, samplingFrequency);

    for (size_t k = 0; k < N; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (size_t n = 0; n < N; ++n) {
            double theta = 2.0 * M_PI * k * n / N;
            sum += at(n) * std::exp(std::complex<double>(0.0, -theta));
        }
        spectrum[k] = sum;
    }

    return spectrum;
}

Signal Signal::IDFT(const std::vector<std::complex<double>> &spectrum) const
{
    const size_t N = spectrum.size();
    Signal reconstructedSignal(N, samplingFrequency);

    for (size_t n = 0; n < N; ++n) {
        std::complex<double> sum(0.0, 0.0);
        for (size_t k = 0; k < N; ++k) {
            double theta = 2.0 * M_PI * k * n / N;
            sum += spectrum[k] * std::exp(std::complex<double>(0.0, theta));
        }
        reconstructedSignal[n] = sum / static_cast<double>(N);
    }

    return reconstructedSignal;
}

/* ------------------------------ */

void Signal::display() const {
    for (const auto& element : *this) {
        std::cout << element << " ";
    }
    std::cout << std::endl;
}