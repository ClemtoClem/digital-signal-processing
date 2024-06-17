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

std::vector<double> Signal::getRealBuffer() const {
    std::vector<double> realBuffer(this->size());
    for (size_t i=0; i<this->size(); i++)
        realBuffer[i] = (*this)[i].real();
    return realBuffer;
}

/* ------------------------------- */

Signal Signal::operator+(const Signal& other) const {
    Signal result(std::max(this->size(), other.size()), this->samplingFrequency);
    for (size_t i = 0; i < result.size(); i++) {
        std::complex<double> val1 = (i < this->size() ? (*this)[i] : 0);
        std::complex<double> val2 = (i < other.size() ? other[i] : 0);
        result[i] = val1 + val2;
    }
    return result;
}

Signal Signal::operator-(const Signal& other) const {
    Signal result(std::max(this->size(), other.size()), this->samplingFrequency);
    for (size_t i = 0; i < result.size(); i++) {
        std::complex<double> val1 = (i < this->size() ? (*this)[i] : 0);
        std::complex<double> val2 = (i < other.size() ? other[i] : 0);
        result[i] = val1 - val2;
    }
    return result;
}

Signal Signal::operator*(const Signal& other) const {
    Signal result(std::max(this->size(), other.size()), this->samplingFrequency);
    for (size_t i = 0; i < result.size(); i++) {
        std::complex<double> val1 = (i < this->size() ? (*this)[i] : 1);
        std::complex<double> val2 = (i < other.size() ? other[i] : 1);
        result[i] = val1 * val2;
    }
    return result;
}

Signal Signal::operator/(const Signal& other) const {
    Signal result(std::max(this->size(), other.size()), this->samplingFrequency);
    for (size_t i = 0; i < result.size(); i++) {
        std::complex<double> val1 = (i < this->size() ? (*this)[i] : 1);
        std::complex<double> val2 = (i < other.size() && other[i] != std::complex<double>(0, 0) ? other[i] : 1);
        result[i] = val1 / val2;
    }
    return result;
}

Signal &Signal::operator+=(const Signal &other){
    for (size_t i = 0; i < std::max(this->size(), other.size()); i++) {
        if (i >= this->size()) {
            this->push_back(other[i]);
        } else {
            (*this)[i] += (i < other.size() ? other[i] : 0);
        }
    }
    return *this;
}

Signal& Signal::operator-=(const Signal& other) {
    for (size_t i = 0; i < std::max(this->size(), other.size()); i++) {
        if (i >= this->size()) {
            this->push_back(-other[i]);
        } else {
            (*this)[i] -= (i < other.size() ? other[i] : 0);
        }
    }
    return *this;
}

Signal& Signal::operator*=(const Signal& other) {
    for (size_t i = 0; i < std::max(this->size(), other.size()); i++) {
        if (i >= this->size()) {
            this->push_back(0);
        } else {
            (*this)[i] *= (i < other.size() ? other[i] : 1);
        }
    }
    return *this;
}

Signal& Signal::operator/=(const Signal& other) {
    for (size_t i = 0; i < std::max(this->size(), other.size()); i++) {
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
    for (size_t i = 0; i < this->size(); i++) {
        result[i] = (*this)[i] + value;
    }
    return result;
}

Signal Signal::operator-(const std::complex<double>& value) const {
    Signal result(this->size(), this->samplingFrequency);
    for (size_t i = 0; i < this->size(); i++) {
        result[i] = (*this)[i] - value;
    }
    return result;
}

Signal Signal::operator*(const std::complex<double>& value) const {
    Signal result(this->size(), this->samplingFrequency);
    for (size_t i = 0; i < this->size(); i++) {
        result[i] = (*this)[i] * value;
    }
    return result;
}

Signal Signal::operator/(const std::complex<double>& value) const {
    Signal result(this->size(), this->samplingFrequency);
    for (size_t i = 0; i < this->size(); i++) {
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
    std::complex<double> maxElem = (*this)[0];
    for (size_t i = 0; i<this->size(); i++) {
        if (maxElem.real() < (*this)[i].real()) maxElem = (*this)[i];
    }
    return maxElem;
}

std::complex<double> Signal::min() const {
    if (this->empty()) {
        throw std::runtime_error("Signal is empty");
    }
    std::complex<double> minElem = (*this)[0];
    for (size_t i = 0; i<this->size(); i++) {
        if (minElem.real() > (*this)[i].real()) minElem = (*this)[i];
    }
    return minElem;
}

std::complex<double> Signal::mean() const {
    if (this->empty()) {
        throw std::runtime_error("Signal is empty");
    }
    std::complex<double> sum = 0;
    for (size_t i = 0; i<this->size(); i++) {
        sum += (*this)[i];
    }
    return sum / static_cast<double>(this->size());
}

/* ------------------------------- */

Signal Signal::DFT(size_t size_zero_padding) const {
    size_t N = size();
    size_t P = N + size_zero_padding;
    Signal spectrum(P, samplingFrequency);

    // Bit-reversal permutation
    size_t n = P;
    size_t bits = 0;
    while (n >>= 1) ++bits;

    vector<size_t> reversed(P);
    for (size_t i = 0; i < P; i++) {
        size_t j = 0;
        for (size_t k = 0; k < bits; ++k) {
            if (i & (1 << k)) j |= (1 << (bits - 1 - k));
        }
        reversed[i] = j;
    }

    // Copy input data to spectrum with bit-reversed order and zero-padding
    size_t i;
    for (i = 0; i < N; i++) {
        spectrum[reversed[i]] = (*this)[i];
    }
    for (i = N; i < P; i++) {
        spectrum[reversed[i]] = 0.0; // Zero-padding
    }

    // FFT algorithm
    for (size_t s = 1; s <= bits; ++s) {
        size_t m = 1 << s;
        std::complex<double> wm = std::exp(std::complex<double>(0, -2.0 * M_PI / m));
        for (size_t k = 0; k < P; k += m) {
            std::complex<double> w = 1;
            for (size_t j = 0; j < m / 2; ++j) {
                std::complex<double> t = w * spectrum[k + j + m / 2];
                std::complex<double> u = spectrum[k + j];
                spectrum[k + j] = u + t;
                spectrum[k + j + m / 2] = u - t;
                w *= wm;
            }
        }
    }

    for (i = 0; i < P; i++) {
        spectrum[i] = (spectrum[i]*2.0)/static_cast<double>(N);
    }

    return spectrum;
}

Signal Signal::IDFT(size_t size_zero_padding) const {
    size_t N = this->size();
    size_t P = N + size_zero_padding;
    Signal reconstructedSignal(P, samplingFrequency);

    // Bit-reversal permutation
    size_t n = P;
    size_t bits = 0;
    while (n >>= 1) bits++;

    vector<size_t> reversed(P);
    for (size_t i = 0; i < P; i++) {
        size_t j = 0;
        for (size_t k = 0; k < bits; ++k) {
            if (i & (1 << k)) j |= (1 << (bits - 1 - k));
        }
        reversed[i] = j;
    }

    // Copy input data to reconstructedSignal with bit-reversed order
    size_t i;
    for (i = 0; i < N; i++) {
        reconstructedSignal[reversed[i]] = (*this)[i];
    }
    for (i = N; i < P; i++) {
        reconstructedSignal[reversed[i]] = 0.0;
    }

    // IFFT algorithm
    for (size_t s = 1; s <= bits; s++) {
        size_t m = 1 << s;
        std::complex<double> wm = std::exp(std::complex<double>(0.0, 2.0 * M_PI / m)); // Note the sign change
        for (size_t k = 0; k < P; k += m) {
            std::complex<double> w = 1;
            for (size_t j = 0; j < m / 2; j++) {
                std::complex<double> t = w * reconstructedSignal[k + j + m / 2];
                std::complex<double> u = reconstructedSignal[k + j];
                reconstructedSignal[k + j] = u + t;
                reconstructedSignal[k + j + m / 2] = u - t;
                w *= wm;
            }
        }
    }

    // Normalize by dividing by N
    for (size_t i = 0; i < P; i++) {
        reconstructedSignal[i] /= N;
    }

    return reconstructedSignal;
}

void Signal::generateFrequencyAxis() {
    const size_t N = size();
    for (size_t k = 0; k < N; ++k) {
        (*this)[k] = k * samplingFrequency / N;
    }
}

/* ------------------------------ */

void Signal::display() const {
    for (const auto& element : *this) {
        std::cout << element << " ";
    }
    std::cout << std::endl;
}