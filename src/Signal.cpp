#include "Signal.hpp"
#include "Spectrum.hpp"

Signal::Signal(const std::string &name) : std::vector<double>(BUFFER_SIZE, 0.0), mName(name) {}

Signal::Signal(size_t size, const std::string &name) : std::vector<double>(size), mName(name) {}

Signal::Signal(const std::vector<double> &values, const std::string &name) : std::vector<double>(values), mName(name) {}

Signal::Signal(const Signal &other) : std::vector<double>(other), mName(other.mName) {}

Signal &Signal::operator=(const Signal &other) {
    if (this != &other) {
        std::vector<double>::operator=(other);
    }
    return *this;
}

/* ------------------------------- */

const std::string &Signal::getName() const {
    return mName;
}

void Signal::setName(const std::string &name) {
    mName = name;
}

bool Signal::hasInfitityValue() const {
    for (size_t i = 0; i < size(); i++) {
        if (((*this)[i] == INFINITY) || ((*this)[i] == -INFINITY)) {
            return true;
        }
    }
    return false;
}

/* ------------------------------- */

Signal &Signal::addValue(double value) {
    std::vector<double>::push_back(value);
    return *this;
}

Signal &Signal::fill(double value)
{
    for (size_t i = 0; i < size(); i++) {
        (*this)[i] = value;
    }
    return *this;
}

/* ------------------------------- */

Signal Signal::operator+(const Signal &other) const {
    Signal output(other.size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = (*this)[i] + other[i];
    }
    return output;
}

Signal Signal::operator-(const Signal &other) const {
    Signal output(other.size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = (*this)[i] - other[i];
    }
    return output;
}

Signal Signal::operator*(const Signal &other) const {
    Signal output(other.size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = (*this)[i] * other[i];
    }
    return output;
}

Signal Signal::operator/(const Signal &other) const {
    Signal output(other.size());
    for (size_t i = 0; i < size(); i++) {
        if (other[i] == 0) output[i] = INFINITY;
        else output[i] = (*this)[i] / other[i];
    }
    return output;
}

Signal &Signal::operator+=(const Signal &other) {
    for (size_t i = 0; i < size(); i++) {
        (*this)[i] += other[i];
    }
    return *this;
}

Signal &Signal::operator-=(const Signal &other) {
    for (size_t i = 0; i < size(); i++) {
        (*this)[i] -= other[i];
    }
    return *this;
}

Signal &Signal::operator*=(const Signal &other) {
    for (size_t i = 0; i < size(); i++) {
        (*this)[i] *= other[i];
    }
    return *this;
}

Signal &Signal::operator/=(const Signal &other) {
    for (size_t i = 0; i < size(); i++) {
        if (other[i] == 0) (*this)[i] = INFINITY;
        else (*this)[i] /= other[i];
    }
    return *this;
}

/* ------------------------------- */

Signal Signal::operator+(double value) const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = (*this)[i] + value;
    }
    return output;
}

Signal Signal::operator-(double value) const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = (*this)[i] - value;
    }
    return output;
}

Signal Signal::operator*(double value) const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = (*this)[i] * value;
    }
    return output;
}

Signal Signal::operator/(double value) const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        if (value == 0) output[i] = INFINITY;
        else output[i] = (*this)[i] / value;
    }
    return output;
}

Signal &Signal::operator+=(double value) {
    for (size_t i = 0; i < size(); i++) {
        (*this)[i] += value;
    }
    return *this;
}

Signal &Signal::operator-=(double value) {
    for (size_t i = 0; i < size(); i++) {
        (*this)[i] -= value;
    }
    return *this;
}

Signal &Signal::operator*=(double value) {
    for (size_t i = 0; i < size(); i++) {
        (*this)[i] *= value;
    }
    return *this;
}

Signal &Signal::operator/=(double value) {
    for (size_t i = 0; i < size(); i++) {
        if (value == 0) (*this)[i] = INFINITY;
        else (*this)[i] /= value;
    }
    return *this;
}

/* ------------------------------- */

Signal operator+(double value, const Signal &input) {
    Signal output(input.size());
    for (size_t i = 0; i < input.size(); i++) {
        output[i] = input[i] + value;
    }
    return output;
}

Signal operator-(double value, const Signal &input) {
    Signal output(input.size());
    for (size_t i = 0; i < input.size(); i++) {
        output[i] = input[i] - value;
    }
    return output;
}

Signal operator*(double value, const Signal &input) {
    Signal output(input.size());
    for (size_t i = 0; i < input.size(); i++) {
        output[i] = input[i] * value;
    }
    return output;
}

Signal operator/(double value, const Signal &input) {
    Signal output(input.size());
    for (size_t i = 0; i < input.size(); i++) {
        if (value == 0) output[i] = INFINITY;
        else output[i] = input[i] / value;
    }
    return output;
}

/* ------------------------------- */

Signal Signal::cos() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::cos((*this)[i]);
    }
    return output;
}

Signal Signal::sin() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::sin((*this)[i]);
    }
    return output;
}

Signal Signal::tan() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::tan((*this)[i]);
    }
    return output;
}

Signal Signal::cosh() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::cosh((*this)[i]);
    }
    return output;
}

Signal Signal::sinh() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::sinh((*this)[i]);
    }
    return output;
}

Signal Signal::tanh() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::tanh((*this)[i]);
    }
    return output;
}

Signal Signal::acos() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::acos((*this)[i]);
    }
    return output;
}

Signal Signal::asin() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::asin((*this)[i]);
    }
    return output;
}

Signal Signal::atan() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::atan((*this)[i]);
    }
    return output;
}

Signal Signal::acosh() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::acosh((*this)[i]);
    }
    return output;
}

Signal Signal::asinh() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::asinh((*this)[i]);
    }
    return output;
}

Signal Signal::atanh() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::atanh((*this)[i]);
    }
    return output;
}

/* ------------------------------- */

Signal Signal::square() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = (*this)[i] * (*this)[i];
    }
    return output;
}

Signal Signal::pow(double exponent) const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::pow((*this)[i], exponent);
    }
    return output;
}

Signal Signal::sqrt() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::sqrt((*this)[i]);
    }
    return output;
}

Signal Signal::log() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::log((*this)[i]);
    }
    return output;
}

Signal Signal::log2() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::log2((*this)[i]);
    }
    return output;
}

Signal Signal::log10() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::log10((*this)[i]);
    }
    return output;
}

Signal Signal::exp() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::exp((*this)[i]);
    }
    return output;
}

/* ------------------------------- */

double Signal::max() const {
    if (this->empty()) {
        return NAN;
    }
    double maxElem = (*this)[0];
    for (size_t i = 0; i<this->size(); i++) {
        if (maxElem < (*this)[i]) maxElem = (*this)[i];
    }
    return maxElem;
}

double Signal::min() const {
    if (this->empty()) {
        return NAN;
    }
    double minElem = (*this)[0];
    for (size_t i = 0; i<this->size(); i++) {
        if (minElem > (*this)[i]) minElem = (*this)[i];
    }
    return minElem;
}

double Signal::mean() const {
    if (this->empty()) {
        return NAN;
    }
    double sum = 0;
    for (size_t i = 0; i<this->size(); i++) {
        sum += (*this)[i];
    }
    return sum / static_cast<double>(this->size());
}

/* ------------------------------- */

bool Signal::operator < (const Signal &input) const {
    for (size_t i = 0; i < size(); i++) {
        if ((*this)[i] >= input[i]) return false;
    }
    return true;
}

bool Signal::operator <= (const Signal &input) const {
    for (size_t i = 0; i < size(); i++) {
        if ((*this)[i] > input[i]) return false;
    }
    return true;
}

bool Signal::operator > (const Signal &input) const {
    for (size_t i = 0; i < size(); i++) {
        if ((*this)[i] <= input[i]) return false;
    }
    return true;
}

bool Signal::operator >= (const Signal &input) const {
    for (size_t i = 0; i < size(); i++) {
        if ((*this)[i] < input[i]) return false;
    }
    return true;
}

bool Signal::operator == (const Signal &input) const {
    for (size_t i = 0; i < size(); i++) {
        if ((*this)[i] != input[i]) return false;
    }
    return true;
}

bool Signal::operator != (const Signal &input) const {
    for (size_t i = 0; i < size(); i++) {
        if ((*this)[i] == input[i]) return false;
    }
    return true;
}

/* ------------------------------- */

void Signal::generateWaveform(WaveformType type, double amplitude, double frequency, double phase, double offset, double delay, double duty_cycle) {
    // Vérification du duty cycle pour les formes d'onde concernées
    if (duty_cycle < 0.0 || duty_cycle > 100.0) {
        std::cerr << "Error: Duty cycle must be between 0 and 100." << std::endl;
        return;
    }
    duty_cycle /= 100;

    // Calcul de la période d'échantillonnage
    double sample_period = 1.0 / SAMPLING_FREQUENCY;
    // Conversion du délai en nombre d'échantillons
    double delay_samples = delay / 1000.0 * SAMPLING_FREQUENCY;

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

        // Stockage de l'échantillon
        (*this)[i] = value;
    }
}


/* ------------------------------- */

Spectrum Signal::DFT(size_t size_zero_padding) const {
    size_t N = size();
    size_t P = N + size_zero_padding;
    std::vector<complexd> spectrum(P, 0.0);

    // Bit-reversal permutation
    size_t n = P;
    size_t bits = 0;
    while (n >>= 1) ++bits;

    std::vector<size_t> reversed(P);
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
        complexd wm = std::exp(complexd(0, -2.0 * M_PI / m));
        for (size_t k = 0; k < P; k += m) {
            complexd w = 1;
            for (size_t j = 0; j < m / 2; ++j) {
                complexd t = w * spectrum[k + j + m / 2];
                complexd u = spectrum[k + j];
                spectrum[k + j] = u + t;
                spectrum[k + j + m / 2] = u - t;
                w *= wm;
            }
        }
    }

    for (i = 0; i < P; i++) {
        spectrum[i] = (spectrum[i] * 2.0) / static_cast<double>(N);
    }

    return Spectrum(spectrum); // Return a Spectrum object
}

/* ------------------------------- */

double Signal::calculateNoiseRMS() const {
    double mean = this->mean();
    double sum_squares = 0.0;
    for (double value : *this) {
        double noise = value - mean;
        sum_squares += noise * noise;
    }
    return std::sqrt(sum_squares / this->size());
}

/* ------------------------------- */

std::ostream &operator<<(std::ostream &out, const Signal &signal) {
    out << signal.getName() << "{";
    for (size_t i = 0; i < signal.size(); i++) {
        out << signal[i];
        if (i < signal.size()-1) out << ",";
    }
    out << "}";
    return out;
}