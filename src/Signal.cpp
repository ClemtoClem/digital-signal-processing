#include "Signal.hpp"
#include <sstream>

Signal::Signal(size_t size, int sample_frequency, complexd default_value)
    : m_sample_frequency(sample_frequency), m_defaul_value(default_value), m_buffer(size, default_value) {
        m_waveform.type = WaveformType::UNUSED;
    }

Signal::~Signal() {}

Signal::Signal(const Signal &other) :
    m_sample_frequency(other.m_sample_frequency),
    m_defaul_value(other.m_defaul_value),
    m_buffer(other.m_buffer) {
        m_waveform.type = WaveformType::UNUSED;
    }

Signal &Signal::operator=(const Signal &other) {
    if (this != &other) {
        m_sample_frequency = other.m_sample_frequency;
        m_defaul_value = other.m_defaul_value;
        m_buffer = other.m_buffer;
    }
    return *this;
}

complexd &Signal::operator[](size_t index) {
    if (index >= m_buffer.size()) {
        std::stringstream ss;
        ss << "Signal::operator[], Index out of range " << index << " >= " << m_buffer.size();
        throw std::out_of_range(ss.str());
    }
    return m_buffer[index];
}

const complexd &Signal::operator[](size_t index) const {
    if (index >= m_buffer.size()) {
        std::stringstream ss;
        ss << "Signal::operator[], Index out of range " << index << " >= " << m_buffer.size();
        throw std::out_of_range(ss.str());
    }
    return m_buffer[index];
}

size_t Signal::size() const {
    return m_buffer.size();
}

int Signal::getSampleFrequency() const {
    return m_sample_frequency;
}

const complexd &Signal::getDefaulValue() const {
    return m_defaul_value;
}

const std::vector<complexd> &Signal::getBuffer() const {
    return m_buffer;
}

std::vector<complexd> &Signal::getBuffer() {
    return m_buffer;
}

Signal Signal::operator-() {
    Signal result(m_buffer.size(), m_sample_frequency, m_defaul_value);
    for (size_t i = 0; i < m_buffer.size(); i++) {
        result.m_buffer[i] = -m_buffer[i];
    }
    return result;
}

Signal Signal::operator+(const Signal &other) {
    if (m_buffer.size() != other.m_buffer.size()) {
        throw std::invalid_argument("Signals must be of the same size");
    }
    Signal result(m_buffer.size(), m_sample_frequency, m_defaul_value);
    for (size_t i = 0; i < m_buffer.size(); i++) {
        result.m_buffer[i] = m_buffer[i] + other.m_buffer[i];
    }
    return result;
}

Signal Signal::operator-(const Signal &other) {
    if (m_buffer.size() != other.m_buffer.size()) {
        throw std::invalid_argument("Signals must be of the same size");
    }
    Signal result(m_buffer.size(), m_sample_frequency, m_defaul_value);
    for (size_t i = 0; i < m_buffer.size(); i++) {
        result.m_buffer[i] = m_buffer[i] - other.m_buffer[i];
    }
    return result;
}

Signal Signal::operator*(const Signal &other) {
    if (m_buffer.size() != other.m_buffer.size()) {
        throw std::invalid_argument("Signals must be of the same size");
    }
    Signal result(m_buffer.size(), m_sample_frequency, m_defaul_value);
    for (size_t i = 0; i < m_buffer.size(); i++) {
        result.m_buffer[i] = m_buffer[i] * other.m_buffer[i];
    }
    return result;
}

Signal Signal::operator/(const Signal &other) {
    if (m_buffer.size() != other.m_buffer.size()) {
        throw std::invalid_argument("Signals must be of the same size");
    }
    Signal result(m_buffer.size(), m_sample_frequency, m_defaul_value);
    for (size_t i = 0; i < m_buffer.size(); i++) {
        if (other.m_buffer[i] == complexd(0, 0)) {
            throw std::invalid_argument("Division by zero");
        }
        result.m_buffer[i] = m_buffer[i] / other.m_buffer[i];
    }
    return result;
}

Signal Signal::operator+(const complexd &c) {
    Signal result(m_buffer.size(), m_sample_frequency, m_defaul_value);
    for (size_t i = 0; i < m_buffer.size(); i++) {
        result.m_buffer[i] = m_buffer[i] + c;
    }
    return result;
}

Signal Signal::operator-(const complexd &c) {
    Signal result(m_buffer.size(), m_sample_frequency, m_defaul_value);
    for (size_t i = 0; i < m_buffer.size(); i++) {
        result.m_buffer[i] = m_buffer[i] - c;
    }
    return result;
}

Signal Signal::operator*(const complexd &c) {
    Signal result(m_buffer.size(), m_sample_frequency, m_defaul_value);
    for (size_t i = 0; i < m_buffer.size(); i++) {
        result.m_buffer[i] = m_buffer[i] * c;
    }
    return result;
}

Signal Signal::operator/(const complexd &c) {
    if (c == complexd(0, 0)) {
        throw std::invalid_argument("Division by zero");
    }
    Signal result(m_buffer.size(), m_sample_frequency, m_defaul_value);
    for (size_t i = 0; i < m_buffer.size(); i++) {
        result.m_buffer[i] = m_buffer[i] / c;
    }
    return result;
}

Signal &Signal::operator+=(const Signal &other) {
    if (m_buffer.size() != other.m_buffer.size()) {
        throw std::invalid_argument("Signals must be of the same size");
    }
    for (size_t i = 0; i < m_buffer.size(); i++) {
        m_buffer[i] += other.m_buffer[i];
    }
    return *this;
}

Signal &Signal::operator-=(const Signal &other) {
    if (m_buffer.size() != other.m_buffer.size()) {
        throw std::invalid_argument("Signals must be of the same size");
    }
    for (size_t i = 0; i < m_buffer.size(); i++) {
        m_buffer[i] -= other.m_buffer[i];
    }
    return *this;
}

Signal &Signal::operator*=(const Signal &other) {
    if (m_buffer.size() != other.m_buffer.size()) {
        throw std::invalid_argument("Signals must be of the same size");
    }
    for (size_t i = 0; i < m_buffer.size(); i++) {
        m_buffer[i] *= other.m_buffer[i];
    }
    return *this;
}

Signal &Signal::operator/=(const Signal &other) {
    if (m_buffer.size() != other.m_buffer.size()) {
        throw std::invalid_argument("Signals must be of the same size");
    }
    for (size_t i = 0; i < m_buffer.size(); i++) {
        if (other.m_buffer[i] == complexd(0, 0)) {
            throw std::invalid_argument("Division by zero");
        }
        m_buffer[i] /= other.m_buffer[i];
    }
    return *this;
}

Signal &Signal::operator+=(const complexd &c) {
    for (size_t i = 0; i < m_buffer.size(); i++) {
        m_buffer[i] += c;
    }
    return *this;
}

Signal &Signal::operator-=(const complexd &c) {
    for (size_t i = 0; i < m_buffer.size(); i++) {
        m_buffer[i] -= c;
    }
    return *this;
}

Signal &Signal::operator*=(const complexd &c) {
    for (size_t i = 0; i < m_buffer.size(); i++) {
        m_buffer[i] *= c;
    }
    return *this;
}

Signal &Signal::operator/=(const complexd &c) {
    if (c == complexd(0, 0)) {
        throw std::invalid_argument("Division by zero");
    }
    for (size_t i = 0; i < m_buffer.size(); i++) {
        m_buffer[i] /= c;
    }
    return *this;
}

Signal Signal::abs() {
    Signal result(m_buffer.size(), m_sample_frequency, m_defaul_value);
    for (size_t i = 0; i < m_buffer.size(); i++) {
        result.m_buffer[i] = std::abs(m_buffer[i]);
    }
    return result;
}

Signal Signal::normalize() {
    Signal result(m_buffer.size(), m_sample_frequency, m_defaul_value);
    complexd cmax = m_buffer[0];

    // Trouver la valeur maximale absolue dans le signal
    for (size_t i = 0; i < m_buffer.size(); i++) {
        if (std::abs(m_buffer[i]) > std::abs(cmax)) {
            cmax = m_buffer[i];
        }
    }

    // Si la valeur maximale est différente de zéro, normaliser
    if (std::abs(cmax) > 0) {
        for (size_t i = 0; i < m_buffer.size(); i++) {
            result.m_buffer[i] = m_buffer[i] / cmax;
        }
    }
    return result;
}

complexd Signal::max() {
    complexd cmax = m_buffer[0];

    // Trouver la valeur maximale absolue dans le signal
    for (size_t i = 0; i < m_buffer.size(); i++) {
        if (std::abs(m_buffer[i]) > std::abs(cmax)) {
            cmax = m_buffer[i];
        }
    }
    return cmax;
}

complexd Signal::min() {
    complexd cmin = m_buffer[0];

    // Trouver la valeur maximale absolue dans le signal
    for (size_t i = 0; i < m_buffer.size(); i++) {
        if (std::abs(m_buffer[i]) < std::abs(cmin)) {
            cmin = m_buffer[i];
        }
    }
    return cmin;
}

void Signal::ceil(int precision) {
    double factor = std::pow(10, precision);
    for (size_t i = 0; i < m_buffer.size(); i++) {
        double real_part = std::ceil(m_buffer[i].real() * factor) / factor;
        double imag_part = std::ceil(m_buffer[i].imag() * factor) / factor;
        m_buffer[i] = complexd(real_part, imag_part);
    }
}

void Signal::setWaveform(const Waveform &wf) {
    m_waveform = wf;
}

Waveform &Signal::getWaveform() {
    return m_waveform;
}

void Signal::generate() {
    double sample_period = 1.0 / m_sample_frequency;
    double delay_samples = m_waveform.delay / 1000.0 * m_sample_frequency;

    for (size_t i = 0; i < m_buffer.size(); i++) {
        double t = (i - delay_samples) * sample_period;
        double value = 0;

        switch (m_waveform.type) {
            case WaveformType::SINUS:
                value = m_waveform.amplitude * std::sin(2 * M_PI * m_waveform.frequency * t + m_waveform.phase) + m_waveform.offset;
                break;
            case WaveformType::COSINUS:
                value = m_waveform.amplitude * std::cos(2 * M_PI * m_waveform.frequency * t + m_waveform.phase) + m_waveform.offset;
                break;
            case WaveformType::SQUARE:
                if (t >= 0) {
                    double period = 1.0 / m_waveform.frequency;
                    double time_in_period = fmod(t, period);
                    value = (time_in_period < m_waveform.duty_cycle * period) ? m_waveform.amplitude : -m_waveform.amplitude;
                    value += m_waveform.offset;
                }
                break;
            case WaveformType::TRIANGLE:
                if (t >= 0) {
                    double period = 1.0 / m_waveform.frequency;
                    double time_in_period = fmod(t, period);
                    double ramp_up_duration = m_waveform.duty_cycle * period;
                    if (time_in_period < ramp_up_duration) {
                        value = (2 * m_waveform.amplitude / ramp_up_duration) * time_in_period - m_waveform.amplitude;
                    } else {
                        value = (-2 * m_waveform.amplitude / (period - ramp_up_duration)) * (time_in_period - ramp_up_duration) + m_waveform.amplitude;
                    }
                    value += m_waveform.offset;
                }
                break;
            case WaveformType::RAMP_UP:
                if (t >= 0) {
                    double period = 1.0 / m_waveform.frequency;
                    double time_in_period = fmod(t, period);
                    value = (2 * m_waveform.amplitude / period) * time_in_period - m_waveform.amplitude + m_waveform.offset;
                }
                break;
            case WaveformType::RAMP_DOWN:
                if (t >= 0) {
                    double period = 1.0 / m_waveform.frequency;
                    double time_in_period = fmod(t, period);
                    value = (-2 * m_waveform.amplitude / period) * time_in_period + m_waveform.amplitude + m_waveform.offset;
                }
                break;
            case WaveformType::POSITIVE_DC:
                value = m_waveform.amplitude + m_waveform.offset;
                break;
            case WaveformType::NEGATIVE_DC:
                value = -m_waveform.amplitude + m_waveform.offset;
                break;
            case WaveformType::UNUSED:
            default:
                value = m_defaul_value.real();
                break;
        }

        m_buffer[i] = complexd(value, 0);
    }
}

Signal Signal::demodulate(int local_oscillator_freq, bool sin_or_cos)
{
    Signal demodulate_signal(m_buffer.size(), m_sample_frequency, m_defaul_value);
    for (size_t i = 0; i < m_buffer.size(); i++) {
        double t = i * m_sample_frequency;
        if (sin_or_cos == true)
            demodulate_signal[i] = std::cos(2 * M_PI * local_oscillator_freq * t) * m_buffer[i];
        else
            demodulate_signal[i] = std::sin(2 * M_PI * local_oscillator_freq * t) * m_buffer[i];
    }
    return demodulate_signal;
}

// Function to perform FFT recursively
static void __fft(complexd *x, complexd *X, size_t N, bool inverse) {
    if (N <= 1) {
        X[0] = x[0];
        return;
    }

    // Split the sequence into even and odd parts
    complexd *even = new complexd[N / 2];
    complexd *odd = new complexd[N / 2];
    for (size_t i = 0; i < N / 2; ++i) {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    // Recursive calls for FFT on even and odd parts
    complexd *Even = new complexd[N / 2];
    complexd *Odd = new complexd[N / 2];
    __fft(even, Even, N / 2, inverse);
    __fft(odd, Odd, N / 2, inverse);

    // Combine results
    double angle = (inverse ? 2 : -2) * M_PI / N;
    complexd w(1);
    complexd wn(std::cos(angle), std::sin(angle));
    for (size_t k = 0; k < N / 2; ++k) {
        X[k] = Even[k] + w * Odd[k];
        X[k + N / 2] = Even[k] - w * Odd[k];
        if (inverse) {
            X[k] /= 2;
            X[k + N / 2] /= 2;
        }
        w *= wn;
    }

    delete[] even;
    delete[] odd;
    delete[] Even;
    delete[] Odd;
}

Signal Signal::fft() {
    Signal result(m_buffer.size(), m_sample_frequency, m_defaul_value);
    __fft(m_buffer.data(), result.m_buffer.data(), m_buffer.size(), false);
    return (result*2*M_PI)/m_buffer.size();
}

Signal Signal::ifft() {
    Signal result(m_buffer.size(), m_sample_frequency, m_defaul_value);
    __fft(m_buffer.data(), result.m_buffer.data(), m_buffer.size(), true);
    return result;
}