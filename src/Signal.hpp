#ifndef __SIGNAL_HPP
#define __SIGNAL_HPP

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <stdexcept>

using complexd = std::complex<double>;

enum class WaveformType {
    UNUSED,
    SINUS,
    COSINUS,
    SQUARE,   // with duty cycle
    TRIANGLE, // with duty cycle
    RAMP_UP,
    RAMP_DOWN,
    POSITIVE_DC,
    NEGATIVE_DC,
};

struct Waveform {
    WaveformType type;
    double amplitude;
    double frequency;
    double phase;
    double offset;
    double delay;  // in ms
    double duty_cycle;
};

class Signal {
public:
    Signal(size_t size, int sample_frequency, complexd default_value = 0);
    Signal(const Signal &other);
    ~Signal();

    Signal &operator = (const Signal &other);

    complexd &operator [] (size_t index);
    const complexd &operator [] (size_t index) const;
    size_t size() const;
    int getSampleFrequency() const;
    const complexd &getDefaulValue() const;
    const std::vector<complexd> &getBuffer() const;
    std::vector<complexd> &getBuffer();

    Signal operator - ();

    Signal operator + (const Signal &other);
    Signal operator - (const Signal &other);
    Signal operator * (const Signal &other);
    Signal operator / (const Signal &other);

    Signal operator + (const complexd &c);
    Signal operator - (const complexd &c);
    Signal operator * (const complexd &c);
    Signal operator / (const complexd &c);

    Signal &operator += (const Signal &other);
    Signal &operator -= (const Signal &other);
    Signal &operator *= (const Signal &other);
    Signal &operator /= (const Signal &other);

    Signal &operator += (const complexd &c);
    Signal &operator -= (const complexd &c);
    Signal &operator *= (const complexd &c);
    Signal &operator /= (const complexd &c);

    Signal abs();
    Signal normalize();
    complexd max();
    complexd min();
    void ceil(int precision);

    /* Génération des formes d'onde */
    void setWaveform(const Waveform &wf);
    Waveform &getWaveform();
    void generate();

    /* Fonction de démodulation avec un oscillateur local */
    Signal demodulate(int local_oscillator_freq, bool sin_or_cos = false);

    /* Effectuer une fft sur le signal */
    Signal fft();
    /* Effectuer une fft inverse sur le signal */
    Signal ifft();

private:
    int m_sample_frequency;
    complexd m_defaul_value;
    std::vector<complexd> m_buffer;
    Waveform m_waveform;
};

#endif // __SIGNAL_HPP