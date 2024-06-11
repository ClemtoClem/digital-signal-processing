#ifndef __SIGNAL_HPP
#define __SIGNAL_HPP

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <stdexcept>
#include <functional>

class BaseFilter;

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
    void resize(size_t new_size);
    void clear();
    int getSampleFrequency() const;
    void setDefaultValue(const complexd &default_value);
    const complexd &getDefaulValue() const;
    const std::vector<complexd> &getBuffer() const;
    std::vector<complexd> &getBuffer();

    Signal operator - ();

    Signal operator + (const Signal &other);
    Signal operator - (const Signal &other);
    Signal operator * (const Signal &other);
    Signal operator / (const Signal &other);

    Signal operator + (const complexd &c);
    friend Signal operator + (const complexd &c, const Signal &sig);
    Signal operator - (const complexd &c);
    friend Signal operator - (const complexd &c, const Signal &sig);
    Signal operator * (const complexd &c);
    friend Signal operator * (const complexd &c, const Signal &sig);
    Signal operator / (const complexd &c);
    friend Signal operator / (const complexd &c, const Signal &sig);

    Signal &operator += (const Signal &other);
    Signal &operator -= (const Signal &other);
    Signal &operator *= (const Signal &other);
    Signal &operator /= (const Signal &other);

    Signal &operator += (const complexd &c);
    Signal &operator -= (const complexd &c);
    Signal &operator *= (const complexd &c);
    Signal &operator /= (const complexd &c);

    Signal &self_square();
    Signal square() const;
    friend Signal square(const Signal &sig) { return sig.square(); }
    
    Signal &self_sqrt();
    Signal sqrt() const;
    friend Signal sqrt(const Signal &sig) { return sig.sqrt(); }
    
    Signal &self_tan();
    Signal tan() const;
    friend Signal tan(const Signal &sig) { return sig.tan(); }

    Signal &self_atan();
    Signal atan() const;
    friend Signal atan(const Signal &sig) { return sig.atan(); }

    Signal &self_cos();
    Signal cos() const;
    friend Signal cos(const Signal &sig) { return sig.cos(); }

    Signal &self_sin();
    Signal sin() const;
    friend Signal sin(const Signal &sig) { return sig.sin(); }

    Signal &self_abs();
    Signal abs() const;
    friend Signal abs(const Signal &sig) { return sig.abs(); }

    Signal normalize() const;
    friend Signal normalize(const Signal &sig) { return sig.normalize(); }

    complexd max() const;
    friend complexd max(const Signal &sig) { return sig.max(); }

    complexd min() const;
    friend complexd min(const Signal &sig) { return sig.min(); }

    void ceil(int precision);


    /* Génération des formes d'onde */
    Signal &setWaveform(const Waveform &wf);
    Signal &setWaveform(WaveformType type, double amplitude, double frequency,
        double phase = 0.0, double offset = 0.0, double delay = 0.0, double duty_cycle = 50.0);
    Waveform &getWaveform();
    void generate();
    // générer la forma du signal à partir fonction qui calcul la valeur d'un echantillon par rapport à sa position en temps
    void generateAbitraryForm(std::function<complexd(double)> &equation);

    /* Effectuer une fft sur le signal */
    Signal fft();
    /* Effectuer une fft inverse sur le signal */
    Signal ifft();

    Signal filter(BaseFilter &filter);

private:
    int m_sample_frequency;
    complexd m_defaul_value;
    std::vector<complexd> m_buffer;
    Waveform m_waveform;
};

#endif // __SIGNAL_HPP