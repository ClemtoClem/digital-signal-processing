#ifndef __GLOBALS_HPP
#define __GLOBALS_HPP

#include <cstdint>
#include <complex>

using complexd = std::complex<double>;

inline unsigned int DECIMATION = 16;
inline double MAX_SAMPLING_FREQUENCY = 125e6;
inline double SAMPLING_FREQUENCY = MAX_SAMPLING_FREQUENCY / DECIMATION;
inline size_t BUFFER_SIZE = 1 << 14;
inline size_t MAX_SIGNAL_SIZE = 1 << 14;

inline void SetMaxSampleFrequency(size_t freq) {
    if (MAX_SAMPLING_FREQUENCY <= 0.0) {
        std::cerr << "Sampling frequency must be greater than zero" << std::endl;
        return;
    }
    MAX_SAMPLING_FREQUENCY = freq;
    SAMPLING_FREQUENCY = MAX_SAMPLING_FREQUENCY/DECIMATION;
}

inline void SetDecimation(unsigned int decimation) {
    // vérifier que la décimation est une puissance de 2
    if (decimation != 0 && (decimation & (decimation - 1)) != 0) {
        std::cerr << "Decimation must be a power of 2" << std::endl;
        return;
    }
    DECIMATION = decimation;
    SAMPLING_FREQUENCY = MAX_SAMPLING_FREQUENCY/DECIMATION;
}

inline void SetBufferSize(size_t size) {
    if (size > MAX_SIGNAL_SIZE) BUFFER_SIZE = MAX_SIGNAL_SIZE;
    else BUFFER_SIZE = size;
}

#endif // __GLOBALS_HPP