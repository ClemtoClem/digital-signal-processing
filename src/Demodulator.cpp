#include "Demodulator.hpp"
#include "Signal.hpp"

Demodulator::Demodulator(double freq_filter, double freq_oscillator):
    _freqFilter(freq_filter), _freqOscillator(freq_oscillator), _filter() {
    _filter.set(4, _freqFilter, 0, FilterGabarit::LOW_PASS, AnalogFilter::BUTTERWORTH);
}

void Demodulator::setup() {
    _filter.setup();
}

void Demodulator::demodulate(Signal &signal, Signal &outputAmplitude, Signal &outputPhase, bool rms) {
    Signal sinus(signal.size());
    Signal cosinus(signal.size());
    sinus.generateWaveform(WaveformType::SINUS, 1.0, _freqOscillator);
    cosinus.generateWaveform(WaveformType::COSINUS, 1.0, _freqOscillator);

    Signal tempA = signal * sinus;
    tempA = _filter.apply(tempA);

    Signal tempPhi = signal * cosinus;
    tempPhi = _filter.apply(tempPhi);

    outputAmplitude.resize(signal.size());
    outputPhase.resize(signal.size());
    for (size_t i = 0; i<tempPhi.size(); i++) {
        if (rms) {
            outputAmplitude[i] = std::sqrt(tempA[i]*tempA[i] + tempPhi[i]*tempPhi[i])*sqrt(2.0);
        } else {
            outputAmplitude[i] = std::sqrt(2.0*tempA[i]*tempA[i] + 2.0*tempPhi[i]*tempPhi[i])*sqrt(2.0);
        }
        outputPhase[i] = (tempA[i] != 0.0)? atan2(tempPhi[i], tempA[i]) : 0.0;
    }
}
