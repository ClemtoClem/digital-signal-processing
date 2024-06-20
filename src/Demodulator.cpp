#include "Demodulator.hpp"
#include "Signal.hpp"

Demodulator::Demodulator(double freq_filter, double freq_oscillator):
    mFreqFilter(freq_filter), mFreqOscillator(freq_oscillator), mFilter() {
    mFilter.set(4, mFreqFilter, 0, FilterGabarit::LOW_PASS, AnalogFilter::BUTTERWORTH);
}

void Demodulator::setup() {
    mFilter.setup();
    std::cout << "Filtre passe bas démodulation :\n";
    mFilter.printCoefficients();
}

void Demodulator::demodulate(Signal &signal, Signal &outputAmplitude, Signal &outputPhase) {
    Signal sinus(signal.size());
    Signal cosinus(signal.size());
    sinus.generateWaveform(WaveformType::SINUS, 1.0, mFreqOscillator);
    cosinus.generateWaveform(WaveformType::COSINUS, 1.0, mFreqOscillator);

    Signal tempA = signal * sinus;
    tempA = mFilter.process(tempA);

    Signal tempPhi = signal * cosinus;
    tempPhi = mFilter.process(tempPhi);

    outputAmplitude.resize(signal.size());
    outputPhase.resize(signal.size());
    for (size_t i = 0; i<tempPhi.size(); i++) {
        outputAmplitude[i] = std::sqrt(2.0*tempA[i]*tempA[i] + 2.0*tempPhi[i]*tempPhi[i]);
        outputPhase[i] = (tempA[i] != 0.0)? atan2(tempPhi[i], tempA[i]) : 0.0;
    }
}
