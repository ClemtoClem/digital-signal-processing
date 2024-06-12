#ifndef __DEMODULATOR_HPP
#define __DEMODULATOR_HPP

#include <iostream>
#include <complex>
#include "Filter.hpp"

class Signal;

class Demodulator {
public:
	Demodulator(int freq_filter, int freq_oscillator, int sample_freq);
    void setup();
	void demodulate(Signal &signal, Signal &outputAmplitude, Signal &outputPhase);
private:
	int mFreqFilter, mFreqOscillator;
	int mSampleFreq;
    IIRFilter mFilter;
};

#endif // __DEMODULATOR_HPP