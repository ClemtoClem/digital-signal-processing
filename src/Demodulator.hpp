#ifndef __DEMODULATOR_HPP
#define __DEMODULATOR_HPP

#include <iostream>
#include <complex>
#include "Filter.hpp"

class Signal;

class Demodulator {
public:
	Demodulator(double freq_filter, double freq_oscillator);
    void setup();
	void demodulate(Signal &signal, Signal &outputAmplitude, Signal &outputPhase);
private:
	double mFreqFilter, mFreqOscillator;
    IIRFilter mFilter;
};

#endif // __DEMODULATOR_HPP