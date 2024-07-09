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

	/**
	 * @brief Demodulate a signal
	 * @param signal Input signal
	 * @param outputAmplitude Output signal containing amplitude of the signal
	 * @param outputPhase Output signal containing phase of the signal
	 * @param rms If true, output amplitude will be the RMS value of the signal
	 */
	void demodulate(Signal &signal, Signal &outputAmplitude, Signal &outputPhase, bool rms = false);
private:
	double _freqFilter, _freqOscillator;
    IIRFilter _filter;
};

#endif // __DEMODULATOR_HPP