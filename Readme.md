# Digital Signal Processing

## General Info

This project was developed with the aim of creating and testing a digital demodulator before implementing it on an embedded system.
This project, programmed in C++, offers the possibility of synthesizing real and complex signals and carrying out operations between them.

In addition, in order to carry out demodulation, an initial implementation of a digital impulse response filter has been put in place. This filter offers the possibility of choosing its template (high pass, low pass, band pass and band stop), its type allowing its coefficients to be calculated (Butterworth, Bessel, Chebyshev, elliptic) and its cutoff frequency.

For the moment only the Butterworth filter is implemented, which is sufficient to carry out amplitude and phase demodulation.
Lately, in order to visualize the signals the project allows complex signals to be exported in csv format and displayed in a python script with matplotlib and numpy.

## Technologies

This project was developed under Linux with the vscode editor.

The code was completed with version g++ 11.4.0.

The Python version used is 3.10.12.

## Example

```cpp
#include <iostream>
#include <cmath>
#include "CSVFile.hpp"
#include "Filter.hpp"
#include "Noise.hpp"

const int SAMPLING_FREQ = 250e6;
const uint32_t SIGNAL_SIZE = 16384;

const double pi = M_PI;

int main() {
    Signal s1(SIGNAL_SIZE, SAMPLING_FREQ, "s1(t)");
    Signal s2(SIGNAL_SIZE, SAMPLING_FREQ, "s2(t)");
    Signal out(SIGNAL_SIZE, SAMPLING_FREQ, "output(t)");
    Signal out_filtered(SIGNAL_SIZE, SAMPLING_FREQ, "output_filtered(t)");

    WaveformType s1_waveform = WaveformType::SINUS;
    double s1_amplitude = 0.8;
    double s1_frequency = 10e3;
    double s1_phase = pi/12;
    double s1_offset = 0.0;
    double s1_delay = 0.0;
    double s1_duty_cycle = 0.0;

    WaveformType s2_waveform = WaveformType::SQUARE;
    double s2_amplitude = 0.2;
    double s2_frequency = 100e3;
    double s2_phase = 0.0;
    double s2_offset = 0.0;
    double s2_delay = 0.0;
    double s2_duty_cycle = 25.0;

    s1.generateWaveform(s1_waveform::SINUS, s1_amplitude, s1_frequency, s1_phase, s1_offset, s1_delay, s1_duty_cycle);
    s2.generateWaveform(s2_waveform, s2_amplitude, s2_frequency, s2_phase, s2_offset, s2_delay, s2_duty_cycle);

    out = (s1 + s2) * 0.8;


    Noise noise;
    noise.setParams(NoiseType::WHITE, SAMPLING_FREQ, 0.05);
    out = noise.process(out);


    IIRFilter iir_filter;

    /* Low pass filter */
    iir_filter.set(4, 20e3, 0, SAMPLING_FREQ, FilterGabarit::LOW_PASS, AnalogFilter::BUTTERWORTH);
    iir_filter.setup();
    iir_filter.printCoefficients();
    out_filtered = iir_filter.process(out);

    CSVFile outFile("./data/test.csv");
    std::vector<Signal> outSignals;
    outSignals.emplace_back(s1);
    outSignals.emplace_back(s2);
    outSignals.emplace_back(out);
    outSignals.emplace_back(out_filtered);
    outFile.writeSignals(outSignals, false); // with time axis


    Signal DFT_s1(SIGNAL_SIZE, SAMPLING_FREQ, "S1(f)");
    Signal DFT_s2(SIGNAL_SIZE, SAMPLING_FREQ, "S2(f)");
    Signal DFT_out(SIGNAL_SIZE, SAMPLING_FREQ, "OUTPUT(f)");
    Signal DFT_out_filtered(SIGNAL_SIZE, SAMPLING_FREQ, "OUTPUT_FILTERED(f)");

    DFT_s1 = s1.DFT(SIGNAL_SIZE); // with zero padding
    DFT_s2 = s2.DFT(); // unless zero padding
    DFT_out = out.DFT(SIGNAL_SIZE);
    DFT_out_filtered = out_filtered.DFT(SIGNAL_SIZE);

    CSVFile outFileDFT("./data/test_DFT.csv");
    std::vector<Signal> outSignalsDFT;
    outSignalsDFT.emplace_back(DFT_s1);
    outSignalsDFT.emplace_back(DFT_s2);
    outSignalsDFT.emplace_back(DFT_out);
    outSignalsDFT.emplace_back(DFT_out_filtered);
    outFileDFT.writeSignals(outSignalsDFT, true); // with frequency axis

    return 0;
}
```