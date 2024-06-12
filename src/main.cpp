#include <iostream>
#include "Signal.hpp"
#include "CSVFile.hpp"
#include "Filter.hpp"
#include "Demodulator.hpp"

const int SAMPLE_FREQUENCY = 250e6/16;
const int SIGNAL_SIZE = 16384;

void test_demodulate() {
    try {
        Signal s1(SIGNAL_SIZE, SAMPLE_FREQUENCY);
        Signal s2(SIGNAL_SIZE, SAMPLE_FREQUENCY);

        Signal signal_output        (SIGNAL_SIZE, SAMPLE_FREQUENCY, "output(t)");
        //Signal signal_input       (SIGNAL_SIZE, SAMPLE_FREQUENCY, "input(t)");
        Signal signal_demAmpli      (SIGNAL_SIZE, SAMPLE_FREQUENCY, "amplitude(t)");
        Signal signal_demPhase      (SIGNAL_SIZE, SAMPLE_FREQUENCY, "phase(t)");
        Signal signalLP             (SIGNAL_SIZE, SAMPLE_FREQUENCY, "signalFiltré(t)");

        Signal DFT_signal_output    (SIGNAL_SIZE, SAMPLE_FREQUENCY, "DFT_output(t)");
        //Signal DFT_signal_input   (SIGNAL_SIZE, SAMPLE_FREQUENCY, "DFT_input(t)");
        Signal DFT_signal_demAmpli  (SIGNAL_SIZE, SAMPLE_FREQUENCY, "DFT_amplitude(t)");
        Signal DFT_signal_demPhase  (SIGNAL_SIZE, SAMPLE_FREQUENCY, "DFT_phase(t)");

        std::cout << "Generate waveform" << std::endl;
        
        double s1_amplitude = 0.8;
        double s1_frequency = 10e3;
        s1.generateWaveform(WaveformType::SINUS, s1_amplitude, s1_frequency);
        
        double s2_amplitude = 0.2;
        double s2_frequency = 100e3;
        s2.generateWaveform(WaveformType::SINUS, s2_amplitude, s2_frequency);
        
        signal_output = s1+s2;
        //DFT_signal_output = signal_output.DFT();

        std::cout << "Filtrage" << std::endl;

        IIRFilter iir_filter;

        /* Filtre pass bas */
        iir_filter.set(4, 20e3, 0, SAMPLE_FREQUENCY, FilterGabarit::LOW_PASS, AnalogFilter::BUTTERWORTH);
        iir_filter.setup();
        std::cout << "Filtre passe bas :\n";
        iir_filter.printCoefficients();
        signalLP = iir_filter.filtering(signal_output);

        double filter_frequency = 3e3;
        double oscillator_frequency = 10e3;
        Demodulator dem(filter_frequency, oscillator_frequency, SAMPLE_FREQUENCY);
        dem.setup();
        dem.demodulate(signal_output, signal_demAmpli, signal_demPhase);

        DFT_signal_output   = signal_output.DFT();
        DFT_signal_demAmpli = signal_demAmpli.DFT();
        DFT_signal_demPhase = signal_demPhase.DFT();

        std::cout << "Sauvgarde des signaux" << std::endl;

        CSVFile outFile("./data/test.csv");
        std::vector<Signal> outSignals;
        outSignals.emplace_back(signal_output);
        outSignals.emplace_back(signalLP);
        outSignals.emplace_back(signal_demAmpli);
        outSignals.emplace_back(signal_demPhase);
        outFile.writeSignals(outSignals, false); // with time axis

        CSVFile outFileDFT("./data/test_DFT.csv");
        std::vector<Signal> outSignalsDFT;
        outSignalsDFT.emplace_back(DFT_signal_output);
        outSignalsDFT.emplace_back(DFT_signal_demAmpli);
        outSignalsDFT.emplace_back(DFT_signal_demPhase);
        outFileDFT.writeSignals(outSignalsDFT, true); // with frequency axis

        std::cout << "Test success finish\n";
    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

int main() {
    test_demodulate();

    return 0;
}