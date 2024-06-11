#include <iostream>
#include "Signal.hpp"
#include "CSV_File.hpp"
#include "Filter.hpp"

const int SAMPLE_FREQUENCY = 250e6/4;
const int SIGNAL_SIZE = 163840;

void initDemodulate(BaseFilter &filter, int freq_filter, int sample_freq) {
    filter.set(4, freq_filter, 0, sample_freq, FilterGabarit::LOW_PASS, AnalogFilter::BUTTERWORTH);
    filter.setup();
    std::cout << "Filtre passe bas démodulation :\n";
    filter.printCoefficients();
}

void demodulateSignal(Signal &signal, BaseFilter &filter, int freq_oscillator, Signal &outputAmplitude, Signal &outputPhase) {
    // Génération d'un signal sinusoïdale et cosinusoïdale et fréquence "freq_oscillator"
    Signal sinus(signal.size(), signal.getSampleFrequency());
    Signal cosinus(signal.size(), signal.getSampleFrequency());
    sinus.setWaveform(WaveformType::SINUS, 1.0, freq_oscillator).generate();
    cosinus.setWaveform(WaveformType::COSINUS, 1.0, freq_oscillator).generate();
    
    Signal tempA = signal * sinus;
    tempA = filter.filtering(tempA);

    Signal tempPhi = signal * cosinus;
    tempPhi = filter.filtering(tempPhi);

    outputAmplitude.resize(signal.size());
    outputPhase.resize(signal.size());
    for (size_t i = 0; i<tempPhi.size(); i++) {
        outputAmplitude[i] = sqrt(2.0*tempA[i]*tempA[i]+2.0*tempPhi[i]*tempPhi[i]);
        outputPhase[i] = (tempA[i] != 0.0)? atan(tempPhi[i]/tempA[i]) : 0.0;
    }
    
    /*
    double a, phi, a2, phi2;
    double To = 1/freq_oscillator;
    for (size_t i=0; i < signal.size(); i++) {
        // multiplier le signal par un sinus
        a = signal[i].real() * sin(To*i);
        // filtrer a
        a = filter1.eqdiff(a);

        // multiplier le signal par un cosinus
        phi = signal[i].real() * cos(To*f);
        // filtrer phi
        phi = filter2.eqdiff(phi);

        // calculer l'amplitude
        outputAmplitude[i] = sqrt(2*a*a + 2*phi*phi);
        // calculer la phase
        outputPhase[i] = (a2 != 0.0)? atan(phi/a) : 0.0;
    }*/

}



void test_signals() {

    Signal signal1(SIGNAL_SIZE, SAMPLE_FREQUENCY);
    Signal signal2(SIGNAL_SIZE, SAMPLE_FREQUENCY);

    signal1.setWaveform(WaveformType::SINUS, 0.8, 1000).generate();
    signal1.setWaveform(WaveformType::COSINUS, 0.8, 1000).generate();

    Signal signal3 = signal1+signal2;

    std::cerr << "signal3 = [";
    for (size_t i = 0; i<signal3.size(); i++) {
        std::cerr << signal3[i] << " ";
    }
    std::cerr << "]\n";
}

void test_filters() {
    try {
        CSV_File csv(SAMPLE_FREQUENCY);
        Signal signal1(SIGNAL_SIZE, SAMPLE_FREQUENCY); //auto signal1 = ccsv.create("s1", SIGNAL_SIZE, 0);
        Signal signal2(SIGNAL_SIZE, SAMPLE_FREQUENCY); //auto signal2 = csv.create("s2", SIGNAL_SIZE, 0);
        Signal signal3(SIGNAL_SIZE, SAMPLE_FREQUENCY); //auto signal3 = csv.create("s3", SIGNAL_SIZE, 0);
        auto signal4 = csv.create("s4=s1*s2*s3", SIGNAL_SIZE, 0);
        auto signalLP = csv.create("signalLP(s4)", SIGNAL_SIZE, 0);
        auto signalHP = csv.create("signalHP(s4)", SIGNAL_SIZE, 0);
        auto signalBP = csv.create("signalBP(s4)", SIGNAL_SIZE, 0);

        CSV_File fft_csv(SAMPLE_FREQUENCY);
        auto fft_s4 = fft_csv.create("fft(s4)", SIGNAL_SIZE, 0);
        auto fft_signalLP  = fft_csv.create("fft_signalLP(s4)", SIGNAL_SIZE, 0);
        auto fft_signalHP = fft_csv.create("fft_signalHP(s4)", SIGNAL_SIZE, 0);
        auto fft_signalBP = fft_csv.create("fft_signalBP(s4)", SIGNAL_SIZE, 0);

        Waveform wf_1 {
            .type = WaveformType::SINUS,
            .amplitude = 0.8,
            .frequency = 10e3,
            .phase = 0.0,
            .offset = 0.0,
            .delay = 0.0,
            .duty_cycle = 0.0,
        };
        
        Waveform wf_2 {
            .type = WaveformType::SINUS,
            .amplitude = 0.2,
            .frequency = 100e4,
            .phase = M_PI/8,
            .offset = 0.0,
            .delay = 0.0,
            .duty_cycle = 0.0
        };
        
        /*Waveform wf_3 {
            .type = WaveformType::SINUS,
            .amplitude = 0.1,
            .frequency = 100000,
            .phase = 0.0,
            .offset = 0.0,
            .delay = 0.0,
            .duty_cycle = 0.0
        };*/
        
        signal1.setWaveform(wf_1);
        signal2.setWaveform(wf_2);
        /*signal3.setWaveform(wf_3);*/

        signal1.generate();
        signal2.generate();
        signal3.generate();

        (*signal4) = signal1 + signal2 /*+ signal3*/;
        (*fft_s4) = signal4->fft();

        IIRFilter iir_filter;

        /* Filtre pass bas */
        iir_filter.set(4, 20e4, 0, SAMPLE_FREQUENCY, FilterGabarit::LOW_PASS, AnalogFilter::BUTTERWORTH);
        iir_filter.setup();
        std::cout << "Filtre passe bas :\n";
        iir_filter.printCoefficients();
        //(*signalLP) = signal4->filter(iir_filter);
        (*signalLP) = iir_filter.filtering((*signal4));
        (*fft_signalLP) = signalLP->fft();

        /* Filtre pass haut */
        iir_filter.set(4, 100e4, 0, SAMPLE_FREQUENCY, FilterGabarit::HIGH_PASS, AnalogFilter::BUTTERWORTH);
        iir_filter.setup();
        std::cout << "Filtre passe haut :\n";
        iir_filter.printCoefficients();
        //(*signalHP) = signal4->filter(iir_filter);
        (*signalHP) = iir_filter.filtering((*signal4));
        (*fft_signalHP) = signalHP->fft();

        /* Filtre pass bande */
        /*filter.set(4, 100, 10000, SAMPLE_FREQUENCY, FilterGabarit::BAND_PASS, AnalogFilter::BUTTERWORTH);
        filter.setup();
        std::cout << "Filtre coupe bande :\n";
        filter.printCoefficients();
        (*signalBP) = signal4->filter(filter);
        //(*signalBP) = filter.filter((*signal4));
        (*fft_signalBP) = signalBP->fft();*/


        csv.write("test.csv");
        fft_csv.write("test_fft.csv");

        std::cout << "Test success finish\n";
    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}


void test_demodulate() {
    try {
        Signal s1(SIGNAL_SIZE, SAMPLE_FREQUENCY);
        Signal s2(SIGNAL_SIZE, SAMPLE_FREQUENCY);

        CSV_File csv(SAMPLE_FREQUENCY);
        auto signal          = csv.create("signal(t)", SIGNAL_SIZE, 0);
        auto signal_demAmpli = csv.create("amplitude(t)", SIGNAL_SIZE, 0);
        auto signal_demPhase = csv.create("phase(t)", SIGNAL_SIZE, 0);
        auto tempA_after    = csv.create("tempA_after", SIGNAL_SIZE, 0);
        auto tempPhi_after  = csv.create("tempPhi_after", SIGNAL_SIZE, 0);
        auto tempA    = csv.create("tempA", SIGNAL_SIZE, 0);
        auto tempPhi  = csv.create("tempPhi", SIGNAL_SIZE, 0);
        auto signalLP = csv.create("signalLP", SIGNAL_SIZE, 0);

        CSV_File fft_csv(SAMPLE_FREQUENCY);
        auto fft_signal          = fft_csv.create("fft_signal(t)", SIGNAL_SIZE, 0);
        auto fft_signal_demAmpli = fft_csv.create("fft_amplitude(t)", SIGNAL_SIZE, 0);
        auto fft_signal_demPhase = fft_csv.create("fft_phase(t)", SIGNAL_SIZE, 0);

        Waveform wf_1 {
            .type = WaveformType::SINUS,
            .amplitude = 0.8,
            .frequency = 10e3,
            .phase = 0.0,
            .offset = 0.0,
            .delay = 0.0,
            .duty_cycle = 0.0,
        };
        s1.setWaveform(wf_1).generate();

        Waveform wf_2 {
            .type = WaveformType::SINUS,
            .amplitude = 0.2,
            .frequency = 100e3,
            .phase = 0.0,
            .offset = 0.0,
            .delay = 0.0,
            .duty_cycle = 0.0,
        };
        s2.setWaveform(wf_2).generate();
        
        (*signal) = s1+s2;
        (*fft_signal) = signal->fft();

        IIRFilter iir_filter;

        /* Filtre pass bas */
        iir_filter.set(4, 1e3, 0, SAMPLE_FREQUENCY, FilterGabarit::LOW_PASS, AnalogFilter::BUTTERWORTH);
        iir_filter.setup();
        std::cout << "Filtre passe bas :\n";
        iir_filter.printCoefficients();
        (*signalLP) = iir_filter.filtering((*signal));


        IIRFilter dem_filter;
        double f_coupure = 3e3;
        double f_oscillator = 10e3;
        initDemodulate(dem_filter, f_coupure, SAMPLE_FREQUENCY);
        demodulateSignal((*signal), dem_filter, f_oscillator, (*signal_demAmpli), (*signal_demPhase));


        (*fft_signal_demAmpli) = signal_demAmpli->fft();
        (*fft_signal_demPhase) = signal_demPhase->fft();

        csv.write("test.csv");
        fft_csv.write("test_fft.csv");

        std::cout << "Test success finish\n";
    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

int main() {
    test_demodulate();

    return 0;
}