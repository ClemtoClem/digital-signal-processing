#include <iostream>
#include "Signal.hpp"
#include "CSV_File.hpp"
#include "Filter.hpp"

const int SAMPLE_FREQUENCY = 240000;
const int SIGNAL_SIZE = 1024;

void test_signals() {
    Waveform wf_1 {
        .type = WaveformType::SINUS,
        .amplitude = 0.8,
        .frequency = 1000,
        .phase = 0.0,
        .offset = 0.0,
        .delay = 0.0,
        .duty_cycle = 0.0,
    };
    
    Waveform wf_2 {
        .type = WaveformType::TRIANGLE,
        .amplitude = 0.8,
        .frequency = 1000,
        .phase = 0.0,
        .offset = 0.0,
        .delay = 0.0,
        .duty_cycle = 0.0
    };

    Signal signal1(SIGNAL_SIZE, SAMPLE_FREQUENCY);
    Signal signal2(SIGNAL_SIZE, SAMPLE_FREQUENCY);

    signal1.setWaveform(wf_1);
    signal1.setWaveform(wf_2);

    signal1.generate();
    signal2.generate();

    Signal signal3 = signal1+signal2;

    std::cerr << "signal3 = [";
    for (size_t i = 0; i<signal3.size(); i++) {
        std::cerr << signal3[i] << " ";
    }
    std::cerr << "]\n";
}

void test_csv_file() {
    try {
        CSV_File csv(SAMPLE_FREQUENCY);
        Signal signal1(SIGNAL_SIZE, SAMPLE_FREQUENCY); //auto signal1 = ccsv.create("s1", SIGNAL_SIZE, 0);
        Signal signal2(SIGNAL_SIZE, SAMPLE_FREQUENCY); //auto signal2 = csv.create("s2", SIGNAL_SIZE, 0);
        auto signal3 = csv.create("s3=s1*s2", SIGNAL_SIZE, 0);
        auto signal_filtre = csv.create("filtre(s3)", SIGNAL_SIZE, 0);
        //auto demodulate = csv.create("demodulate(s3)", SIGNAL_SIZE, 0);

        CSV_File fft_csv(SAMPLE_FREQUENCY);
        auto fft_s3 = fft_csv.create("fft(s3)", SIGNAL_SIZE, 0);
        auto fft_filter = fft_csv.create("filter_FIR", SIGNAL_SIZE, 0);
        //auto fft_dem_s3 = fft_csv.create("fft(demodulate(s3))", SIGNAL_SIZE, 0);

        Waveform wf_1 {
            .type = WaveformType::SINUS,
            .amplitude = 1.0,
            .frequency = 1000,
            .phase = 0.0,
            .offset = 0.0,
            .delay = 0.0,
            .duty_cycle = 0.0,
        };
        
        Waveform wf_2 {
            .type = WaveformType::SINUS,
            .amplitude = 0.2,
            .frequency = 10000,
            .phase = 0.0,
            .offset = 0.0,
            .delay = 0.0,
            .duty_cycle = 0.0
        };
        
        signal1.setWaveform(wf_1);
        signal2.setWaveform(wf_2);

        signal1.generate();
        signal2.generate();

        (*signal3) = signal1 + signal2;
        //(*demodulate) = signal1->demodulate(10000);
        //(*fft_dem_s3) = demodulate->fft();

        IIRFilter filter1;
        filter1.butterworthLowPass(4, 1000, SAMPLE_FREQUENCY);
        (*signal_filtre) = filter1.filter((*signal3));

        (*fft_s3) = signal3->fft();


        csv.write("test.csv");
        fft_csv.write("test_fft.csv");
    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

int main() {
    test_csv_file();

    return 0;
}