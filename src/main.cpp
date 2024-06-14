#include <iostream>
#include "Signal.hpp"
#include "CSVFile.hpp"
#include "Filter.hpp"
#include "Noise.hpp"
#include "Demodulator.hpp"

const int MAX_SAMPLING_FREQ = 250e6;
const uint32_t MAX_SIGNAL_SIZE = 16384;

int parseExponentNumber(const std::string &expStr) {
    double number;
    std::stringstream ss(expStr);
    
    ss >> number;
    
    if (ss.fail()) {
        std::cerr << "Error: Failed to parse the number." << std::endl;
        return 0;
    }

    // Convert to an integer
    int integerNumber = static_cast<int>(round(number));
    return integerNumber;
}

int test_demodulate(int argc, char *argv[]) try {
    std::string data_filename = "./data/test";
    double s1_frequency = 10e3;
    double s2_frequency = 100e3;
    double offset       = 0;
    double s1_amplitude = 0.8;
    double s2_amplitude = 0.2;
    double phase        = 0;
    int    decimation   = 32;
    double dem_filter_freq = 3e3;
    uint32_t signal_size = MAX_SIGNAL_SIZE;
    
    if (argc > 1) {
        // Parcourir chaque argument à partir du deuxième (le premier étant le nom du programme)
        for (int i = 1; i < argc; i++) {
            std::string param = argv[i];

            // Vérifier si l'argument est "help"
            if (param == "help") {
                std::cerr << "\033[4;0mHelp message\033[0m" << std::endl;
                std::cerr << "Details:" << std::endl;
                std::cerr << "  This application generated numeric signals, in order to perform" << std::endl;
                std::cerr << "  demodulation and fft calculation." << std::endl;
                std::cerr << "  The time and frequency signals are written to a csv file." << std::endl;
                std::cerr << "Usage: " << argv[0] << std::endl;
                std::cerr << "  amplitude1=<valeur>\t\ta1=<valeur>\t\t" << std::endl;
                std::cerr << "  frequency1=<valeur>\t\tf1=<valeur>\t\t" << std::endl;
                std::cerr << "  amplitude2=<valeur>\t\ta2=<valeur>\t\t" << std::endl;
                std::cerr << "  frequency2=<valeur>\t\tf2=<valeur>\t\t" << std::endl;
                std::cerr << "  offset=<valeur>\t\toff=<valeur>\t\t" << std::endl;
                std::cerr << "  phase=<value>\t\tdec=<value>\t\t" << std::endl;
                std::cerr << "  decimation=<value>\t\tph=<value>\t\t" << std::endl;
                std::cerr << "  buffsize=<value>\t\tbs=<value>\t\t" << std::endl;
                std::cerr << "  dem_filter_freq=<value>\tdemff=<value>\t\t" << std::endl;
                std::cerr << "  tfilename=<string>\t\tf=<string>\t\t";
                return 0;
            } else {
                /// Parse other arguments
                size_t pos = param.find('=');
                if (pos != std::string::npos) {
                    std::string name = param.substr(0, pos);
                    std::string value = param.substr(pos + 1);
                    
                    // Compare the parameter name to retrieve the value
                    if (name == "amplitude1" || name == "a1") {
                        s1_amplitude = std::stod(value);
                    } else if (name == "frequency1" || name == "f1") {
                        s1_frequency = parseExponentNumber(value);
                    } else if (name == "amplitude2" || name == "a2") {
                        s2_amplitude = std::stod(value);
                    } else if (name == "frequency2" || name == "f2") {
                        s2_frequency = parseExponentNumber(value);
                    } else if (name == "offset" || name == "off") {
                        offset = std::stod(value);
                    } else if (name == "phase" || name == "ph") {
                        phase = std::stod(value);
                    } else if (name == "decimation" || name == "dec") {
                        decimation = std::stoi(value);
                    } else if (name == "dem_filter_freq" || name == "demff") {
                        dem_filter_freq = std::stod(value);
                    } else if (name == "buffsize" || name == "bs") {
                        signal_size = std::stoul(value);
                        if (signal_size > MAX_SIGNAL_SIZE)
                            signal_size = MAX_SIGNAL_SIZE;
                    } else if (name == "filename" || name == "f") {
                        data_filename = value;
                    }
                } else {
                    std::cerr << "Invalid argument format: " << param << std::endl;
                    return 1;
                }
            }
        }
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - */

    const uint32_t SIGNAL_SIZE = signal_size;
    const int SAMPLING_FREQ = (double) MAX_SAMPLING_FREQ / decimation;

    std::cerr << "- - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
    std::cerr << "-------- Sampling -------" << std::endl;
    std::cerr << "Decimation  : " << decimation << std::endl;
    std::cerr << "Sample freq : " << SAMPLING_FREQ << std::endl;
    std::cerr << "Size buffer : " << SIGNAL_SIZE << std::endl;
    std::cerr << "---- Generate signal ----" << std::endl;
    std::cerr << "Amplitude 1 : " << s1_amplitude << std::endl;
    std::cerr << "Frequency 1 : " << s1_frequency << std::endl;
    std::cerr << "Amplitude 2 : " << s2_amplitude << std::endl;
    std::cerr << "Frequency 2 : " << s2_frequency << std::endl;
    std::cerr << "Offset      : " << offset << std::endl;
    std::cerr << "Phase       : " << phase << std::endl;
    std::cerr << "------ Demodulation -----" << std::endl;
    std::cerr << "Freq filter : " << dem_filter_freq << std::endl;
    std::cerr << "- - - - - - - - - - - - - - - - - - - - - - - " << std::endl;

    /* - - - - - - - - - - - - - - - - - - - - - - - */

    // signal generation and processing code here

    Signal s1(SIGNAL_SIZE, SAMPLING_FREQ);
    Signal s2(SIGNAL_SIZE, SAMPLING_FREQ);

    Signal signal_output        (SIGNAL_SIZE, SAMPLING_FREQ, "output(t)");
    //Signal signal_input       (SIGNAL_SIZE, SAMPLING_FREQ, "input(t)");
    Signal signal_demAmpli      (SIGNAL_SIZE, SAMPLING_FREQ, "amplitude(t)");
    Signal signal_demPhase      (SIGNAL_SIZE, SAMPLING_FREQ, "phase(t)");
    Signal signalLP             (SIGNAL_SIZE, SAMPLING_FREQ, "signalFiltré(t)");

    Signal DFT_signal_output    (SIGNAL_SIZE, SAMPLING_FREQ, "DFT_output(t)");
    //Signal DFT_signal_input   (SIGNAL_SIZE, SAMPLING_FREQ, "DFT_input(t)");
    Signal DFT_signal_demAmpli  (SIGNAL_SIZE, SAMPLING_FREQ, "DFT_amplitude(t)");
    Signal DFT_signal_demPhase  (SIGNAL_SIZE, SAMPLING_FREQ, "DFT_phase(t)");


    /* - - - - - - - - - - - - - - - - - - - - - - - */
    /* synthétisation de la forme du signal */

    std::cout << "Generate waveform" << std::endl;
    
    s1.generateWaveform(WaveformType::SINUS, s1_amplitude, s1_frequency);
    s2.generateWaveform(WaveformType::SINUS, s2_amplitude, s2_frequency);
    
    signal_output = s1+s2;

    std::cerr << "- - - - - - - - - - - - - - - - - - - - - - - " << std::endl;

    /* - - - - - - - - - - - - - - - - - - - - - - - */
    /* Bruiter le signal */

    Noise noise;
    noise.setParams(NoiseType::WHITE, SAMPLING_FREQ, 0.05);
    signal_output = noise.process(signal_output);

    /* - - - - - - - - - - - - - - - - - - - - - - - */
    /* Test du filtre IIR passe bas de butterworth */

    std::cout << "Filtrage" << std::endl;

    IIRFilter iir_filter;

    /* Filtre pass bas */
    iir_filter.set(4, 20e3, 0, SAMPLING_FREQ, FilterGabarit::LOW_PASS, AnalogFilter::BUTTERWORTH);
    iir_filter.setup();
    std::cout << "Filtre passe bas :\n";
    iir_filter.printCoefficients();
    signalLP = iir_filter.process(signal_output);

    std::cerr << "- - - - - - - - - - - - - - - - - - - - - - - " << std::endl;

    /* - - - - - - - - - - - - - - - - - - - - - - - */
    /* Initialisation de la démodulation */
    
    std::cout << "Demodulation :\n";

    double oscillator_frequency = s1_frequency;
    Demodulator dem(dem_filter_freq, oscillator_frequency, SAMPLING_FREQ);
    dem.setup();

    std::cerr << "- - - - - - - - - - - - - - - - - - - - - - - " << std::endl;

    /* - - - - - - - - - - - - - - - - - - - - - - - */
    /* Démodulation du signal */

    dem.demodulate(signal_output, signal_demAmpli, signal_demPhase);

    /* - - - - - - - - - - - - - - - - - - - - - - - */
    /* Calcul des transformées de fourier discrètes des signaux avec buff_size zero padding */

    std::cerr << "Calcul des transformées de fourier discrètes" << std::endl;

    DFT_signal_output   = signal_output.DFT(SIGNAL_SIZE);
    DFT_signal_demAmpli = signal_demAmpli.DFT(SIGNAL_SIZE);
    DFT_signal_demPhase = signal_demPhase.DFT(SIGNAL_SIZE);

    std::cerr << "- - - - - - - - - - - - - - - - - - - - - - - " << std::endl;

    /* - - - - - - - - - - - - - - - - - - - - - - - */
    /* Sauvgarde des signaux */

    std::cout << "Sauvgarde des signaux" << std::endl;

    std::cout << data_filename << ".csv" << std::endl;
    CSVFile outFile(data_filename+".csv");
    std::vector<Signal> outSignals;
    outSignals.emplace_back(signal_output);
    outSignals.emplace_back(signalLP);
    outSignals.emplace_back(signal_demAmpli);
    outSignals.emplace_back(signal_demPhase);
    outFile.writeSignals(outSignals, false); // with time axis

    std::cout << data_filename << "_DFT.csv" << std::endl;
    CSVFile outFileDFT(data_filename+"_DFT.csv");
    std::vector<Signal> outSignalsDFT;
    outSignalsDFT.emplace_back(DFT_signal_output);
    outSignalsDFT.emplace_back(DFT_signal_demAmpli);
    outSignalsDFT.emplace_back(DFT_signal_demPhase);
    outFileDFT.writeSignals(outSignalsDFT, true); // with frequency axis

    std::cerr << "- - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
    std::cout << "Test success finish" << std::endl;

    return 0;
} catch (std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
}

int main(int argc, char *argv[]) {
    return test_demodulate(argc, argv);
}