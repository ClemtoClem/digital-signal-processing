#include <iostream>
#include "globals.hpp"
#include "Signal.hpp"
#include "Spectrum.hpp"
#include "CSVFile.hpp"
#include "Filter.hpp"
#include "Noise.hpp"
#include "Demodulator.hpp"


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
    double dem_filter_freq = 3e3;
    
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
                        SetDecimation(std::stoi(value));
                    } else if (name == "dem_filter_freq" || name == "demff") {
                        dem_filter_freq = std::stod(value);
                    } else if (name == "buffsize" || name == "bs") {
                        SetBufferSize(std::stoull(value));
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

    std::cerr << "- - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
    std::cerr << "-------- Sampling -------" << std::endl;
    std::cerr << "Decimation  : " << DECIMATION << std::endl;
    std::cerr << "Sample freq : " << SAMPLING_FREQUENCY << std::endl;
    std::cerr << "Size buffer : " << BUFFER_SIZE << std::endl;
    std::cerr << "---- Generate signal ----" << std::endl;
    std::cerr << "Amplitude 1 : " << s1_amplitude << std::endl;
    std::cerr << "Frequency 1 : " << s1_frequency << std::endl;
    std::cerr << "Amplitude 2 : " << s2_amplitude << std::endl;
    std::cerr << "Frequency 2 : " << s2_frequency << std::endl;
    std::cerr << "Offset      : " << offset << std::endl;
    std::cerr << "Phase       : " << phase << std::endl;
    std::cerr << "------ Demodulation -----" << std::endl;
    std::cerr << "Freq filter : " << dem_filter_freq << std::endl;
    std::cerr << "- - - - - - - - - - - - - - - - - - - - - - - " << std::endl << std::endl;

    /* - - - - - - - - - - - - - - - - - - - - - - - */

    // signal generation and processing code here

    Signal s1;
    Signal s2;

    Signal signal_output        ("output(t)");
    Signal signal_demAmpli      ("amplitude(t)");
    Signal signal_demPhase      ("phase(t)");
    Signal signalLP             ("signalFiltré(t)");

    Spectrum spectrum_output    ("OUTPUT(f)");
    Spectrum spectrum_demAmpli  ("AMPLITUDE(f)");


    /* - - - - - - - - - - - - - - - - - - - - - - - */
    /* synthétisation de la forme du signal */

    std::cerr << "Generate waveform" << std::endl;
    
    s1.generateWaveform(WaveformType::SINUS, s1_amplitude, s1_frequency);
    s2.generateWaveform(WaveformType::SINUS, s2_amplitude, s2_frequency);
    
    signal_output = s1+s2;
    std::cerr << std::endl;

    /* Bruiter le signal */

    WhiteNoise noise;
    noise.setGain(0.05);
    signal_output = noise.apply(signal_output);

    /* - - - - - - - - - - - - - - - - - - - - - - - */
    /* Test du filtre IIR passe bas de butterworth */

    std::cerr << "Filtrage" << std::endl;

    IIRFilter iir_filter;

    /* Filtre pass bas */
    iir_filter.set(4, 20e3, 0, FilterGabarit::LOW_PASS, AnalogFilter::BUTTERWORTH);
    iir_filter.setup();
    std::cerr << "Filtre passe bas :\n";
    iir_filter.printCoefficients();
    signalLP = iir_filter.apply(signal_output);


    std::cerr << "----------- Démodulation du signal -----------" << std::endl;

    /* Initialisation de la démodulation */
    double oscillator_frequency = s1_frequency;
    Demodulator dem(dem_filter_freq, oscillator_frequency);
    dem.setup();

    dem.demodulate(signal_output, signal_demAmpli, signal_demPhase);
    std::cerr << std::endl;

    std::cerr << "--------------- Calcul des FFT ---------------" << std::endl;
    
    size_t low_index, high_index;
    double T_per_period = 1.0 / s1_frequency;
	float rising_time = signal_demAmpli.getRisingTime(low_index, high_index);
    high_index = low_index + (rising_time/T_per_period) * 5;

    signal_output.FFT(spectrum_output);
    signal_demAmpli.FFT(spectrum_demAmpli, high_index);

    std::cerr << "Rising time: " << rising_time << std::endl;
    std::cerr << std::endl;


    std::cerr << "----------------------------------------------" << std::endl;
    // Calculer et afficher le niveau de bruit RMS
    double noiseRMS = signal_output.calculateNoiseRMS();
    std::cerr << "Niveau de bruit RMS: " << noiseRMS << std::endl;
    std::cerr << std::endl;
    

    /* Sauvgarde des signaux */
    std::cerr << "------------- Write data to file -------------" << std::endl;
	std::cerr << "File: ./data/test_temporel.csv" << std::endl;
	CSVFile outFile1("./data/test_temporel.csv");
	std::vector<Signal> outSpSig = {signal_output, signal_demAmpli, signal_demPhase};
	outFile1.writeSignals(outSpSig, true);

	std::cerr << "File: ./data/test_FFT_signal.csv" << std::endl;
	CSVFile outFile6("./data/test_FFT_signal.csv");
	std::vector<Spectrum> outSpectrum = {spectrum_output};
	outFile6.writeSpectrums(outSpectrum, true);

	std::cerr << "File: ./data/test_FFT_amplitude.csv" << std::endl;
	CSVFile outFile4("./data/test_FFT_amplitude.csv");
	std::vector<Spectrum> outSpectrum2 = {spectrum_demAmpli};
	outFile4.writeSpectrums(outSpectrum2, true);
	std::cerr << "----------------------------------------------" << std::endl << std::endl;

    std::cerr << "Test success finish" << std::endl;

    return 0;
} catch (std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
}

int main(int argc, char *argv[]) {
    return test_demodulate(argc, argv);
}