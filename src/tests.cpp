#include "tests.hpp"
#include "CSVFile.hpp"
#include "utils.hpp"
#include <stdexcept>



int test_simple(const std::vector<std::string> &args) {
	int result;
	
	try {

		WindowType window_type = WindowType::Rectangular;
		WaveformType waveform_type = WaveformType::SINUS;
		int    frequency = 10000;
		double  amplitude = 1.0;
		double  offset = 0;
		double  phase  = 0;
		double  dutycycle = 0.5;
		bool   hasSetDecimation = false;
		int    points_per_period = 100;

		if (args.size() >= 1) {
			for (auto param : args) {

				// Vérifier si l'argument est "help"
				if (param == "help") {
					std::cerr << "\033[4;0mHelp message\033[0m" << std::endl;
					std::cerr << "Details:" << std::endl;
					std::cerr << "  This application acquires on analog signal in channel 1 or channel 2, a signal" << std::endl;
					std::cerr << "  generated on analog output 1, in order to perform fft calculation." << std::endl;
					std::cerr << "  The time and frequency signals are written to a csv file." << std::endl;
					std::cerr << "Possibilities of utilisation : " << std::endl;
					std::cerr << "  waveform=<string>\t\twf=<string>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is " << waveformTypeToString(waveform_type) << "." << std::endl;
					std::cerr << "  amplitude=<value_>\t\tamp=<value>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is " << amplitude << "." << std::endl;
					std::cerr << "  frequency=<integer>\tf=<integer>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is " << frequency << "." << std::endl;
					std::cerr << "  offset=<valeur>\t\toff=<valeur>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is " << offset << "." << std::endl;
					std::cerr << "  phase=<value>\t\t\tph=<value>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is " << phase << "." << std::endl;
					std::cerr << "  dutycycle=<value>\tdc=<value>\t\t" << std::endl;
					std::cerr << "  window=<string>\twin=<string>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is " << windowTypeToString(window_type) << "." << std::endl;
					std::cerr << "  points_per_period=<integer>\tppp=<integer>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is "<< points_per_period << "." << std::endl;
					std::cerr << "  decimation=<integer>\t\tdec=<integer>\t\t" << std::endl;
					std::cerr << "      note: the 'points_per_period' and 'decimation' arguments allow decimation to be defined differently." << std::endl;
					std::cerr << "            The 'points_per_period' argument defines the number of points per period, and the 'decimation' argument defines the decimation factor." << std::endl;
					std::cerr << "  buffsize=<integer>\t\tbs=<integer>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is " << BUFFER_SIZE << "." << std::endl;
					return 0;
				} else {
					// Parse other arguments
					size_t pos = param.find('=');
					if (pos != std::string::npos) {
						std::string name = param.substr(0, pos);
						std::string value = param.substr(pos + 1);

						for (size_t i = 0; i < name.size(); ++i) {
							if (name[i] <= 'Z' && name[i] >= 'A') {
								name[i] = name[i] + 32;
							}
						}

						if (name == "waveform" || name == "wf") {
                            waveform_type = stringToWaveformType(value);
                        } else if (name == "amplitude" || name == "amp") {
                            amplitude = std::stof(value);
                        } else if (name == "frequency" || name == "f") {
                            frequency = convertToInteger(value);
                        } else if (name == "offset" || name == "off") {
                            offset = std::stof(value);
                        } else if (name == "phase" || name == "ph") {
                            phase = std::stof(value);
                        } else if (name == "dutycycle" || name == "dc") {
                            dutycycle = std::stof(value);
                        } else if (name == "window" || name == "win") {
                            window_type = stringToWindowType(value);
                        } else if (name == "points_per_period" || name == "ppp") {
                            points_per_period = std::stoi(value);
                        } else if (name == "decimation" || name == "dec") {
							hasSetDecimation = true;
                            SetDecimation(std::stoi(value));
                        } else if (name == "buffsize" || name == "bs") {
                            SetBufferSize(std::stoi(value));
                        } else {
                            std::cerr << "Invalid argument format: " << param << std::endl;
                            return 1;
                        }
                    } else {
                        std::cerr << "Invalid argument format: " << param << std::endl;
                        return 1;
                    }
				}
			}
		}

		if (!hasSetDecimation) {
			SetDecimation(calculateDecimation(frequency, points_per_period));
		}
		
		Signal signal("signal(t)");
		Signal windowedSignal;
		Spectrum spectrum("SIGNAL(f)");

		signal.generateWaveform(waveform_type, amplitude, frequency, phase, offset, 0, dutycycle*100);


		/* Initialisation de la fenêtre */
		Window window;
		if (!window.set(window_type, BUFFER_SIZE)) {
			std::cerr << "Error: Unable to set the window function." << std::endl;
			return 1;
		}
		window.setup();

		std::cerr << "------ Start ------" << std::endl;
		std::cerr << "waveform=" << waveformTypeToString(waveform_type) << std::endl;
		std::cerr << "amplitude=" << amplitude << std::endl;
		std::cerr << "frequency=" << frequency << std::endl;
		std::cerr << "offset=" << offset << std::endl;
		std::cerr << "phase=" << phase << std::endl;
		std::cerr << "dutycycle=" << dutycycle << std::endl;
		std::cerr << "window=" << windowTypeToString(window_type) << std::endl;
		double ppp = SAMPLING_FREQUENCY / static_cast<double>(frequency);
		std::cerr << "points_per_period=" << ppp << std::endl;
		std::cerr << "decimation=" << DECIMATION << std::endl;
		std::cerr << "buffsize=" << BUFFER_SIZE << std::endl;

		windowedSignal = window.apply(signal);
		windowedSignal.FFT(spectrum);

		std::cerr << "---- Save data ----" << std::endl;
		CSVFile::setTimeToFilepath();
		std::string filename = "test.csv";
		std::cerr << "File: " << filename << std::endl;
		CSVFile outFile1(filename);
		std::vector<Signal> outSig = {signal};
		outFile1.writeSignals(outSig, true);

		filename = "test_FFT.csv";
		std::cerr << "File: " << filename << std::endl;
		CSVFile outFile2(filename);
		std::vector<Spectrum> outSpec = {spectrum};
		outFile2.writeSpectrums(outSpec, true);
		std::cerr << "------- End -------" << std::endl;

		result = 0;
	} catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        result = 1;
    }

	return result;
}

int test_spectrum(const std::vector<std::string> &args)
{
	int result;
	
	try {
		SetBufferSize(MAX_BUFFER_SIZE);
		SetDecimation(16);

		std::vector<int> frequencies  = {10000, 100000};
		std::vector<double> amplitudes = {0.8, 0.2};
		double offset = 0;
		double phase  = 0;
		double noise = 0.0;
		WindowType window_type  = WindowType::Hamming;
		bool hasSetDecimation = false;
		int points_per_period = 100;

		if (args.size() >= 1) {
			for (auto param : args) {

				// Vérifier si l'argument est "help"
				if (param == "help") {
					std::cerr << "\033[4;0mHelp message\033[0m" << std::endl;
					std::cerr << "Details:" << std::endl;
					std::cerr << "  This application acquires on analog input 1, a signal" << std::endl;
					std::cerr << "  generated on analog output 1, in order to perform" << std::endl;
					std::cerr << "  fft calculation." << std::endl;
					std::cerr << "  The time and frequency signals are written to a csv file." << std::endl;
					std::cerr << "Possibilities of utilisation : " << std::endl;
					std::cerr << "  amplitudes=<value_1,value_2,...,value_n>\t\tamp=<value_1,value_2,...,value_n>\t\t" << std::endl;
					std::cerr << "  frequencies=<integer_1,integer_2,...,integer_n>\tf=<integer_1,integer_2,...,integer_n>\t\t" << std::endl;
					std::cerr << "      note: The number of frequencies and amplitudes entered must be the same." << std::endl;
					std::cerr << "            Values ​​must be separated by a comma." << std::endl;
					std::cerr << "  offset=<valeur>\t\toff=<valeur>\t\t" << std::endl;
					std::cerr << "      note: This argument is optional, and if not entered, the default value is " << offset << "." << std::endl;
					std::cerr << "  phase=<value>\t\t\tph=<value>\t\t" << std::endl;
					std::cerr << "      note: This argument is optional, and if not entered, the default value is " << phase << "." << std::endl;
					std::cerr << "  noise=<value>\tns=<value>\t\t" << std::endl;
					std::cerr << "      note: This argument is optional, and if not entered, the default value is " << noise << "." << std::endl;
					std::cerr << "  window=<string>\twin=<string>\t\t" << std::endl;
					std::cerr << "      note: This argument is optional, and if not entered, the default value is " << windowTypeToString(window_type) << "." << std::endl;
					std::cerr << "  points_per_period=<integer>\tppp=<integer>\t\t" << std::endl;
					std::cerr << "      note: This argument is optional, and if not entered, the default value is " << points_per_period << "." << std::endl;
					std::cerr << "  decimation=<integer>\t\tdec=<integer>\t\t" << std::endl;
					std::cerr << "      note: The 'points_per_period' and 'decimation' arguments allow decimation to be defined differently." << std::endl;
					std::cerr << "            The 'points_per_period' argument defines the number of points per period, and the 'decimation' argument defines the decimation factor." << std::endl;
					std::cerr << "  buffsize=<integer>\t\tbs=<integer>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is " << BUFFER_SIZE << std::endl;
					return 0;
				} else {
					// Parse other arguments
					size_t pos = param.find('=');
					if (pos != std::string::npos) {
						std::string name = param.substr(0, pos);
						std::string value = param.substr(pos + 1);

						for (size_t i = 0; i < name.size(); ++i) {
							if (name[i] <= 'Z' && name[i] >= 'A') {
								name[i] = name[i] + 32;
							}
						}

						if (name == "amplitudes" || name == "amp") {
							std::vector<std::string> values = split(value, ",");
							amplitudes.clear();
							for (const std::string& v : values) {
								amplitudes.push_back(std::stof(v));
							}
						} else if (name == "frequencies" || name == "f") {
							std::vector<std::string> values = split(value, ",");
							frequencies.clear();
							for (const std::string& v : values) {
								frequencies.push_back(convertToInteger(v));
							}
						} else if (name == "offset" || name == "off") {
							offset = std::stof(value);
						} else if (name == "phase" || name == "ph") {
							phase = std::stof(value);
						} else if (name == "noise" || name == "ns") {
							noise = std::stof(value);
						} else if (name == "window" || name == "win") {
							window_type = stringToWindowType(value);
						} else if (name == "points_per_period" || name == "ppp") {
							points_per_period = convertToInteger(value);
						} else if (name == "decimation" || name == "dec") {
							hasSetDecimation = true;
							SetDecimation(std::stoi(value));
						} else if (name == "buffsize" || name == "bs") {
							SetBufferSize(convertToInteger(value));
						} else {
							std::cerr << "Invalid argument format: " << param << std::endl;
							return 1;
						}
					} else {
						std::cerr << "Invalid argument format: " << param << std::endl;
						return 1;
					}
				}
			}
		}
		if (amplitudes.size() == 0 || frequencies.size() == 0 || amplitudes.size() != frequencies.size()) {
			std::cerr << "Error: Amplitudes and frequencies must be specified and have the same size." << std::endl;
			return 1;
		}
		
		int max_frequency = 0;
		for (double f : frequencies) {
			if (f > max_frequency) {
				max_frequency = f;
			}
		}
		
		if (!hasSetDecimation) {
			SetDecimation(calculateDecimation(max_frequency, points_per_period));
		}

		/* - - - - - - - - - - - - - - - - - - - - - - - */
		/* Affichage des paramètres*/
		std::cerr << "+--------- Parameters ---------+" << std::endl;
		std::cerr << "|           Sampling           |" << std::endl;
		std::cerr << "| fs         : " << std::setw(15) << SAMPLING_FREQUENCY << " |" << std::endl;
		std::cerr << "| decimation : " << std::setw(15) << DECIMATION << " |" << std::endl;
		std::cerr << "| buffsize   : " << std::setw(9)  << BUFFER_SIZE << "/16384 |" << std::endl;
		double T_per_period = 1.0 / max_frequency;
		double ppp = T_per_period * SAMPLING_FREQUENCY;
		std::cerr << "| Samples per period : " << std::setw(7) << ppp << " |" << std::endl;
		std::cerr << "+------------------------------+" << std::endl;
		std::cerr << "|          Generator           |" << std::endl;
		std::cerr << "| amplitudes :                 |" << std::endl;
		for (double a : amplitudes) {
			std::cerr << "|         " << std::setw(15) << a << " V   |" << std::endl;
		}
		std::cerr << "| frequencies:                 |" << std::endl;
		for (int f : frequencies) {
			std::cerr << "|         " << std::setw(15) << f << " Hz  |" << std::endl;
		}
		std::cerr << "| noise      : " << std::setw(11) << noise       << "     |" << std::endl;
		std::cerr << "| offset     : " << std::setw(11) << offset      << " V   |" << std::endl;
		std::cerr << "| phase      : " << std::setw(11) << phase       << " rad |" << std::endl;
		std::cerr << "| window     : " << std::setw(11) << windowTypeToString(window_type) << "     |" << std::endl;
		std::cerr << "+------------------------------+" << std::endl;

		/* - - - - - - - - - - - - - - - - - - - - - - - */

		Signal signal("signal(t)");
		Signal windowedSignal;
		Spectrum spectrum("spectrum(f)");

		WhiteNoise whiteNoise;
		whiteNoise.setGain(noise);

		/* Signal arbitraire formé d'une somme de signaux sinusoïdaux */
		double t;
		for (size_t k = 0; k < BUFFER_SIZE; ++k) {
			t = (k / SAMPLING_FREQUENCY) + offset;
			signal[k] = 0;
			for (size_t i = 0; i < amplitudes.size(); ++i) {
				signal[k] += amplitudes[i] * sin(2 * M_PI * frequencies[i] * t + phase) + offset;
			}
			signal[k] = whiteNoise.apply(signal[k]);
		}

		/* - - - - - - - - - - - - - - - - - - - - - - - */

		/* Initialisation de la fenêtre */
		Window window;
		if (!window.set(window_type, BUFFER_SIZE)) {
			std::cerr << "Error: Unable to set the window function." << std::endl;
			return 1;
		}
		window.setup();
		
		/* - - - - - - - - - - - - - - - - - - - - - - - */
		/* Echantillonnage du signal */
		windowedSignal = window.apply(signal);
		std::cerr << "| Windowing done               |" << std::endl;
		std::cerr << "+------------------------------+" << std::endl << std::endl;

		/* Calcul des transformées de fourier discrètes des signaux avec BUFFER_SIZE zero padding */
		
		std::cerr << "+------- Calculate FFT --------+" << std::endl;
		windowedSignal.FFT(spectrum);
		std::cerr << "| End of FFT calculation       |" << std::endl;
		std::cerr << "+------------------------------+" << std::endl << std::endl;

		/* - - - - - - - - - - - - - - - - - - - - - - - */
		/* Sauvgarde des résultats */
		CSVFile::setTimeToFilepath();
		std::string filename = "test_spectrum";

		std::cerr << "+--------------- Write results ---------------+" << std::endl;
		std::cerr << "| File: " << std::setw(27) << filename << "_signals.csv |" << std::endl;
		CSVFile outFile1(filename+"_signals.csv");
		std::vector<Signal> outSpSig = {signal};
		outFile1.writeSignals(outSpSig, true);

		std::cerr << "| File: " << std::setw(23) << filename << "_FFT.csv |" << std::endl;
		CSVFile outFile6(filename+"_FFT.csv");
		std::vector<Spectrum> outSpectrum = {spectrum};
		outFile6.writeSpectrums(outSpectrum, true);
		std::cerr << "+----------------------------------------------+" << std::endl << std::endl;

		result = 0;
	} catch (const std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		result = 1;
	}

	return result;
}

/* - - - - - - - - - - - - - - - - - - - - - - - */

int test_demodulation(const std::vector<std::string> &args) {
	int result;
	
	try {
		srand( time( NULL ) );

		SetBufferSize(MAX_BUFFER_SIZE);
		SetDecimation(16);

		std::vector<int> frequencies  = {10000, 100000};
		std::vector<double> amplitudes = {0.8, 0.2};
		double offset  = 0;
		double phase   = 0;
		double noise   = 0.0;
		WindowType window_type  = WindowType::Hamming;
		double dem_filter_freq = 3e3;
		bool hasSetDecimation = false;
		int points_per_period = 100;
		
		if (args.size() >= 1) {
			for (auto param : args) {

				// Vérifier si l'argument est "help"
				if (param == "help") {
					std::cerr << "\033[4;0mHelp message\033[0m" << std::endl;
					std::cerr << "Details:" << std::endl;
					std::cerr << "  This application acquires on analog input 1, a signal" << std::endl;
					std::cerr << "  generated on analog output 1, in order to perform" << std::endl;
					std::cerr << "  demodulation and fft calculation." << std::endl;
					std::cerr << "  The time and frequency signals are written to a csv file." << std::endl;
					std::cerr << "Possibilities of utilisation : " << std::endl;
					std::cerr << "  amplitudes=<value_1,value_2,...,value_n>\t\tamp=<value_1,value_2,...,value_n>\t\t" << std::endl;
					std::cerr << "  frequencies=<integer_1,integer_2,...,integer_n>\tf=<integer_1,integer_2,...,integer_n>\t\t" << std::endl;
					std::cerr << "      note: The number of frequencies and amplitudes entered must be the same." << std::endl;
					std::cerr << "            Values ​​must be separated by a comma." << std::endl;
					std::cerr << "  offset=<valeur>\t\toff=<valeur>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is " << offset << "." << std::endl;
					std::cerr << "  phase=<value>\t\t\tph=<value>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is " << phase << "." << std::endl;
					std::cerr << "  noise=<value>\tns=<value>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is 0.0" << std::endl;
					std::cerr << "  window=<string>\twin=<string>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is Hamming." << std::endl;
					std::cerr << "  dem_filter_freq=<value>\tdemff=<value>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is 3000." << std::endl;
					std::cerr << "  points_per_period=<integer>\tppp=<integer>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is 100." << std::endl;
					std::cerr << "  decimation=<integer>\t\tdec=<integer>\t\t" << std::endl;
					std::cerr << "      note: the 'points_per_period' and 'decimation' arguments allow decimation to be defined differently." << std::endl;
					std::cerr << "            The 'points_per_period' argument defines the number of points per period, and the 'decimation' argument defines the decimation factor." << std::endl;
					std::cerr << "  buffsize=<integer>\t\tbs=<integer>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is 1024." << std::endl;
					return 0;
				} else {
					// Parse other arguments
					size_t pos = param.find('=');
					if (pos != std::string::npos) {
						std::string name = param.substr(0, pos);
						std::string value = param.substr(pos + 1);

						for (size_t i = 0; i < name.size(); ++i) {
							if (name[i] <= 'Z' && name[i] >= 'A') {
								name[i] = name[i] + 32;
							}
						}

						if (name == "amplitudes" || name == "amp") {
							std::vector<std::string> values = split(value, ",");
							amplitudes.clear();
							for (const std::string& v : values) {
								amplitudes.push_back(std::stof(v));
							}
						} else if (name == "frequencies" || name == "f") {
							std::vector<std::string> values = split(value, ",");
							frequencies.clear();
							for (const std::string& v : values) {
								frequencies.push_back(convertToInteger(v));
							}
						} else if (name == "offset" || name == "off") {
							offset = std::stof(value);
						} else if (name == "phase" || name == "ph") {
							phase = std::stof(value);
						} else if (name == "noise" || name == "ns") {
							noise = std::stof(value);
						} else if (name == "window" || name == "win") {
							window_type = stringToWindowType(value);
						} else if (name == "dem_filter_freq" || name == "demff") {
							dem_filter_freq = convertToInteger(value);
						} else if (name == "points_per_period" || name == "ppp") {
							points_per_period = convertToInteger(value);
						} else if (name == "decimation" || name == "dec") {
							hasSetDecimation = true;
							SetDecimation(std::stoi(value));
						} else if (name == "buffsize" || name == "bs") {
							SetBufferSize(convertToInteger(value));
						} else {
							std::cerr << "Invalid argument: " << param << std::endl;
							return 1;
						}
					} else {
						std::cerr << "Invalid argument format: " << param << std::endl;
						return 1;
					}
				}
			}
		}
		if (amplitudes.size() == 0 || frequencies.size() == 0 || amplitudes.size() != frequencies.size()) {
			std::cerr << "Error: Amplitudes and frequencies must be specified and have the same size." << std::endl;
			return 1;
		}
		
		int max_frequency = frequencies[0];
		int min_frequency = frequencies[0];
		for (double f : frequencies) {
			if (f > max_frequency) {
				max_frequency = f;
			}
			if (f < min_frequency) {
				min_frequency = f;
			}
		}
		
		if (!hasSetDecimation) {
			SetDecimation(calculateDecimation(max_frequency, points_per_period));
		}

		/* - - - - - - - - - - - - - - - - - - - - - - - */
		/* affichage des paramètres*/
		std::cerr << "+--------- Parameters ---------+" << std::endl;
		std::cerr << "|           Sampling           |" << std::endl;
		std::cerr << "| fs         : " << std::setw(15) << SAMPLING_FREQUENCY << " |" << std::endl;
		std::cerr << "| decimation : " << std::setw(15) << DECIMATION << " |" << std::endl;
		std::cerr << "| buffsize   : " << std::setw(9)  << BUFFER_SIZE << "/16384 |" << std::endl;
		double T_per_period = 1.0 / max_frequency;
		double ppp = T_per_period * SAMPLING_FREQUENCY;
		std::cerr << "| Samples per period : " << std::setw(7) << ppp << " |" << std::endl;
		std::cerr << "+------------------------------+" << std::endl;
		std::cerr << "|          Generator           |" << std::endl;
		std::cerr << "| amplitudes :                 |" << std::endl;
		for (double a : amplitudes) {
			std::cerr << "|         " << std::setw(12) << a << " V   |" << std::endl;
		}
		std::cerr << "| frequencies:                 |" << std::endl;
		for (int f : frequencies) {
			std::cerr << "|             " << std::setw(12) << f << " Hz  |" << std::endl;
		}
		std::cerr << "| noise      : " << std::setw(11) << noise       << "     |" << std::endl;
		std::cerr << "| offset     : " << std::setw(11) << offset      << " V   |" << std::endl;
		std::cerr << "| phase      : " << std::setw(11) << phase       << " rad |" << std::endl;
		std::cerr << "+------------------------------+" << std::endl;
		std::cerr << "|            Window            |" << std::endl;
		std::cerr << "| type       : " << std::setw(11) << windowTypeToString(window_type) << "     |" << std::endl;
		std::cerr << "+------------------------------+" << std::endl;
		std::cerr << "|         Demodulation         |" << std::endl;
		std::cerr << "| Freq filter: " << std::setw(11) << dem_filter_freq << " Hz  |" << std::endl;
		std::cerr << "+------------------------------+" << std::endl << std::endl;

		Signal signal("signal(t)");
		Signal windowedSignal;
		Signal signal_demAmpli("amplitude(t)");
		Signal signal_demPhase("phase(t)");

		Spectrum spectrum("SIGNAL(f)");
		Spectrum spectrum_demAmpli("AMPLITUDE(f)");

		/* - - - - - - - - - - - - - - - - - - - - - - - */

		WhiteNoise whiteNoise;
		whiteNoise.setGain(noise);

		/* Signal arbitraire formé d'une somme de signaux sinusoïdaux */
		double t;
		for (size_t k = 0; k < BUFFER_SIZE; ++k) {
			t = (k / SAMPLING_FREQUENCY) + offset;
			signal[k] = 0;
			for (size_t i = 0; i < amplitudes.size(); ++i) {
				signal[k] += amplitudes[i] * sin(2 * M_PI * frequencies[i] * t + phase) + offset;
			}
			signal[k] = whiteNoise.apply(signal[k]);
		}
		
		/* - - - - - - - - - - - - - - - - - - - - - - - */
		/* Initialisation de la démodulation */
		double freq_oscillator = static_cast<double>(min_frequency); // fréquence de l'oscillateur
		Demodulator dem(dem_filter_freq, freq_oscillator);
		dem.setup();

		/* Initialisation de la fenêtre */
		Window window;
		if (!window.set(window_type, BUFFER_SIZE)) {
			std::cerr << "Error: Unable to set the window function." << std::endl;
			return 1;
		}
		window.setup();
		
		windowedSignal = window.apply(signal);

		/* - - - - - - - - - - - - - - - - - - - - - - - */
		/* Démodulation du signal */
		dem.apply(signal, signal_demAmpli, signal_demPhase, true);

		/* - - - - - - - - - - - - - - - - - - - - - - - */
		/* Calcul des transformées de fourier discrètes des signaux avec BUFFER_SIZE zero padding */
		
		std::cerr << "+------- Calculate FFT --------+" << std::endl;

		windowedSignal.FFT(spectrum);		
		
		double rising_time = 3/dem_filter_freq;
		size_t index = rising_time * SAMPLING_FREQUENCY;

		std::cerr << "| Rising time: " << std::setw(11) << rising_time << " s  |" << std::endl;
		signal.FFT(spectrum_demAmpli, index);

		/* - - - - - - - - - - - - - - - - - - - - - - - */
		/* Sauvgarde des signaux */

		CSVFile::setTimeToFilepath();
		std::string filename = "test_demodulation";

		std::cerr << "+--------------- Write results ---------------+" << std::endl;
		std::cerr << "| File: " << std::setw(27) << filename << "_signals.csv |" << std::endl;
		CSVFile outFile1(filename+"_signals.csv");
		std::vector<Signal> outSpSig = {signal, windowedSignal, signal_demAmpli, signal_demPhase};
		outFile1.writeSignals(outSpSig, true);

		std::cerr << "| File: " << std::setw(23) << filename << "_FFT_signal.csv |" << std::endl;
		CSVFile outFile2(filename+"_FFT_signal.csv");
		std::vector<Spectrum> outSpectrum = {spectrum};
		outFile2.writeSpectrums(outSpectrum, true, false);

		std::cerr << "| File: " << std::setw(20) << filename << "_FFT_amplitude.csv |" << std::endl;
		CSVFile outFile3(filename+"_FFT_amplitude.csv");
		std::vector<Spectrum> outSpectrum2 = {spectrum_demAmpli};
		outFile3.writeSpectrums(outSpectrum2, true, false);
		std::cerr << "+----------------------------------------------+" << std::endl << std::endl;

		result = 0;
	} catch (const std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		result = 1;
	}
	
	return result;
}

/* - - - - - - - - - - - - - - - - - - - - - - - */

int test_demodulation2(const std::vector<std::string> &args) {
	int result;
	
	try {
		SetBufferSize(MAX_BUFFER_SIZE);

		int   frequency1        = 10000;
		int   frequency_min     = 100;
		int   frequency_max     = 10000;
		int   nb_frequency      = 21;
		double amplitude1        = 0.8;
		double amplitude2        = 0.2;
		double offset            = 0;
		double phase             = 0;
		double noise            = 0.0;
		WindowType window_type  = WindowType::Hamming;
		double dem_filter_freq  = 1e3;
		bool hasSetDecimation   = false;
		int points_per_period   = 100;
		
		if (args.size() >= 1) {
			for (auto param : args) {

				// Vérifier si l'argument est "help"
				if (param == "help") {
					std::cerr << "\033[4;0mHelp message\033[0m" << std::endl;
					std::cerr << "Details:" << std::endl;
					std::cerr << "  " << std::endl;
					std::cerr << "Possibilities of utilisation : " << std::endl;
					std::cerr << "  frequency1=<valeur>\t\tf1=<valeur>\t\t" << std::endl;
					std::cerr << "  frequency_min=<valeur>\t\tfmin=<valeur>\t\t" << std::endl;
					std::cerr << "  frequency_max=<valeur>\t\tfmax=<valeur>\t\t" << std::endl;
					std::cerr << "  nb_frequency=<valeur>\t\tnf=<valeur>\t\t" << std::endl;
					std::cerr << "  amplitude1=<valeur>\t\ta1=<valeur>\t\t" << std::endl;
					std::cerr << "  amplitude2=<valeur>\t\ta2=<valeur>\t\t" << std::endl;
					std::cerr << "  offset=<valeur>\t\toff=<valeur>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is 0.0" << std::endl;
					std::cerr << "  phase=<value>\t\t\tph=<value>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is 0.0" << std::endl;
					std::cerr << "  noise=<value>\t\t\tn=<value>\t\t" << std::endl;
					std::cerr << "      note: this argument is optional, and if not entered, the default value is 0.0" << std::endl;
					std::cerr << "  window_type=<value>\t\twt=<value>\t\t" << std::endl;
					std::cerr << "  points_per_period=<value>\t\tppp=<value>\t\t" << std::endl;
					std::cerr << "      details: Number of points per period of the sampled signal" << std::endl;
					std::cerr << "  decimation=<integer>\t\tdec=<integer>\t\t" << std::endl;
					std::cerr << "      details : 125MHz quartz frequency decimation factor" << std::endl;
					std::cerr << "      note: decimation must be a power of 2" << std::endl;
					std::cerr << "  buffsize=<value>\t\tbs=<value>\t\t" << std::endl;
					std::cerr << "  dem_filter_freq=<value>\tdemff=<value>\t\t" << std::endl;
					return 0;
				} else {
					// Parse other arguments
					size_t pos = param.find('=');
					if (pos != std::string::npos) {
						std::string name = param.substr(0, pos);
						std::string value = param.substr(pos + 1);
						
						// Compare the parameter name to retrieve the value
						if (name == "frequency1" || name == "f1") {
							frequency1 = convertToInteger(value);
						} else if (name == "frequency_min" || name == "fmin") {
							frequency_min = convertToInteger(value);
						} else if (name == "frequency_max" || name == "fmax") {
							frequency_max = convertToInteger(value);
						} else if (name == "nb_frequency" || name == "nf") {
							nb_frequency = std::stod(value);
						} else if (name == "amplitude1" || name == "a1") {
							amplitude1 = std::stod(value);
						} else if (name == "amplitude2" || name == "a2") {
							amplitude2 = std::stod(value);
						} else if (name == "offset" || name == "off") {
							offset = std::stod(value);
						} else if (name == "phase" || name == "ph") {
							phase = std::stod(value);
						} else if (name == "noise" || name == "n") {
							noise = std::stod(value);
						} else if (name == "window_type" || name == "wt") {
							window_type = stringToWindowType(value);
						} else if (name == "decimation" || name == "d") {
							hasSetDecimation = true;
							SetDecimation(std::stoi(value));
						} else if (name == "points_per_period" || name == "ppp") {
							points_per_period = convertToInteger(value);
						} else if (name == "dem_filter_freq" || name == "demff") {
							dem_filter_freq = convertToInteger(value);
						} else if (name == "buffsize" || name == "bs") {
							SetBufferSize(convertToInteger(value));
						} else {
							std::cerr << "Invalid argument format: " << param << std::endl;
							return 1;
						}
					} else {
						std::cerr << "Invalid argument format: " << param << std::endl;
						return 1;
					}
				}
			}
		}

		// on calcul la décimation si elle n'a pas été définie
		if (!hasSetDecimation)
			SetDecimation(calculateDecimation(frequency1, points_per_period));

		/* - - - - - - - - - - - - - - - - - - - - - - - */
		/* séquence de fréquences */
		std::vector<int> freq_sequence = {100, 200, 300, 400, 600, 800, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 10000};

		/* - - - - - - - - - - - - - - - - - - - - - - - */
		/* affichage des paramètres*/
		std::cerr << "+--------------- Parameters ---------------+" << std::endl;
		std::cerr << "|                 Sampling                 |" << std::endl;
		std::cerr << "| fs         : " << std::setw(23) << SAMPLING_FREQUENCY << " Hz  |" << std::endl;
		std::cerr << "| decimation : " << std::setw(27) << DECIMATION  << " |" << std::endl;
		std::cerr << "| buffsize   : " << std::setw(21) << BUFFER_SIZE << "/16384 |" << std::endl;
		double T_per_period = 1.0 / frequency1;
		double ppp = T_per_period * SAMPLING_FREQUENCY;
		std::cerr << "| Samples per period : " << std::setw(19) << ppp << " |" << std::endl;
		std::cerr << "+------------------------------------------+" << std::endl;
		std::cerr << "|                Generator                 |" << std::endl;
		std::cerr << "| frequency1 : " << std::setw(23) << frequency1  << " Hz  |" << std::endl;
		std::cerr << "| amplitude1 : " << std::setw(23) << amplitude1  << " V   |" << std::endl;
		std::cerr << "| sequence   :                             |"    << std::endl;
		for (auto& f : freq_sequence) {
			std::cerr << "|              " << std::setw(23) << f       << " Hz  |" << std::endl;
		}
		std::cerr << "| amplitude2 : " << std::setw(23) << amplitude2  << " V   |" << std::endl;
		std::cerr << "| offset     : " << std::setw(23) << offset      << " V   |" << std::endl;
		std::cerr << "| phase      : " << std::setw(23) << phase       << " rad |" << std::endl;
		std::cerr << "| gain noise : " << std::setw(23) << noise       << "     |" << std::endl;
		std::cerr << "| window     : " << std::setw(23) << windowTypeToString(window_type) << "     |" << std::endl;
		std::cerr << "+------------------------------------------+"    << std::endl;
		std::cerr << "|               Demodulation               |"    << std::endl;
		std::cerr << "| Freq filter: " << std::setw(23) << dem_filter_freq << " Hz  |" << std::endl;
		std::cerr << "+------------------------------------------+"    << std::endl << std::endl;

		/* - - - - - - - - - - - - - - - - - - - - - - - */

		std::vector<Signal> outSig;
		std::vector<Signal> outAmp;
		std::vector<Signal> outPha;
		std::vector<Spectrum> outSpAmp;
		std::vector<Spectrum> outSpPha, outSpSig;

		Signal signal;
		Signal signal_ampli;
		Signal signal_phase;

		Spectrum FFT_sig;
		Spectrum FFT_ampli;


		WhiteNoise whiteNoise;
		whiteNoise.setGain(noise);

		/* - - - - - - - - - - - - - - - - - - - - - - - */
		/* Initialisation de la démodulation */
		Demodulator dem(dem_filter_freq, frequency1);
		dem.setup();

		std::cerr << "+---------------- Process  ----------------+" << std::endl;
		for (int freq : freq_sequence) {
			int frequency2 = frequency1+freq;
			std::cerr << "\r| frequency2 : " << std::setw(23) << frequency2 << " Hz  |" << std::flush;

			std::string fstr = std::to_string(frequency2);
			signal.setName("signal_" + fstr + "(t)");
			signal.resize(BUFFER_SIZE, 0.0);
			signal_ampli.setName("amplitude_" + fstr + "(t)");
			signal_ampli.resize(BUFFER_SIZE, 0.0);
			signal_phase.setName("phase_" + fstr + "(t)");
			signal_phase.resize(BUFFER_SIZE, 0.0);

			/* Signal arbitraire formé d'une somme de signaux sinusoïdaux */
			double t;
			for (size_t k = 0; k < BUFFER_SIZE; ++k) {
				t = (k / SAMPLING_FREQUENCY) + offset;
				signal[k] = whiteNoise.apply(amplitude1 * sin(2*M_PI*frequency1*t + phase) + amplitude2 * sin(2*M_PI*frequency2*t + phase));
			}

			/* - - - - - - - - - - - - - - - - - - - - - - - */
			/* Démodulation du signal */
			//std::cerr << "---- Demodulate ----" << std::endl;
			dem.apply(signal, signal_ampli, signal_phase);

			outSig.push_back(signal);
			outAmp.push_back(signal_ampli);
			outPha.push_back(signal_phase);
		}
		std::cerr << "\n| End of process                           |" << std::endl;
		std::cerr << "+------------------------------------------+" << std::endl << std::endl;

		/* - - - - - - - - - - - - - - - - - - - - - - - */
		/* Calcul de la durée de rising time max */
		std::cerr << "+------------ Rising time max -------------+" << std::endl;
		size_t max_high_index = 0;
		for (size_t idSig = 0; idSig < outSig.size(); idSig++) {
			size_t low_index, high_index;
			// on prend une marge en multipliant par 2 le temps en régime transitoire pour s'assurer que la dft est effectuée seulement sur le régime stable
			double rising_time = outSig[idSig].getRisingTime(low_index, high_index);
			high_index = low_index + (rising_time/T_per_period) * 5; // régime transitoire à 5 tau
			if (high_index > max_high_index) {
				max_high_index = high_index;
			}
			
		}
		std::cerr << "| Rising Time : " << std::setw(23) << max_high_index*SAMPLING_FREQUENCY << " s  |" << std::endl;
		std::cerr << "+------------------------------------------+" << std::endl << std::endl;

		/* - - - - - - - - - - - - - - - - - - - - - - - */
		/* Calcul des transformées de fourier discrètes des signaux avec BUFFER_SIZE zero padding */
		std::cerr << "+------------- Calculate FFT --------------+" << std::endl;
		for (size_t idSig = 0; idSig < outSig.size(); idSig++) {
			std::cerr << "\r| Frequency : " << std::setw(23) << freq_sequence[idSig] << " Hz   |" << std::flush;
			FFT_sig.setName("SIGNAL_" + std::to_string(freq_sequence[idSig]) + "(f)");
			FFT_ampli.setName("AMPLITUDE_" + std::to_string(freq_sequence[idSig]) + "(f)");
			outAmp[idSig].FFT(FFT_ampli, max_high_index);
			outSig[idSig].FFT(FFT_sig);
			outSpAmp.push_back(FFT_ampli);
			outSpSig.push_back(FFT_sig);
		}
		std::cerr << "\n| End of FFT                               |" << std::endl;
		std::cerr << "+------------------------------------------+" << std::endl << std::endl;

		/* - - - - - - - - - - - - - - - - - - - - - - - */
		/* Sauvgarde des signaux */


		CSVFile::setTimeToFilepath();
		std::string filename = "test_demodulation2";

		std::cerr << "+--------------- Write results ---------------+" << std::endl;
		std::cerr << "| File: " << std::setw(26) << filename << "_signals.csv |" << std::endl;
		CSVFile outFile1(filename+"_signals.csv");
		outFile1.writeSignals(outSig, true);

		std::cerr << "| File: " << std::setw(23) << filename << "_amplitudes.csv |" << std::endl;
		CSVFile outFile2(filename+"_amplitudes.csv");
		outFile2.writeSignals(outAmp, true);

		std::cerr << "| File: " << std::setw(27) << filename << "_phases.csv |" << std::endl;
		CSVFile outFile3(filename+"_phases.csv");
		outFile3.writeSignals(outPha, true);

		std::cerr << "| File: " << std::setw(22) << filename << "_FFT_signals.csv |" << std::endl;
		CSVFile outFile6(filename+"_FFT_signals.csv");
		outFile6.writeSpectrums(outSpSig, true);

		std::cerr << "| File: " << std::setw(19) << filename << "_FFT_amplitudes.csv |" << std::endl;
		CSVFile outFile4(filename+"_FFT_amplitudes.csv");
		outFile4.writeSpectrums(outSpAmp, true);
		std::cerr << "+----------------------------------------------+" << std::endl << std::endl;

		/* - - - - - - - - - - - - - - - - - - - - - - - */

		result = 0;
	} catch (const std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		result = 1;
	}
	return result;
}
