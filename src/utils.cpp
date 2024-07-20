#include "utils.hpp"

bool stringToBool(const std::string &str) {
    if (str == "true" || str == "True" || str == "TRUE" || str == "1") {
        return true;
    } else if (str == "false" || str == "False" || str == "FALSE" || str == "0") {
        return false;
    } else {
        throw std::invalid_argument("Invalid boolean value: " + str);
    }
}

std::string waveformTypeToString(WaveformType waveform) {
	switch (waveform) {
		case WaveformType::SINUS:
			return "sinus";
		case WaveformType::COSINUS:
			return "cosinus";
		case WaveformType::SQUARE:
			return "square";
		case WaveformType::TRIANGLE:
			return "triangle";
		case WaveformType::RAMP_UP:
			return "ramp_up";
		case WaveformType::RAMP_DOWN:
			return "ramp_down";
		case WaveformType::POSITIVE_DC:
			return "dc";
		case WaveformType::NEGATIVE_DC:
			return "ndc";
		default:
			return "sine";
	}
}

WaveformType stringToWaveformType(const std::string &str) {
    if (str == "sinus" || str == "SINUS") {
        return WaveformType::SINUS;
	} else if (str == "cosinus" || str == "COSINUS") {
        return WaveformType::COSINUS;
    } else if (str == "square" || str == "SQUARE") {
        return WaveformType::SQUARE;
    } else if (str == "triangle" || str == "TRIANGLE") {
        return WaveformType::TRIANGLE;
    } else if (str == "ramp_up" || str == "RAMP_UP") {
        return WaveformType::RAMP_UP;
    } else if (str == "ramp_down" || str == "RAMP_DOWN") {
        return WaveformType::RAMP_DOWN;
    } else if (str == "dc" || str == "DC" || str == "POSITIVE_DC" || str == "positive_dc") {
        return WaveformType::POSITIVE_DC;
	} else if (str == "ndc" || str == "NDC" || str == "NEGATIVE_DC" || str == "negative_dc") {
        return WaveformType::NEGATIVE_DC;
    } else {
        throw std::invalid_argument("Invalid waveform type: " + str);
    }
}

// for string delimiter
std::vector<std::string> split(std::string s, std::string delimiter) {
	size_t pos_start = 0, pos_end, delim_len = delimiter.length();
	std::string token;
	std::vector<std::string> res;

	while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
		token = s.substr (pos_start, pos_end - pos_start);
		pos_start = pos_end + delim_len;
		res.push_back (token);
	}

	res.push_back (s.substr (pos_start));
	return res;
}

int calculateDecimation(double f, int N_per_period) {
	// Calcul du nombre d'échantillons par période
	double T_per_period = 1.0 / f;
	int N = static_cast<int>(round(T_per_period * MAX_SAMPLING_FREQUENCY)); // Nombre total d'échantillons par période

	// Trouver la décimation la plus grande puissance de 2 inférieure ou égale à 16384 et ne dépassant pas N
	int decimation = 1;
	while (N > 2*N_per_period && decimation <= 16384) {
		decimation *= 2;
		N = static_cast<int>(round(T_per_period * static_cast<double>(MAX_SAMPLING_FREQUENCY)/static_cast<double>(decimation))); // Nombre total d'échantillons par période
	}

	return decimation;
}

uint32_t getTimeDelay(int decimation) {	
	/* Find optimal decimation setting */
	if (decimation < 8) {
		return 131;
	} else if(decimation < 64) {
		return 1048;
	} else if(decimation < 1024) {
		return 8388;
	} else if(decimation < 8192) {
		return 134200;
	} else if(decimation < 65536) {
		return 1000000;
	} else if(decimation == 65536) {
		return 8589000;
	} else {
		return 1000000;
	}
}

int convertToInteger(const std::string& str) {
    size_t length = str.length();
    if (length == 0) {
        throw std::invalid_argument("Empty string is not a valid number.");
    }

    double number = 0.0;
    int exponent = 0;

    // Identify and separate the numeric part and the suffix
    size_t pos = 0;
    bool hasExponent = false;
    while (pos < length && (isdigit(str[pos]) || str[pos] == '.' || str[pos] == '-' || str[pos] == '+')) {
        if (str[pos] == 'e' || str[pos] == 'E') {
            hasExponent = true;
            break;
        }
        pos++;
    }

    if (pos > 0) {
        number = std::stod(str.substr(0, pos));
    } else {
        throw std::invalid_argument("Invalid numeric part.");
    }

    if (hasExponent) {
        // Handle scientific notation
        size_t exponentPos = pos;
        pos++;
        while (pos < length && (isdigit(str[pos]) || str[pos] == '-' || str[pos] == '+')) {
            pos++;
        }
        exponent = std::stoi(str.substr(exponentPos + 1, pos - exponentPos - 1));
    } else if (pos < length) {
        // Handle suffix notation
        char suffix = str[pos];
        switch (suffix) {
            case 'M': exponent = 6; break;
            case 'k':
            case 'K': exponent = 3; break;
            case 'm': exponent = -3; break;
            case 'u': exponent = -6; break; // micro
            case 'n': exponent = -9; break;
            case 'p': exponent = -12; break;
            default: throw std::invalid_argument("Unknown suffix.");
        }
        pos++;
    }

    if (pos != length) {
        throw std::invalid_argument("Invalid format.");
    }

    // Apply the exponent
    double result = number * std::pow(10, exponent);

    // Convert to integer
    return static_cast<int>(result);
}

// Fonction pour calculer le PGCD de deux nombres
int gcd(int a, int b) {
	while (b != 0) {
		int t = b;
		b = a % b;
		a = t;
	}
	return a;
}

// Fonction pour calculer le PGCD de plusieurs nombres
int gcd_multiple(const std::vector<int>& numbers) {
	// Initialiser le GCD avec le premier élément du vecteur
	int result = numbers[0];
	// Utiliser std::accumulate pour appliquer gcd sur tout le vecteur
	result = std::accumulate(numbers.begin(), numbers.end(), result, gcd);
	return result;
}

// Fonction pour générer une période complète du signal
std::vector<double> generate_waveform_with_n_sinus(size_t buffer_size, std::vector<int> frequencies, std::vector<double> amplitudes, int &fundamental_frequency) {
	// Trouver la période fondamentale
	fundamental_frequency = gcd_multiple(frequencies);
	double period = 1.0 / static_cast<double>(fundamental_frequency);

	// Initialiser le buffer
	std::vector<double> waveform(buffer_size, 0.0f);

	// Générer le signal pour une période complète
	double t, dt = period / BUFFER_SIZE;
	for (size_t i = 0; i < buffer_size; ++i) {
		t = static_cast<double>(i) * dt;
		for (size_t j = 0; j < amplitudes.size(); ++j) {
			waveform[i] += amplitudes[j] * sin(2 * M_PI * frequencies[j] * t);
		}
	}

	// Normaliser le signal si nécessaire
	double max_amplitude = waveform[0];
	for (size_t i = 0; i < buffer_size; ++i) {
		if (std::abs(waveform[i]) > max_amplitude) {
			max_amplitude = std::abs(waveform[i]);
		}
	}
	if (max_amplitude > 1.0) {
		for (size_t i = 0; i < buffer_size; ++i) {
			waveform[i] /= max_amplitude;
		}
	}

	return waveform;
}