#include "CSVFile.hpp"
#include "globals.hpp"

int calculateDecimation(double fs, double f, int N_per_period) {
    // Calcul du nombre d'échantillons par période
    double T_per_period = 1.0 / f;
    int N = static_cast<int>(round(T_per_period * fs)); // Nombre total d'échantillons par période

    // Trouver la décimation la plus grande puissance de 2 inférieure ou égale à 16384 et ne dépassant pas N
    int decimation = 1;
    while (decimation < N_per_period && decimation <= 16384) {
        decimation *= 2;
    }

    return decimation;
}

CSVFile::CSVFile(const std::string &filename) : _filename(filename) {}

std::vector<Signal> CSVFile::readSignals() {
    _fileStream.open(_filename, std::ios::in);
    if (!_fileStream.is_open()) {
        throw std::ios_base::failure("Failed to open file");
    }

    std::vector<Signal> signals;
    std::string line;

    while (std::getline(_fileStream, line))  {
        std::stringstream _fileStream(line);
        std::string name;
        std::getline(_fileStream, name, ','); // Read the signal name
        std::string token;
        std::vector<double> values;
        while (std::getline(_fileStream, token, ',')) {
            values.push_back(std::atol(token.c_str()));
        }
        signals.emplace_back(values, name);
    }
    _fileStream.close();

    return signals;
}

void CSVFile::writeSignals(const std::vector<Signal> &signals, bool axis) {
    if (signals.empty()) {
        throw std::invalid_argument("No signals provided");
    }

    _fileStream.open(_filename, std::ios::out | std::ios::trunc);
    if (!_fileStream.is_open()) {
        throw std::ios_base::failure("Failed to open file");
    }

    size_t nbSignals = signals.size();

    // Write the header
    if (axis) {
        _fileStream << "time,";
    }
    for (size_t s = 0; s < nbSignals; s++) {
        _fileStream << signals[s].getName() << ((s < nbSignals-1)?",":""); // Write the signal name
    }
    _fileStream << "\n";

    // Find the maximum size of all signals
    size_t maxSize = 0;
    for (const auto &signal : signals) {
        if (signal.size() > maxSize) {
            maxSize = signal.size();
        }
    }

    // Write the time and signal values
    double axis_value;
    for (size_t i = 0; i < maxSize; ++i) {
        if (axis) { // time
            axis_value = i / SAMPLING_FREQUENCY;
            _fileStream << axis_value;
            if (nbSignals != 0) {
                _fileStream << ",";
            }
        }
        for (size_t s = 0; s < nbSignals; s++) {
            if (i < signals[s].size()) {
                _fileStream << signals[s][i];
            } else {
                _fileStream << "0"; // Use 0 for mi_fileStreaming values
            }
            if (s < nbSignals-1) {
                _fileStream << ",";
            }
        }
        _fileStream << "\n";
    }
    _fileStream.close();
}

void CSVFile::appendSignals(const std::vector<Signal> &signals) {
    if (signals.empty()) {
        throw std::invalid_argument("No signals provided");
    }

    std::fstream inFile(_filename, std::ios::in);
    if (!inFile.is_open()) {
        throw std::ios_base::failure("Failed to open file for reading");
    }

    std::vector<std::string> lines;
    std::string line;

    // Read all lines from the file
    while (std::getline(inFile, line)) {
        lines.push_back(line);
    }
    inFile.close();

    // Re-open the file in write mode to overwrite it
    _fileStream.open(_filename, std::ios::out);
    if (!_fileStream.is_open()) {
        throw std::ios_base::failure("Failed to open file for writing");
    }

    // Write the updated header
    line = lines[0];
    for (const auto &signal : signals) {
        line += "," + signal.getName();
    }
    _fileStream << line << "\n";

    // Write the updated data lines
    for (size_t i = 1; i < lines.size(); ++i) {
        line = lines[i];
        size_t pos = 0;
        for (const auto &signal : signals) {
            pos = line.find(",", pos);
            if (pos != std::string::npos) {
                line.insert(pos + 1, std::to_string(signal[i - 1]) + ",");
            }
        }
        _fileStream << line << "\n";
    }

    _fileStream.close();
}


std::vector<Spectrum> CSVFile::readSpectrums() {
    _fileStream.open(_filename, std::ios::in);
    if (!_fileStream.is_open()) {
        throw std::ios_base::failure("Failed to open file");
    }

    std::vector<Spectrum> spectrums;
    std::string line;

    while (std::getline(_fileStream, line))  {
        std::stringstream _fileStream(line);
        std::string name;
        std::getline(_fileStream, name, ','); // Read the signal name
        std::string token;
        std::vector<complexd> values;
        while (std::getline(_fileStream, token, ',')) {
            values.push_back(parseComplex(token.c_str()));
        }
        spectrums.emplace_back(values, name); // samplingFrequency is set to 0 for this example
    }
    _fileStream.close();

    return spectrums;
}

void CSVFile::writeSpectrums(const std::vector<Spectrum> &spectrums, bool axis, double freqOffset) {
    if (spectrums.empty()) {
        throw std::invalid_argument("No spectrum provided");
    }

    _fileStream.open(_filename, std::ios::out | std::ios::trunc);
    if (!_fileStream.is_open()) {
        throw std::ios_base::failure("Failed to open file");
    }

    size_t nbSpectrums = spectrums.size();

    // Find the maximum size of all spectrums
    size_t N = spectrums[0].size();
    for (const auto &spectrum : spectrums) {
        if (spectrum.size() < N) {
            N = spectrum.size();
        }
    }

    // Write the header
    if (axis) {
        _fileStream << "freq,";
    }

    for (size_t s = 0; s < nbSpectrums; s++) {
        _fileStream << spectrums[s].getName();
        if (s < nbSpectrums - 1) {
            _fileStream << ",";
        }
    }
    _fileStream << "\n";

    // Write the frequency and spectrums values
    double axis_value;
    for (size_t k = 0; k < N; k++) {
        if (axis) { // frequency
            axis_value = ((static_cast<double>(k) * SAMPLING_FREQUENCY) / static_cast<double>(N)) -freqOffset;
            _fileStream << axis_value << ",";
        }

        for (size_t s = 0; s < nbSpectrums; s++) {
            if (k < spectrums[s].size()) {
                const complexd &c = spectrums[s][k];
                _fileStream << c.real();
                if (c.imag() != 0) {
                    _fileStream << (c.imag() >= 0 ? "+" : "") << c.imag() << "j";
                }
            } else {
                _fileStream << "0"; // Use 0 for missing values
            }
            if (s < nbSpectrums - 1) {
                _fileStream << ",";
            }
        }
        _fileStream << "\n";
    }
    _fileStream.close();
}

complexd CSVFile::parseComplex(const std::string &str) {
    double real = 0.0, imag = 0.0;
    size_t pos = str.find_first_of("+-", 1);
    if (pos != std::string::npos) {
        real = std::stod(str.substr(0, pos));
        imag = std::stod(str.substr(pos));
    } else {
        real = std::stod(str);
    }
    return complexd(real, imag);
}

void CSVFile::appendSpectrums(const std::vector<Spectrum> &spectrums) {
    _fileStream.open(_filename, std::ios::in);
    if (!_fileStream.is_open()) {
        throw std::ios_base::failure("Failed to open file");
    }

    if (spectrums.empty()) {
        throw std::invalid_argument("No spectrums provided");
    }

    // Read the existing content
    std::vector<std::string> lines;
    std::string line;
    size_t maxLines = 0;

    // Read existing lines
    while (std::getline(_fileStream, line)) {
        lines.push_back(line);
        if (lines.size() > 1) { // Ignore the header line for maxLines calculation
            ++maxLines;
        }
    }

    // Close and reopen the file to write
    _fileStream.close();
    _fileStream.open(_filename, std::ios::app);
    if (!_fileStream.is_open()) {
        throw std::ios_base::failure("Failed to reopen file");
    }

    // Proce_fileStream the header
    std::stringstream _fileStreamHeader(lines[0]);
    std::string header;
    std::getline(_fileStreamHeader, header, '\n');
    std::string newHeader = header;

    for (const auto &spectrum : spectrums) {
        newHeader += "," + spectrum.getName();
    }
    _fileStream << newHeader << "\n";

    size_t spectrumMaxSize = 0;
    for (const auto &spectrum : spectrums) {
        if (spectrum.size() > spectrumMaxSize) {
            spectrumMaxSize = spectrum.size();
        }
    }

    // Proce_fileStream each line, appending new spectrum data or filling with zeros
    for (size_t i = 1; i <= maxLines; ++i) {
        std::stringstream _fileStreamLine(lines[i]);
        std::string lineData;
        std::getline(_fileStreamLine, lineData, '\n');
        _fileStream << lineData;

        for (const auto &spectrum : spectrums) {
            _fileStream << ",";
            if (i - 1 < spectrum.size()) {
                const complexd &c = spectrum[i - 1];
                _fileStream << c.real();
                if (c.imag() != 0) {
                    _fileStream << (c.imag() >= 0 ? "+" : "") << c.imag() << "j";
                }
            } else {
                _fileStream << "0";
            }
        }
        _fileStream << "\n";
    }

    // Add any remaining data from the new spectrums that exceed existing lines
    for (size_t i = maxLines; i < spectrumMaxSize; ++i) {
        double freq = i * SAMPLING_FREQUENCY / spectrumMaxSize;
        _fileStream << freq;

        for (const auto &spectrum : spectrums) {
            _fileStream << ",";
            if (i < spectrum.size()) {
                const complexd &c = spectrum[i];
                _fileStream << c.real();
                if (c.imag() != 0) {
                    _fileStream << c;
                }
            } else {
                _fileStream << "0";
            }
        }
        _fileStream << "\n";
    }

    _fileStream.close();
}
