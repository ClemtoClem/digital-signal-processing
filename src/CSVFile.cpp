#include "CSVFile.hpp"

CSVFile::CSVFile(const std::string &filename) : mFilename(filename) {}

std::vector<Signal> CSVFile::readSignals(double samplingFrequency) {
    mFileStream.open(mFilename, std::ios::in);
    if (!mFileStream.is_open()) {
        throw std::ios_base::failure("Failed to open file");
    }

    std::vector<Signal> signals;
    std::string line;

    while (std::getline(mFileStream, line))  {
        std::stringstream ss(line);
        std::string name;
        std::getline(ss, name, ','); // Read the signal name
        std::string token;
        std::vector<double> values;
        while (std::getline(ss, token, ',')) {
            values.push_back(std::atol(token.c_str()));
        }
        signals.emplace_back(values, samplingFrequency, name); // samplingFrequency is set to 0 for this example
    }
    mFileStream.close();

    return signals;
}

void CSVFile::writeSignals(const std::vector<Signal> &signals, bool axis) {
    if (signals.empty()) {
        throw std::invalid_argument("No signals provided");
    }

    mFileStream.open(mFilename, std::ios::out);
    if (!mFileStream.is_open()) {
        throw std::ios_base::failure("Failed to open file");
    }

    size_t nbSignals = signals.size();

    // Write the header
    if (axis) mFileStream << "time,";
    
    for (size_t s = 0; s < nbSignals; s++) {
        mFileStream << signals[s].getName() << ((s < nbSignals-1)?",":""); // Write the signal name
    }
    mFileStream << "\n";

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
            axis_value = i / signals[0].getSamplingFrequency(); 
            mFileStream << axis_value;
        }
        if (nbSignals != 0) mFileStream << ",";
        for (size_t s = 0; s < nbSignals; s++) {
            if (i < signals[s].size()) {
                mFileStream << signals[s][i];
            } else {
                mFileStream << "0"; // Use 0 for missing values
            }
            if (s < nbSignals-1) mFileStream << ",";
        }
        mFileStream << "\n";
    }
    mFileStream.close();
}

std::vector<Spectrum> CSVFile::readSpectrums(double samplingFrequency) {
    mFileStream.open(mFilename, std::ios::in);
    if (!mFileStream.is_open()) {
        throw std::ios_base::failure("Failed to open file");
    }

    std::vector<Spectrum> spectrums;
    std::string line;

    while (std::getline(mFileStream, line))  {
        std::stringstream ss(line);
        std::string name;
        std::getline(ss, name, ','); // Read the signal name
        std::string token;
        std::vector<complexd> values;
        while (std::getline(ss, token, ',')) {
            values.push_back(parseComplex(token.c_str()));
        }
        spectrums.emplace_back(values, samplingFrequency, name); // samplingFrequency is set to 0 for this example
    }
    mFileStream.close();

    return spectrums;
}

void CSVFile::writeSpectrums(const std::vector<Spectrum> &spectrums, bool axis) {
    if (spectrums.empty()) {
        throw std::invalid_argument("No spectrum provided");
    }

    mFileStream.open(mFilename, std::ios::out);
    if (!mFileStream.is_open()) {
        throw std::ios_base::failure("Failed to open file");
    }

    size_t nbSpectrums = spectrums.size();

    // Write the header
    if (axis) mFileStream << "freq,";
    
    for (size_t s = 0; s < nbSpectrums; s++) {
        mFileStream << spectrums[s].getName() << ((s<nbSpectrums-1)?",":""); // Write the spectrum name
    }
    mFileStream << "\n";

    // Find the maximum size of all spectrums
    size_t maxSize = 0;
    for (const auto &spectrum : spectrums) {
        if (spectrum.size() > maxSize) {
            maxSize = spectrum.size();
        }
    }

    // Write the time and spectrums values
    double axis_value;
    for (size_t i = 0; i < maxSize; ++i) {
        if (axis) { // frequency
            axis_value = i * spectrums[0].getSamplingFrequency() / maxSize;
            mFileStream << axis_value;
        }
        if (nbSpectrums != 0) mFileStream << ",";
        for (size_t s = 0; s < nbSpectrums; s++) {
            if (i < spectrums[s].size()) {
                const complexd &c = spectrums[s][i];
                mFileStream << c.real();
                if (c.imag()!= 0) {
                    mFileStream << (c.imag() >= 0 ? "+" : "") << c.imag() << "j";
                }
            } else {
                mFileStream << "0"; // Use 0 for missing values
            }
            if (s < nbSpectrums-1) mFileStream << ",";
        }
        mFileStream << "\n";
    }
    mFileStream.close();
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