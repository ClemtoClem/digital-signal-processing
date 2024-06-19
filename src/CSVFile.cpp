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
            if (nbSignals != 0) mFileStream << ",";
        }
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

void CSVFile::appendSignals(const std::vector<Signal> &signals) {
    if (signals.empty()) {
        throw std::invalid_argument("No signals provided");
    }

    std::fstream inFile(mFilename, std::ios::in);
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
    mFileStream.open(mFilename, std::ios::out);
    if (!mFileStream.is_open()) {
        throw std::ios_base::failure("Failed to open file for writing");
    }

    // Write the updated header
    line = lines[0];
    for (const auto &signal : signals) {
        line += "," + signal.getName();
    }
    mFileStream << line << "\n";

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
        mFileStream << line << "\n";
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
            if (nbSpectrums != 0) mFileStream << ",";
        }
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

void CSVFile::appendSpectrums(const std::vector<Spectrum> &spectrums) {
    mFileStream.open(mFilename, std::ios::in);
    if (!mFileStream.is_open()) {
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
    while (std::getline(mFileStream, line)) {
        lines.push_back(line);
        if (lines.size() > 1) { // Ignore the header line for maxLines calculation
            ++maxLines;
        }
    }

    // Close and reopen the file to write
    mFileStream.close();
    mFileStream.open(mFilename, std::ios::app);
    if (!mFileStream.is_open()) {
        throw std::ios_base::failure("Failed to reopen file");
    }

    // Process the header
    std::stringstream ssHeader(lines[0]);
    std::string header;
    std::getline(ssHeader, header, '\n');
    std::string newHeader = header;

    for (const auto &spectrum : spectrums) {
        newHeader += "," + spectrum.getName();
    }
    mFileStream << newHeader << "\n";

    size_t spectrumMaxSize = 0;
    for (const auto &spectrum : spectrums) {
        if (spectrum.size() > spectrumMaxSize) {
            spectrumMaxSize = spectrum.size();
        }
    }

    // Process each line, appending new spectrum data or filling with zeros
    for (size_t i = 1; i <= maxLines; ++i) {
        std::stringstream ssLine(lines[i]);
        std::string lineData;
        std::getline(ssLine, lineData, '\n');
        mFileStream << lineData;

        for (const auto &spectrum : spectrums) {
            mFileStream << ",";
            if (i - 1 < spectrum.size()) {
                const complexd &c = spectrum[i - 1];
                mFileStream << c.real();
                if (c.imag() != 0) {
                    mFileStream << (c.imag() >= 0 ? "+" : "") << c.imag() << "j";
                }
            } else {
                mFileStream << "0";
            }
        }
        mFileStream << "\n";
    }

    // Add any remaining data from the new spectrums that exceed existing lines
    for (size_t i = maxLines; i < spectrumMaxSize; ++i) {
        double freq = i * spectrums[0].getSamplingFrequency() / spectrumMaxSize;
        mFileStream << freq;

        for (const auto &spectrum : spectrums) {
            mFileStream << ",";
            if (i < spectrum.size()) {
                const complexd &c = spectrum[i];
                mFileStream << c.real();
                if (c.imag() != 0) {
                    mFileStream << c;
                }
            } else {
                mFileStream << "0";
            }
        }
        mFileStream << "\n";
    }

    mFileStream.close();
}
