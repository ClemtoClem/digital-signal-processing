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
        std::vector<std::complex<double>> values;
        while (std::getline(ss, token, ',')) {
            values.push_back(parseComplex(token));
        }
        signals.emplace_back(values, samplingFrequency, name); // sample_frequency is set to 0 for this example
    }
    mFileStream.close();

    return signals;
}

void CSVFile::writeSignals(const std::vector<Signal> &signals, bool time_or_freq_axis) {
    if (signals.empty()) {
        throw std::invalid_argument("No signals provided");
    }

    mFileStream.open(mFilename, std::ios::out);
    if (!mFileStream.is_open()) {
        throw std::ios_base::failure("Failed to open file");
    }

    // Write the header
    mFileStream << ((time_or_freq_axis)? "freq":"time");
    for (auto &signal : signals) {
        mFileStream << "," << signal.getName(); // Write the signal name
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
        if (!time_or_freq_axis) // time
            axis_value = i / signals[0].getSamplingFrequency(); 
        else // frequency
            axis_value = i * signals[0].getSamplingFrequency() / maxSize;
        mFileStream << axis_value;
        for (const auto &signal : signals) {
            mFileStream << ",";
            if (i < signal.size()) {
                const std::complex<double> &c = signal[i];
                mFileStream << c.real();
                if (c.imag()!= 0) {
                    mFileStream << (c.imag() >= 0 ? "+" : "") << c.imag() << "j";
                }
            } else {
                mFileStream << "0"; // Use 0 for missing values
            }
        }
        mFileStream << "\n";
    }
    mFileStream.close();
}


void CSVFile::appendSignal(const std::vector<Signal> &signals) {
    mFileStream.open(mFilename, std::ios::in);
    if (!mFileStream.is_open()) {
        throw std::ios_base::failure("Failed to open file");
    }

    if (signals.empty()) {
        throw std::invalid_argument("No signals provided");
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

    for (const auto &signal : signals) {
        newHeader += "," + signal.getName();
    }
    mFileStream << newHeader << "\n";

    size_t signalMaxSize = 0;
    for (const auto &signal : signals) {
        if (signal.size() > signalMaxSize) {
            signalMaxSize = signal.size();
        }
    }

    // Process each line, appending new signal data or filling with zeros
    for (size_t i = 1; i <= maxLines; ++i) {
        std::stringstream ssLine(lines[i]);
        std::string lineData;
        std::getline(ssLine, lineData, '\n');
        mFileStream << lineData;

        for (const auto &signal : signals) {
            mFileStream << ",";
            if (i - 1 < signal.size()) {
                mFileStream << signal[i - 1];
            } else {
                mFileStream << "0";
            }
        }
        mFileStream << "\n";
    }

    // Add any remaining data from the new signals that exceed existing lines
    for (size_t i = maxLines; i < signalMaxSize; ++i) {
        double time = i / signals[0].getSamplingFrequency();
        mFileStream << time;

        for (const auto &signal : signals) {
            mFileStream << ",";
            if (i < signal.size()) {
                mFileStream << signal[i];
            } else {
                mFileStream << "0";
            }
        }
        mFileStream << "\n";
    }

    mFileStream.close();
}

std::complex<double> CSVFile::parseComplex(const std::string &str) {
    double real = 0.0, imag = 0.0;
    size_t pos = str.find_first_of("+-", 1);
    if (pos != std::string::npos) {
        real = std::stod(str.substr(0, pos));
        imag = std::stod(str.substr(pos));
    } else {
        real = std::stod(str);
    }
    return std::complex<double>(real, imag);
}