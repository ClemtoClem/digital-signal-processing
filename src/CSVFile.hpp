#ifndef __CSVFILE_HPP
#define __CSVFILE_HPP

#include <fstream>
#include <vector>
#include <complex>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "Signal.hpp"

enum class Axis {
    No,
    Time,
    Freq
};

class CSVFile {
public:
    CSVFile(const std::string &filename);

    // Method to read signals from a CSV file
    std::vector<Signal> readSignals(double samplingFrequency);

    // Method to write signals to a CSV file
    void writeSignals(const std::vector<Signal> &signals, Axis axis = Axis::Time);

    // Method to append a single signal to an existing CSV file
    void appendSignal(const std::vector<Signal> &signals);

private:
    // Utility function to parse a complex number from a string
    std::complex<double> parseComplex(const std::string &str);

    // Overload for writing complex numbers to a stream
    friend std::ostream &operator<<(std::ostream &os, const std::complex<double> &c) {
        os << c.real() << (c.imag() >= 0 ? "+" : "") << c.imag() << "j";
        return os;
    }

    std::string mFilename;
    std::fstream mFileStream;
};

#endif // __CSVFILE_HPP