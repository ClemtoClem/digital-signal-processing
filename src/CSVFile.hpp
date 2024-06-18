#ifndef __CSVFILE_HPP
#define __CSVFILE_HPP

#include <fstream>
#include <vector>
#include <complex>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "globals.hpp"
#include "Signal.hpp"
#include "Spectrum.hpp"

class CSVFile {
public:
    CSVFile(const std::string &filename);

    // Method to read signals from a CSV file
    std::vector<Signal> readSignals(double samplingFrequency = SAMPLING_FREQUENCY);

    // Method to write signals to a CSV file
    void writeSignals(const std::vector<Signal> &signals, bool axis = true);

    // Method to read spectrums from a CSV file
    std::vector<Spectrum> readSpectrums(double samplingFrequency = SAMPLING_FREQUENCY);

    // Method to write spectrums to a CSV file
    void writeSpectrums(const std::vector<Spectrum> &spectrums, bool axis = true);

private:
    // Utility function to parse a complex number from a string
    complexd parseComplex(const std::string &str);

    // Overload for writing complex numbers to a stream
    friend std::ostream &operator<<(std::ostream &os, const complexd &c) {
        os << c.real() << (c.imag() >= 0 ? "+" : "") << c.imag() << "j";
        return os;
    }

    std::string mFilename;
    std::fstream mFileStream;
};

#endif // __CSVFILE_HPP