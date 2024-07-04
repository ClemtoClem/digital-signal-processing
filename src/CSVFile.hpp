#ifndef CSVFILE_HPP
#define CSVFILE_HPP

#include <fstream>
#include <vector>
#include <complex>
#include <sstream>
#include <string>
#include <iostream>
#include "Signal.hpp"
#include "Spectrum.hpp"

class CSVFile {
private:
    std::string _filename;
    std::fstream _fileStream;
public:
    CSVFile(const std::string &filename);

    // Method to read signals from a CSV file
    std::vector<Signal> readSignals();

    // Method to write signals to a CSV file
    void writeSignals(const std::vector<Signal> &signals, bool axis = true);

    // Method to append signals to a CSV file
    void appendSignals(const std::vector<Signal> &signals);

    // Method to read spectrums from a CSV file
    std::vector<Spectrum> readSpectrums();

    // Method to write spectrums to a CSV file
    void writeSpectrums(const std::vector<Spectrum> &spectrums, bool axis = true, double freqOffset = 0);

    // Method to append spectrums to a CSV file
    void appendSpectrums(const std::vector<Spectrum> &spectrums);

private:
    // Utility function to parse a complex number from a string
    complexd parseComplex(const std::string &str);

    // Overload for writing complex numbers to a stream
    friend std::ostream &operator<<(std::ostream &os, const complexd &c) {
        os << c.real() << (c.imag() >= 0 ? "+" : "") << c.imag() << "j";
        return os;
    }
};

#endif // CSVFILE_HPP