#ifndef CSVFILE_HPP
#define CSVFILE_HPP

#include <fstream>
#include <vector>
#include <complex>
#include <sstream>
#include <string>
#include <iostream>
#include <filesystem>
#include "Signal.hpp"
#include "Spectrum.hpp"



// Method to generate a frequency axis
std::vector<double> frequency_axis(size_t n, double d = 1.0);

class CSVFile {
public:
	static std::string DEFAULT_FILEPATH;
	static std::string FILEPATH;

	static void setFilePath(const std::string &path) {
		FILEPATH = path;
		std::filesystem::create_directories(FILEPATH);
	}

	static void setTimeToFilepath(bool withDate = true, bool withTime = true) {
		// Get the current date and time
		time_t now = time(0);
		tm *ltm = localtime(&now);
		std::stringstream ss;
		ss << DEFAULT_FILEPATH;
		if (withDate) {
			ss << (1900 + ltm->tm_year) << "-" << (1 + ltm->tm_mon) << "-" << ltm->tm_mday;
		}
		if (withTime) {
			ss << "_" << ltm->tm_hour << ":" << ltm->tm_min << ":" << ltm->tm_sec << "/";
		} else if (withDate) {
			ss << "/";
		}
		FILEPATH = ss.str();

		// créer le répertoire s'il n'existe pas
		std::filesystem::create_directories(FILEPATH);
	}

	CSVFile(const std::string &filename);

	// Method to read signals from a CSV file
	std::vector<Signal> readSignals();

	// Method to write signals to a CSV file
	void writeSignals(const std::vector<Signal> &signals, bool axis = true);

	// Method to write signal to a CSV file
	void writeSignal(const Signal &signal, bool axis = true);

	// Method to append signals to a CSV file
	void writeSignalsToEnd(const std::vector<Signal> &signals);
	
	// Method to add a signal to a CSV file
	void writeSignalToEnd(const Signal &signal);

	// Method to read spectrums from a CSV file
	std::vector<Spectrum> readSpectrums();

	// Method to write spectrums to a CSV file
	void writeSpectrums(const std::vector<Spectrum> &spectrums, bool axis = true, bool withNegativeFrequencies = false);

	// Method to append spectrums to a CSV file
	void writeSpectrumsToEnds(const std::vector<Spectrum> &spectrums, bool withNegativeFrequencies = false);

	// Method to add a spectrum to a CSV file
	void writeSpectrumToEnd(const Spectrum &spectrum, bool withNegativeFrequencies = false);

private:
	// Utility function to parse a complex number from a string
	complexd parseComplex(const std::string &str);

	void ensureSameSize(const std::vector<Signal>& signals);
    void ensureSameSize(const std::vector<Spectrum>& spectrums);

	// Overload for writing complex numbers to a stream
	friend std::ostream &operator<<(std::ostream &os, const complexd &c) {
		os << c.real() << (c.imag() >= 0 ? "+" : "") << c.imag() << "j";
		return os;
	}

	std::string _filename;
    std::fstream _fileStream;
};

#endif // CSVFILE_HPP