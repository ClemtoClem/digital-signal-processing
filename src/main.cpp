#include <iostream>
#include "globals.hpp"
#include "Signal.hpp"
#include "Spectrum.hpp"
#include "CSVFile.hpp"
#include "Filter.hpp"
#include "Noise.hpp"
#include "Demodulator.hpp"
#include "tests.hpp"

int main(int argc, char *argv[]){
	int res = 0;

	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " <name/module_name> [optional arguments]" << std::endl;
		std::cerr << "Write 'help' to see available tests or modules." << std::endl;
		return 1;
	}

	// Parser les arguments
	std::vector<std::string> args(argv + 1, argv + argc);
	std::string name = args[0];

	// Supprimer le premier argument (name)
	args.erase(args.begin());

	if (name == "test_simple") {
		res |= test_simple(args);
	} else if (name == "test_spectrum") {
		res |= test_spectrum(args);
	} else if (name == "test_demodulation") {
		res |= test_demodulation(args);
	} else if (name == "test_demodulation2") {
		res |= test_demodulation2(args);
	} else if (name == "help") {
		std::cout << "Available tests:" << std::endl;
		std::cout << "\ttest_simple <optional arguments>" << std::endl;
		std::cout << "\ttest_spectrum <optional arguments>" << std::endl;
		std::cout << "\ttest_demodulation <optional arguments>" << std::endl;
		std::cout << "\ttest_demodulation2 <optional arguments>" << std::endl;
	} else {
		std::cerr << "Unknown test name: " << name << std::endl;
		std::cerr << "Write 'help' to see available tests." << std::endl;
		return 1;
	}

	return res;
}
