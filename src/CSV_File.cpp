#include <algorithm>
#include "CSV_File.hpp"

bool CSV_File::read(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }

    m_names.clear();
    m_signals.clear();

    std::string line;
    {
        std::getline(file, line);
        std::istringstream iss(line);
        iss >> sample_frequency;
    }

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string name;
        size_t size;
        iss >> name >> size;
        std::shared_ptr<Signal> signal = std::make_shared<Signal>(size, sample_frequency);
        for (size_t i = 0; i < size; i++) {
            double real, imag;
            iss >> real >> imag;
            (*signal)[i] = complexd(real, imag);
        }
        m_names.push_back(name);
        m_signals.push_back(signal);
    }

    file.close();
    return true;
}

bool CSV_File::write(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }

    file << sample_frequency << std::endl;
    for (size_t i = 0; i < m_signals.size(); i++) {
        file << m_names[i] << " " << m_signals[i]->size() << " ";
        for (size_t j = 0; j < m_signals[i]->size(); j++) {
            file << (*m_signals[i])[j].real() << " " << (*m_signals[i])[j].imag();
            if (j < m_signals[i]->size() -1) file << " ";
        }
        if (i < m_signals.size()-1) file << std::endl;
    }

    file.close();
    return true;
}

size_t CSV_File::count() const {
    return m_signals.size();
}

void CSV_File::add(const std::string& name, std::shared_ptr<Signal> signal) {
    m_names.push_back(name);
    m_signals.push_back(signal);
}

std::shared_ptr<Signal> CSV_File::create(const std::string& name, size_t size, complexd defaultValue) {
    m_names.push_back(name);
    auto new_signal = std::make_shared<Signal>(size, sample_frequency, defaultValue);
    m_signals.push_back(new_signal);
    return new_signal;
}

std::shared_ptr<Signal> CSV_File::get(size_t row) {
    if (row >= m_signals.size()) {
        throw std::out_of_range("Row index out of range");
    }
    return m_signals[row];
}

std::shared_ptr<Signal> CSV_File::get(const std::string &name) {
    auto it = std::find(m_names.begin(), m_names.end(), name);
    if (it == m_names.end()) {
        throw std::invalid_argument("Name dosen't exist");
    }
    size_t index = std::distance(m_names.begin(), it);
    return m_signals[index];
}