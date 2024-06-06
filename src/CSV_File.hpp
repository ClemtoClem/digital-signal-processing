#ifndef __CSV_FILE_HPP
#define __CSV_FILE_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <memory>
#include <stdexcept>
#include "Signal.hpp"

class CSV_File {
public:
    CSV_File(int sample_frequency) : sample_frequency(sample_frequency) {}

    // Méthode pour lire les signaux dans un fichier csv, retourne true si réussie
    bool read(const std::string& filename);

    // Méthode écrie les signaux dans un fichie csv, retourne true si réussie
    bool write(const std::string& filename);

    // Méthode pour renvoyer le nombre de signaux
    size_t count() const;

    // Méthode pour ajouter un signal avec son nom
    void add(const std::string& name, std::shared_ptr<Signal> signal);

    // Méthode pour créer un nouveau signal avec son nom
    std::shared_ptr<Signal> create(const std::string& name, size_t size, complexd defaultValue = 0);

    // Méthode pour lire le signal à une ligne donnée
    std::shared_ptr<Signal> get(size_t row);

    // Méthode pour lire un signal à partir de son nom, retourne nullptr si ne nom n'existe pas
    std::shared_ptr<Signal> get(const std::string &name);

private:
    int sample_frequency;
    std::vector<std::string> m_names;
    std::vector<std::shared_ptr<Signal>> m_signals;
};



#endif // __CSV_FILE_HPP