#ifndef __NOISE_HPP
#define __NOISE_HPP

#include <vector>
#include <random>
#include <cmath>

class Signal; // Assurez-vous que la classe Signal est déjà définie

class Noise {
protected:
    double gain;
public:
    Noise() : gain(1.0) {}

    virtual ~Noise() = default;

    void setGain(double gain) {
        this->gain = gain;
    }

    virtual void compute() = 0; // Calculer les coefficients
    virtual Signal process(const Signal &input) = 0;
    virtual double process(double input) = 0;
};

class WhiteNoise : public Noise {
private:
    std::default_random_engine generator;
    std::uniform_real_distribution<double> uniformDistribution;
public:
    WhiteNoise();
    void compute() override;
    Signal process(const Signal &input) override;
    double process(double input) override;
};

class PinkNoise : public Noise {
private:
    double gain;
    std::vector<double> b; // Numerator coefficients
    std::vector<double> a; // Denominator coefficients
    std::vector<double> x_hist; // Input history
    std::vector<double> y_hist; // Output history
public:
    PinkNoise();
    void compute() override;
    Signal process(const Signal &input) override;
    double process(double input) override;
};

class BrownNoise : public Noise {
private:
    double previous;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution;
public:
    BrownNoise();
    void compute() override;
    Signal process(const Signal &input) override;
    double process(double input) override;
};

#endif // __NOISE_HPP