#include "Spectrum.hpp"
#include "Signal.hpp"

Spectrum::Spectrum(const std::string &name) : std::vector<complexd>(BUFFER_SIZE, 0.0), mName(name) {}

Spectrum::Spectrum(size_t size, const std::string &name) : std::vector<complexd>(size), mName(name) {}

Spectrum::Spectrum(const std::vector<complexd> &values, const std::string &name) : std::vector<complexd>(values), mName(name) {}

Spectrum::Spectrum(const Spectrum &other) : std::vector<complexd>(other),mName(other.mName) {}

Spectrum &Spectrum::operator=(const Spectrum &other) {
    if (this != &other) {
        std::vector<complexd>::operator=(other);
    }
    return *this;
}

/* ------------------------------- */

const std::string &Spectrum::getName() const {
    return mName;
}

void Spectrum::setName(const std::string &name) {
    mName = name;
}

bool Spectrum::hasInfitityValue() const {
    for (size_t i = 0; i < size(); i++) {
        if (((*this)[i].real() == INFINITY) || ((*this)[i].real() == -INFINITY)) {
            return true;
        }
    }
    return false;
}

/* ------------------------------- */

Spectrum &Spectrum::addValue(complexd value) {
    std::vector<complexd>::push_back(value);
    return *this;
}

Spectrum &Spectrum::fill(complexd value)
{
    for (size_t i = 0; i < size(); i++) {
        (*this)[i] = value;
    }
    return *this;
}

/* ------------------------------- */

Spectrum Spectrum::operator+(const Spectrum &other) const {
    Spectrum output(other.size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = (*this)[i] + other[i];
    }
    return output;
}

Spectrum Spectrum::operator-(const Spectrum &other) const {
    Spectrum output(other.size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = (*this)[i] - other[i];
    }
    return output;
}

Spectrum Spectrum::operator*(const Spectrum &other) const {
    Spectrum output(other.size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = (*this)[i] * other[i];
    }
    return output;
}

Spectrum Spectrum::operator/(const Spectrum &other) const {
    Spectrum output(other.size());
    for (size_t i = 0; i < size(); i++) {
        if (std::abs(other[i]) == 0) output[i] = {INFINITY, INFINITY};
        else output[i] = (*this)[i] / other[i];
    }
    return output;
}

Spectrum &Spectrum::operator+=(const Spectrum &other) {
    for (size_t i = 0; i < size(); i++) {
        (*this)[i] += other[i];
    }
    return *this;
}

Spectrum &Spectrum::operator-=(const Spectrum &other) {
    for (size_t i = 0; i < size(); i++) {
        (*this)[i] -= other[i];
    }
    return *this;
}

Spectrum &Spectrum::operator*=(const Spectrum &other) {
    for (size_t i = 0; i < size(); i++) {
        (*this)[i] *= other[i];
    }
    return *this;
}

Spectrum &Spectrum::operator/=(const Spectrum &other) {
    for (size_t i = 0; i < size(); i++) {
        if (std::abs(other[i]) == 0) (*this)[i] = {INFINITY, INFINITY};
        else (*this)[i] /= other[i];
    }
    return *this;
}

/* ------------------------------- */

Spectrum Spectrum::operator+(complexd value) const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = (*this)[i] + value;
    }
    return output;
}

Spectrum Spectrum::operator-(complexd value) const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = (*this)[i] - value;
    }
    return output;
}

Spectrum Spectrum::operator*(complexd value) const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = (*this)[i] * value;
    }
    return output;
}

Spectrum Spectrum::operator/(complexd value) const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        if (std::abs(value) == 0) output[i] = {INFINITY, INFINITY};
        else output[i] = (*this)[i] / value;
    }
    return output;
}

Spectrum &Spectrum::operator+=(complexd value) {
    for (size_t i = 0; i < size(); i++) {
        (*this)[i] += value;
    }
    return *this;
}

Spectrum &Spectrum::operator-=(complexd value) {
    for (size_t i = 0; i < size(); i++) {
        (*this)[i] -= value;
    }
    return *this;
}

Spectrum &Spectrum::operator*=(complexd value) {
    for (size_t i = 0; i < size(); i++) {
        (*this)[i] *= value;
    }
    return *this;
}

Spectrum &Spectrum::operator/=(complexd value) {
    for (size_t i = 0; i < size(); i++) {
        if (std::abs(value) == 0) (*this)[i] = {INFINITY, INFINITY};
        else (*this)[i] /= value;
    }
    return *this;
}

/* ------------------------------- */

Spectrum operator+(complexd value, const Spectrum &input) {
    Spectrum output(input.size());
    for (size_t i = 0; i < input.size(); i++) {
        output[i] = input[i] + value;
    }
    return output;
}

Spectrum operator-(complexd value, const Spectrum &input) {
    Spectrum output(input.size());
    for (size_t i = 0; i < input.size(); i++) {
        output[i] = input[i] - value;
    }
    return output;
}

Spectrum operator*(complexd value, const Spectrum &input) {
    Spectrum output(input.size());
    for (size_t i = 0; i < input.size(); i++) {
        output[i] = input[i] * value;
    }
    return output;
}

Spectrum operator/(complexd value, const Spectrum &input) {
    Spectrum output(input.size());
    for (size_t i = 0; i < input.size(); i++) {
        if (std::abs(value) == 0) output[i] = {INFINITY, INFINITY};
        else output[i] = input[i] / value;
    }
    return output;
}

/* ------------------------------- */

Spectrum Spectrum::cos() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::cos((*this)[i]);
    }
    return output;
}

Spectrum Spectrum::sin() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::sin((*this)[i]);
    }
    return output;
}

Spectrum Spectrum::tan() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::tan((*this)[i]);
    }
    return output;
}

Spectrum Spectrum::cosh() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::cosh((*this)[i]);
    }
    return output;
}

Spectrum Spectrum::sinh() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::sinh((*this)[i]);
    }
    return output;
}

Spectrum Spectrum::tanh() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::tanh((*this)[i]);
    }
    return output;
}

Spectrum Spectrum::acos() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::acos((*this)[i]);
    }
    return output;
}

Spectrum Spectrum::asin() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::asin((*this)[i]);
    }
    return output;
}

Spectrum Spectrum::atan() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::atan((*this)[i]);
    }
    return output;
}

Spectrum Spectrum::acosh() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::acosh((*this)[i]);
    }
    return output;
}

Spectrum Spectrum::asinh() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::asinh((*this)[i]);
    }
    return output;
}

Spectrum Spectrum::atanh() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::atanh((*this)[i]);
    }
    return output;
}

/* ------------------------------- */

Spectrum Spectrum::square() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = (*this)[i] * (*this)[i];
    }
    return output;
}

Spectrum Spectrum::pow(complexd exponent) const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::pow((*this)[i], exponent);
    }
    return output;
}

Spectrum Spectrum::sqrt() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::sqrt((*this)[i]);
    }
    return output;
}

Spectrum Spectrum::log() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::log((*this)[i]);
    }
    return output;
}

Spectrum Spectrum::log10() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::log10((*this)[i]);
    }
    return output;
}

Spectrum Spectrum::exp() const {
    Spectrum output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::exp((*this)[i]);
    }
    return output;
}

/* ------------------------------- */

Signal Spectrum::abs() const {
    Signal output(size());
    for (size_t i = 0; i < size(); i++) {
        output[i] = std::abs((*this)[i]);
    }
    return output;
}

complexd Spectrum::max() const
{
    if (this->empty()) {
        return NAN;
    }
    complexd maxElem = (*this)[0];
    for (size_t i = 0; i<this->size(); i++) {
        if (std::abs(maxElem) < std::abs((*this)[i])) maxElem = (*this)[i];
    }
    return maxElem;
}

complexd Spectrum::min() const {
    if (this->empty()) {
        return NAN;
    }
    complexd minElem = (*this)[0];
    for (size_t i = 0; i<this->size(); i++) {
        if (std::abs(minElem) > std::abs((*this)[i])) minElem = (*this)[i];
    }
    return minElem;
}

complexd Spectrum::mean() const {
    if (this->empty()) {
        return NAN;
    }
    complexd sum = 0;
    for (size_t i = 0; i<this->size(); i++) {
        sum += (*this)[i];
    }
    return sum / static_cast<complexd>(this->size());
}

/* ------------------------------- */

bool Spectrum::operator < (const Spectrum &input) const {
    for (size_t i = 0; i < size(); i++) {
        if (std::abs((*this)[i]) >= std::abs(input[i])) return false;
    }
    return true;
}

bool Spectrum::operator <= (const Spectrum &input) const {
    for (size_t i = 0; i < size(); i++) {
        if (std::abs((*this)[i]) > std::abs(input[i])) return false;
    }
    return true;
}

bool Spectrum::operator > (const Spectrum &input) const {
    for (size_t i = 0; i < size(); i++) {
        if (std::abs((*this)[i]) <= std::abs(input[i])) return false;
    }
    return true;
}

bool Spectrum::operator >= (const Spectrum &input) const {
    for (size_t i = 0; i < size(); i++) {
        if (std::abs((*this)[i]) < std::abs(input[i])) return false;
    }
    return true;
}

bool Spectrum::operator == (const Spectrum &input) const {
    for (size_t i = 0; i < size(); i++) {
        if ((*this)[i] != input[i]) return false;
    }
    return true;
}

bool Spectrum::operator != (const Spectrum &input) const {
    for (size_t i = 0; i < size(); i++) {
        if ((*this)[i] == input[i]) return false;
    }
    return true;
}

/* ------------------------------- */

Signal Spectrum::IDFT(size_t size_zero_padding) const {
    size_t N = this->size();
    size_t P = N + size_zero_padding;
    Spectrum reconstructedSignal(P);

    // Bit-reversal permutation
    size_t n = P;
    size_t bits = 0;
    while (n >>= 1) bits++;

    std::vector<size_t> reversed(P);
    for (size_t i = 0; i < P; i++) {
        size_t j = 0;
        for (size_t k = 0; k < bits; ++k) {
            if (i & (1 << k)) j |= (1 << (bits - 1 - k));
        }
        reversed[i] = j;
    }

    // Copy input data to reconstructedSignal with bit-reversed order
    size_t i;
    for (i = 0; i < N; i++) {
        reconstructedSignal[reversed[i]] = (*this)[i];
    }
    for (i = N; i < P; i++) {
        reconstructedSignal[reversed[i]] = 0.0;
    }

    // IFFT algorithm
    for (size_t s = 1; s <= bits; s++) {
        size_t m = 1 << s;
        complexd wm = std::exp(complexd(0.0, -2.0 * M_PI / m)); // Note the sign change
        for (size_t k = 0; k < P; k += m) {
            complexd w = 1;
            for (size_t j = 0; j < m / 2; j++) {
                complexd t = w * reconstructedSignal[k + j + m / 2];
                complexd u = reconstructedSignal[k + j];
                reconstructedSignal[k + j] = u + t;
                reconstructedSignal[k + j + m / 2] = u - t;
                w *= wm;
            }
        }
    }

    // Normalize by dividing by P
    for (size_t i = 0; i < P; i++) {
        reconstructedSignal[i] /= P;
    }

    // Ensure the signal is real by taking the real part of the complex values
    Signal realSignal(P); // Assuming sampling frequency is handled elsewhere
    for (size_t i = 0; i < P; i++) {
        realSignal[i] = reconstructedSignal[i].real();
    }

    return realSignal;
}

/* ------------------------------- */

std::ostream &operator<<(std::ostream &out, const Spectrum &spectrum) {
    out << spectrum.getName() << "{";
    for (size_t i = 0; i < spectrum.size(); i++) {
        out << spectrum[i].real() << (spectrum[i].imag() >= 0 ? "+" : "") << spectrum[i].imag() << "j";
        if (i < spectrum.size()-1) out << ",";
    }
    out << "}";
    return out;
}