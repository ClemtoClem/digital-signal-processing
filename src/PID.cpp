#include "PID.hpp"

PID::PID() : is_set(false)
{
}

PID::~PID()
{
}

bool PID::set(double kp, double ki, double kd, double dt, double max, double min)
{
    if (dt <= 0.0) {
        std::cerr << "Impossible to create a PID regulator with a loop interval time <= 0" << std::endl;
        return false;
    }

    // Set values
    this->kp = kp;
    this->ki = ki;
    this->kd = kd;
    this->dt = dt;
    this->max = max;
    this->min = min;

    // Reset the PID
    reset();

    // Set the PID as set
    is_set = true;
    return true;
}

double PID::process(double desired_value, double current_value)
{
    if (!is_set) {
        std::cerr << "PID regulator is not set" << std::endl;
        return 0.0;
    }
    
    // Calculate error
    double error = current_value - desired_value;
    
    // Proportional term
    double p_term = kp * error;

    // Integral term
    integral += error * dt;
    double i_term = ki * integral;

    // Derivative term
    double derivative = (error - last_error) / dt;
    double d_term = kd * derivative;

    // Calculate the control output
    double output = p_term + i_term + d_term;

    // Limit the output to the max and min values
    if (output > max) output = max;
    else if (output < min) output = min;

    // save the last error
    last_error = error;

    return output;
}

void PID::reset() {
    integral = 0.0;
    last_error = 0.0;
}
