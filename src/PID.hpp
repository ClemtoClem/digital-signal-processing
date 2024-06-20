#ifndef __PID_HPP
#define __PID_HPP

#include <iostream>

// Implement a PID controller
class PID {
private:
    double kp;
    double ki;
    double kd;
    double dt;
    double integral;
    double max;
    double min;
    double last_error;

    bool is_set;
public:
    
    PID();

    ~PID();

    /**
     * Set the PID controller parameters
     * @param[in] kp Proportional gain
     * @param[in] ki Integral gain
     * @param[in] kd Derivative gain
     * @param[in] dt Loop interval time
     * @param[in] max Maximum value of the manipulated variable
     * @param[in] min Minimum value of the manipulated variable
     */
    bool set(double kp, double ki, double kd, double dt, double max, double min);

    // Returns the manipulated variable given current and current prosess value
    double process(double desired_value, double current_value);

    // Reset the PID controller
    void reset();
};

#endif // __PID_HPP