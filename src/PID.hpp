#ifndef __PID_HPP
#define __PID_HPP

#include <iostream>

// Implement a PID controller
class PID {
private:
    double kp;
    double ki;
    double kd;
    double outMin;
    double outMax;
    double integral;
    double last_error;
    double deltaTime;

    // Additional parameters for self-tuning
    double alpha; // Adaptation rate for kp
    double beta;  // Adaptation rate for ki
    double gamma; // Adaptation rate for kd

    bool is_set;

public:
    
    PID();

    ~PID();

    /**
     * Set the PID controller parameters
     * @param[in] kp Proportional gain
     * @param[in] ki Integral gain
     * @param[in] kd Derivative gain
     * @param[in] outMin Minimum value of the manipulated variable
     * @param[in] outMax Maximum value of the manipulated variable
     * @param[in] alpha Adaptation rate for kp
     * @param[in] beta Adaptation rate for ki
     * @param[in] gamma Adaptation rate for kd
     */
    bool set(double kp, double ki, double kd, double outMin = -1, double outMax = 1, double alpha = 0.0, double beta = 0.0, double gamma = 0.0);

    /**
     * Set the PID controller parameters
     * @param[in] kp Proportional gain
     * @param[in] ki Integral gain
     * @param[in] kd Derivative gain
     */
    bool setGains(double kp, double ki, double kd);

    /**
     * Set the PID controller output limits
     * @param[in] outMin Minimum value of the manipulated variable
     * @param[in] outMax Maximum value of the manipulated variable
     */
    bool setOutputLimits(double outMin, double outMax);

    double getKp() { return kp; }
    double getKi() { return ki; }
    double getKd() { return kd; }
    double getOutMax() { return outMax; }
    double getOutMin() { return outMin; }

    // Returns the manipulated variable given current and current process value
    double process(double desired_value, double input);

    // Reset the PID controller
    void reset();
};


#endif // __PID_HPP