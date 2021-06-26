/*  ~ Bare-bones implementation of pid.h ~ */ 

#include<iostream>
#include "pid.h"

int main()
{
    /* PID params */
    float kp = 0.01; 
    float ki = 1.001;
    float kd = 0.0;
    float set_point = 10;
    float dt = 0.01; // Sample time

    // Instantiate PID class
    PID control(kp, ki, kd);
    control.setOutputLimits(100,-100);

    /* In each timestep */
    float process_value = 1; // Plant feedback
    float error = set_point - process_value;
    control.update(error, dt); // Run PID

    // Display output
    std::cout << "PID output is: " << control.output; 

    return 0;
}
