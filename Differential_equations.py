# https://youtu.be/p_di4Zn4wz4?t=1464
import numpy as np

g = 9.8  # acceleration
L = 2  # the length of the line
mu = 0.1  # a generic term encapsulating air resistance


# the problem is that the angle depend on velocity, and velocity depend on acceleration and acceleration depend on angle
# we need to know the angle at specific point of time
#

# definition of ODE
def get_theta_double_dot(theta, theta_dot) :
    # calculate the acceleration based on the angle analysis of gravitational pull
    # and the velocity interaction with air resistance
    return -mu * theta_dot - (g / L) * np.sin(theta)


# Solution to the deferential equation
def theta(t, theta_0, theta_dot_0) :
    # the angle depend on  time initial angle and initial velocity (velocity is the derivative)
    # Initializing changing variables
    theta = theta_0  # initial angle
    # these variables are passed to the first iteration, then the loop will update them
    theta_dot = theta_dot_0  # initial velocity (used (along with the initial angle )to estimate the acceleration
    # acceleration in turn is used to  update the velocity, which updates the angle little by little
    # till we reach the time we are tyring to estimate the velocity at.
    # if big changes is used it will lose accuracy, it will use its predictive power to simulate reality
    # so apparently the changes that actually happens in reality is modeled with smooth changes
    delta_t = 0.01  # some time step
    for time in np.arange(0, t, delta_t) :
        # during these small step
        # the acceleration don't change depending on time but on the velocity and the angle
        theta_double_dot = get_theta_double_dot(theta, theta_dot)
        # the acceleration don't depend on time but on the velocity and the angle
        theta += theta_dot * delta_t  # each step changes the angle a little bit based on the velocity
        theta_dot += theta_double_dot * delta_t  # and change the velocity based on acceleration
    return theta  # the angle = the derivative multiplied by small step in time


for i in np.arange(0, 2000, 10) :
    print("When time is", i, ' seconds, the angle is',
          180 * theta(i, np.pi / 2, 0) / np.pi)  # initial angle of 90 and initial velocity of 0
