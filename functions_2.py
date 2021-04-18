from pylab import *

ydumb = zeros(2);
fReturn = zeros(2);
k1 = zeros(2, float)
k2 = zeros(2, float);
k3 = zeros((2), float)
k4 = zeros((2), float)


def f(t, y):
    # Force function: component 0 is position, component 1 is velocity
    fReturn[0] = y[1]  # dx/dt  = v
    fReturn[1] = -y[0]  # dv/dt = -x
    return fReturn


def f_pendulum(t, y):
    # Force function for the pendulum
    fReturn[0] = y[1]  # d theta/dt  = omega
    fReturn[1] = -sin(y[0])  # d omega/dt = -sin(theta)
    return fReturn

def f_pendulumDD(t, y):
    # Force function for the danped driven pendulum
    FDr = 0.5 # Drive amplitude
    wDr = 1.5 # Drive Frequency
    G = 0.01 # Damping
    fReturn[0] = y[1]  # d theta/dt  = omega
    fReturn[1] = -sin(y[0]) + FDr*sin(wDr*t) - G*y[1] # d omega/dt = -sin(theta)
    return fReturn


def H(x, v):  # energy of the harmonic oscillator
    x = array(x)
    v = array(v)
    return (1 / 2) * (x * x + v * v)


def rk4(t, y, h, n, f):
    # take one 4th order RK step
    # evolve y from time t to time t+h
    # n is number of variables to evolve
    # k1,k2,k3,k4 are four estimates for Delta y
    k1 = h * f(t, y)  # note that we use vector (array) notation instead of using loops
    # see text for how to do this the long way with loops
    ydumb = y + k1 / 2  # estimate for midpoint using Euler
    k2 = h * f(t + h / 2, ydumb)
    ydumb = y + k2 / 2  # another estimate for midpoint
    k3 = h * f(t + h / 2, ydumb)
    ydumb = y + k3  # estmate for ending point of interval
    k4 = h * f(t + h, ydumb)
    ynew = y + (k1 + 2 * (k2 + k3) + k4) / 6
    return ynew


def evolve(y_0, times, h):
    # evolves a the harmonic oscillator from an initial state to some final time
    x = [y_0[0], ]
    v = [y_0[1], ]
    t_0 = []
    y = y_0
    for t in times[1:]:  # Time loop
        y = rk4(t, y, h, 2, f)
        x.append(y[0])
        v.append(y[1])
        if v[-1] < 0 and v[-2] > 0:
            t_0.append((v[-1] * (t - h) - v[-2] * t) / (v[-1] - v[-2]))

    return [x, v, t_0]


def generate_period(x_0, times, h):
    # finds the period of a pendulum from some initial angle x_0
    x = [x_0, ]
    v = [0, ]
    t_0 = []
    y = [x_0, 0]
    for t in times[1:]:  # Time loop
        y = rk4(t, y, h, 2, f_pendulum)
        x.append(y[0])
        v.append(y[1])
        if v[-1] < 0 and v[-2] > 0:
            t_0.append((v[-1] * (t - h) - v[-2] * t) / (v[-1] - v[-2]))

    return sum(gradient(t_0)) / len(gradient(t_0))

def evolvePend(y_0, times, h, f_p):
    # evolves a pendulum function from an initial state to some final time
    x = [y_0[0], ]
    v = [y_0[1], ]
    t_0 = []
    y = y_0
    for t in times[1:]:  # Time loop
        y = rk4(t, y, h, 2, f_p)
        x.append(y[0])
        v.append(y[1])
        if v[-1] < 0 and v[-2] > 0:
            t_0.append((v[-1] * (t - h) - v[-2] * t) / (v[-1] - v[-2]))

    return [x, v, t_0]

