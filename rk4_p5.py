""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    Simplified and adapted by Lev Kaplan 2019"""

# rk4.py 4th order Runge Kutta


from pylab import *
from functions import *

def f_pendulumD(t, y):
    # Force function for the driven pendulum
    global wDr
    FDr = 0.1 # Drive amplitude
    fReturn[0] = y[1]  # d theta/dt  = omega
    fReturn[1] = -sin(y[0]) + FDr*sin(wDr*t)# d omega/dt = -sin(theta)
    return fReturn

pltparams = {
    'text.usetex': True
}
rcParams.update(pltparams)


a = 0.     # evolve from time a to time b in n steps
b = 100.
n = 10000
times, h = linspace(a, b, n, retstep=True)

y_0 = [0,0]   #initialize position and velocity
x_t = [y_0[0], ]
v_t = [y_0[1], ]

drives = np.linspace(0.75,0.75,1)

amps = []

for wDr in drives:
    # find the position and velocities as functions of time and the times when v=0
    x1_t, v1_t, t1_0 = evolvePend(y_0, times, h, f_pendulumD)
    amps.append(max(x1_t)-min(x1_t))

#"""
#plot position and energy of the harmonic oscillator
figure("Driven Pendulum Response")
ylabel("$x(t)$")
xlabel("Time")
title("Driven Pendulum at 0.75 of Resosnant Frequency")
plot(times, 0.2*cos(times), "--", label = "Resonant Frequency")
plot(times, 0.2*cos(wDr*times), "--", label = "Drive Frequency")
plot(times, x1_t, label = "Pendulum Response")
legend()
show()

"""
figure("Driven Pendulum")
ylabel("Amplitude")
xlabel("Frequency")
plot(drives, amps)
show()
"""
