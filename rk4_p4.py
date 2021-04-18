""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    Simplified and adapted by Lev Kaplan 2019"""

# rk4.py 4th order Runge Kutta


from pylab import *
from scipy.special import ellipk
from functions import *

pltparams = {
    'text.usetex': True
}
rcParams.update(pltparams)


a = 0.     # evolve from time a to time b in n steps
b = 550.
n = 10000
times, h = linspace(a, b, n, retstep=True)

y_0 = [1,0]   #initialize position and velocity
x_t = [y_0[0], ]
v_t = [y_0[1], ]

# find the position and velocities as functions of time and the times when v=0
x1_t, v1_t, t1_0 = evolvePend(y_0, times, h, f_pendulumDD)


#plot position and energy of the harmonic oscillator
figure("Damped Driven Pendulum")
ylabel("$x(t)$")
xlabel("Time")
title("Damped Driven Response v. Time")
plot(times, 0.5*cos(1.5*times), "--", label = "Driving Force")
plot(times, x1_t, label = "Pendulum Response")
legend()
show()