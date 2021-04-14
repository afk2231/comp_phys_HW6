""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    Simplified and adapted by Lev Kaplan 2019"""

# rk4.py 4th order Runge Kutta


from pylab import *
from functions import evolve, H

pltparams = {
    'text.usetex': True
}
rcParams.update(pltparams)

a = 0.     # evolve from time a to time b in n steps
b = 100.
n = 10000
times, h = linspace(a, b, n, retstep=True)

y_0 = [1,0]   #initialize position and velocity
x_t = [y_0[0], ]
v_t = [y_0[1], ]

# find the position and velocities as functions of time and the times when v=0
x1_t, v1_t, t1_0 = evolve(y_0, times, h)


#plot position and energy of the harmonic oscillator
figure("position and energy")
subplot(211)
ylabel("$x(t)$")
plot(times, x_t)
plot(times, cos(times), "--")
subplot(212)
ylabel("$E(t)$")
xlabel("Times (a.u.)")
plot(times, H(x_t, v_t))
plot(times, 1/2*ones_like(times))
show()
