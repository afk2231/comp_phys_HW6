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

# find the position and velocities as functions of time and the times when v=0 for initial position equal to 1
x1_t, v1_t, t1_0 = evolve([1, 0], times, h)

# find the position and velocities as functions of time and the times when v=0 for initial position equal to 2
x2_t, v2_t, t2_0 = evolve([2, 0], times, h)

#plot the position and energy of the harmonic oscillator
figure("position and energy")
subplot(211)
ylabel("$x(t)$")
plot(times, x1_t)
plot(times, cos(times), "--")
subplot(212)
ylabel("$E(t)$")
xlabel("Times (a.u.)")
plot(times, H(x1_t, v1_t))
plot(times, 1/2*ones_like(times))
show()

# plot the interpolated times when v=0 at different initial positions
figure("interpolated t")
ylabel("$t(v=0)$")
xlabel("Number of intersections")
plot(t1_0, '.', label="Interpolated times of zero velocity for $x_0 = 1$")
plot(t2_0, '.', label="Interpolated times of zero velocity for $x_0 = 2$")
x = linspace(0, len(t1_0), len(t1_0) + 1)
plot(x, (2 * pi) * (x + 1), label="$y(x) = 2\\pi (x+1)$")
legend()
show()

# calculate the average period of these times
period = sum(gradient(t1_0))/len(gradient(t1_0))
print("Average calculated period to be {}".format(period))
print("Difference from expected period is {}".format(period - 2 * pi))

# plot the difference between the times when v=0
figure("difference in interpolated t")
ylabel("$\\Delta t(v=0)$")
xlabel("Number of intersections")
ylim(-1e-7,1e-7)
plot(array(t1_0) - array(t2_0), '.', label="Difference in interpolated times of zero velocity")
show()
