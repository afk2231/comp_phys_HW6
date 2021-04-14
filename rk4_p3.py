""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    Simplified and adapted by Lev Kaplan 2019"""

# rk4.py 4th order Runge Kutta


from pylab import *
from scipy.special import ellipk
from functions import generate_period

pltparams = {
    'text.usetex': True
}
rcParams.update(pltparams)


a = 0.     # evolve from time a to time b in n steps
b = 100.
n = int(1e4)
times, h = linspace(a, b, n, retstep=True)
thetas = linspace(0.01, 0.99 * pi, 300)
big_ts = []

for theta in thetas:
    big_ts.append(generate_period(theta, times, h))

figure("Period")
ylabel("$T(\\theta_0)$")
xlabel("$\\theta_0$")
plot(thetas, big_ts, label="Calculated periods")
plot(thetas, 4 * ellipk((sin(thetas / 2))**2), "--", label="Analytical value of the periods")
legend()
show()
