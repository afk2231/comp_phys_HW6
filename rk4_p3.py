""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    Simplified and adapted by Lev Kaplan 2019"""

# rk4.py 4th order Runge Kutta


from pylab import *
from scipy.special import ellipk

pltparams = {
    'text.usetex': True
}
rcParams.update(pltparams)

ydumb = zeros(2);       y = zeros(2)
fReturn = zeros(2);     k1 = zeros(2, float)
k2 = zeros(2, float);   k3 = zeros((2), float) 
k4 = zeros((2), float)

def f( t, y):      # Force function: component 0 is position, component 1 is velocity
    fReturn[0] = y[1]       # dx/dt  = v                                     
    fReturn[1] = -y[0]      # dv/dt = -x   
    return fReturn

def f_pendulum(t, y):      # Force function: component 0 is position, component 1 is velocity
    fReturn[0] = y[1]       # d theta/dt  = omega
    fReturn[1] = -sin(y[0])      # d omega/dt = -sin(theta)
    return fReturn

def H(x,v):
    x = array(x)
    v = array(v)
    return (1/2)*(x*x + v*v)


def rk4(t,y,h,n,f):          # take one 4th order RK step
                          # evolve y from time t to time t+h
                          # n is number of variables to evolve
                          # k1,k2,k3,k4 are four estimates for Delta y
    k1 = h*f(t, y)       # note that we use vector (array) notation instead of using loops
                         # see text for how to do this the long way with loops
    ydumb = y + k1/2     #estimate for midpoint using Euler
    k2 = h*f(t+h/2, ydumb)
    ydumb = y + k2/2     #another estimate for midpoint
    k3 = h*f(t+h/2, ydumb)
    ydumb = y + k3       #estmate for endingn poitn of interval
    k4 = h*f(t+h, ydumb)
    ynew = y + (k1 + 2*(k2 + k3) + k4)/6
    return ynew

def evolve(y_0, times):
    x = [y_0[0], ]
    v = [y_0[1], ]
    t_0 = []
    y = y_0
    for t in times[1:]:  # Time loop
        y = rk4(t, y, h, 2)
        x.append(y[0])
        v.append(y[1])
        if v[-1] < 0 and v[-2] > 0:
            t_0.append((v[-1] * (t - h) - v[-2] * t) / (v[-1] - v[-2]))

    return [x, v, t_0]

def generate_period(x_0, times):
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

    return sum(gradient(t_0))/len(gradient(t_0))

a = 0.     # evolve from time a to time b in n steps
b = 100.
n = int(1e4)
times, h = linspace(a, b, n, retstep=True)
thetas = linspace(0.01, 0.99 * pi, 300)
big_ts = []

for theta in thetas:
    big_ts.append(generate_period(theta, times))

figure("Period")
ylabel("$T(\\theta_0)$")
xlabel("$\\theta_0$")
plot(thetas, big_ts, label="Calculated periods")
plot(thetas, 4 * ellipk((sin(thetas / 2))**2), "--", label="Analytical value of the periods")
legend()
show()
