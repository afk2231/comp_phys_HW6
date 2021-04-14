""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    Simplified and adapted by Lev Kaplan 2019"""

# rk4.py 4th order Runge Kutta


from pylab import *		 

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

def H(x,v):
    x = array(x)
    v = array(v)
    return (1/2)*(x*x + v*v)


def rk4(t,h,n):          # take one 4th order RK step          
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


a = 0.     # evolve from time a to time b in n steps
b = 100.
n = 10000
times, h = linspace(a, b, n, retstep=True)

y[0] = 1;   y[1] = 0    #initialize position and velocity
x_t = [y[0], ]
v_t = [y[1], ]
t_0 = []


for t in times[1:]:                                              # Time loop
    y = rk4(t,h,2)
    x_t.append(y[0])
    v_t.append(y[1])
    if v_t[-1] < 0 and v_t[-2] > 0:
        t_0.append((v_t[-1] * (t - h) - v_t[-2] * t)/(v_t[-1] - v_t[-2]))

y[0] = 2;   y[1] = 0    #initialize position and velocity
x_tp = [y[0], ]
v_tp = [y[1], ]
t_0p = []


for t in times[1:]:                                              # Time loop
    y = rk4(t, h, 2)
    x_tp.append(y[0])
    v_tp.append(y[1])
    if v_tp[-1] < 0 and v_tp[-2] > 0:
        t_0p.append((v_tp[-1] * (t - h) - v_tp[-2] * t)/(v_tp[-1] - v_tp[-2]))




print(y)   #print out final position and velocity
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

figure("interpolated t")
ylabel("$t(v=0)$")
xlabel("Number of intersections")
plot(t_0, '.', label="Interpolated times of zero velocity for $x_0 = 1$")
plot(t_0p, '.', label="Interpolated times of zero velocity for $x_0 = 2$")
x = linspace(0, len(t_0), len(t_0) + 1)
plot(x, (2 * pi) * (x + 1), label="$y(x) = 2\\pi (x+1)$")
legend()
show()

period = sum(gradient(t_0))/len(gradient(t_0))
print("Average calculated period to be {}".format(period))
print("Difference from expected period is {}".format(period - 2 * pi))

figure("difference in interpolated t")
ylabel("$\\Delta t(v=0)$")
xlabel("Number of intersections")
ylim(-1e-7,1e-7)
plot(array(t_0) - array(t_0p), '.', label="Difference in interpolated times of zero velocity")
show()

                                