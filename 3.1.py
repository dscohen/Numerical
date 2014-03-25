from scipy import *
from math import *

def rootfind(f, fdiv, a, b, nprobe, tol):
    #set increment for b and a
    interv = (b-a)/float(nprobe)
    current_probe = f(a)
    #set subinterval a1, b1
    a1 = a
    b1 = interv+a
    roots = []
    #whie b1 doesn't pass border
    while (b1 <= b):
        next_probe = f(b1)
        #means root is within this area
        if ((current_probe * next_probe) <= 0):
            #perform newton's method
            fa1 = current_probe
            fb1 = next_probe
            x_old = a1
            x_new = x_old - f(x_old)/fdiv(x_old)
            while((abs(x_new - x_old) > tol*(1+(abs(x_new)))) or (f(x_new) > tol)):
                x_new = x_old - f(x_old)/fdiv(x_old)
                if (abs(f(x_new)) >= 0.5*abs(f(x_old))):
                    a1, b1 = bisect3(f, a1, b1, fa1, fb1)
                    fa1 = f(a1)
                    fb1 = f(b1)
                    x_old = a1
                    x_new = x_old - f(x_old)/fdiv(x_old)
            roots.append(x_new)
        #move the interval
        a1 = b1
        b1 += interv
        current_probe = next_probe
    print roots



def bisect3(f, a, b, fa, fb):
    p = (a+b)/2
    for i in xrange(3):
        fp = f(p)
        if (fa*fp < 0):
            b = p
        else:
            a = p
            fa = fp
        p = (a+b)/2
    return a, b
    
def func(x):
    x1 = x/4.0
    x = 2*((e**x1 + e**(-x1))/2.)-x
    return x
    
def funcderiv(x):
    x1 = x/4.0
    x = -e**-x1 + e**x1 - 1
    return x

