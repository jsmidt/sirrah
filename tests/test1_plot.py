# This comes from the Dullemond lectures
# http://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2012/Chapter_4.pdf

from pylab import *

a, b = loadtxt('outputs/test1-000000.txt',unpack=True)
plot(a,b)
a, b = loadtxt('outputs/test1-000100.txt',unpack=True)
plot(a,b)
a, b = loadtxt('outputs/test1-000200.txt',unpack=True)
plot(a,b)
a, b = loadtxt('outputs/test1-000300.txt',unpack=True)
plot(a,b)
show()
