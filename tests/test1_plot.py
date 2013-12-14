# This comes from the Dullemond lectures
# http://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2012/Chapter_4.pdf

from pylab import *

a, b = loadtxt('outputs/output.txt',unpack=True)
plot(a,b)
show()
