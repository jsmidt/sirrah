# This comes from the Dullemond lectures
# http://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2012/Chapter_4.pdf

from pylab import *

a, b = loadtxt('outputs/test2-000000.txt',unpack=True)
plot(a,b)
a, b = loadtxt('outputs/test2-000700.txt',unpack=True)
plot(a,b)
title('Compare with figure 4.7')
show()
