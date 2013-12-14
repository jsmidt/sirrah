SIRRAH
====== 

A *S*pec*I*al-*R*elativistic *RA*diation *H*ydrodynamics Code

About
-----

As the name suggests, SIRRAH is a publicly available special-relativistic
radiation hydrodynamics. It was designed in my free time (see disclaimer below)
mainly for studying supernova blasts, but the attempt is made to make the code
general. All contributions are welcome.

Running
-------

To run the code, edit the Makefile in the src/ directory to fit your complier
and point to fortranlib by astrofrog [2] which is a build dependancy.  Then run
the traditional 'make' command to compile. (This can also be run in the base
directory) After the code is compiled, an executable 'sirrah' is created in the
base directory. To run type:

    ./sirrah my_infile.txt

where `my_infile.txt` is a file containing your initial parameters. The output
of the run is stored in an outputs directory created by the code.

Please see the tests/ directory for examples. This directory contains both
example parameter files, profile files and python plotting scripts to give you
a feel how to run your own problem. For example, the first test problem can be run as:

    ./sirrah tests/test1_infile.txt

To see a plot of the solution, run the accompanying python script:

    python tests/test1_plot.py 

Most of these tests were taken from the Dullemond lectures [1] and so plots can
be compared against his results. 

License
-------

This code is licensed with the BSD license to ensure maximal availability for
the public while retaining recognition of the original author (Joseph Smidt
<josephsmidt@gmail.com>) and copyright as of 2013. Please see the licence for
more details.

Disclaimer
----------

No Department of Energy resources (including email, machines or DOE Codes for
reference) were used for the creation of this code. Just me, my personal
computer, the Dullemond lectures [1] coding during my free time.

[1] http://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2012/index.shtml
[2] https://github.com/astrofrog/fortranlib
