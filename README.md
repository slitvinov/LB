<h2>Lattice Boltzmann methods</h2>

To build
<par>
$ make
gfortran -O2 -g bgk2.f -o bgk2
</par>

To run
<par>
$ ./bgk2 < bgk2.inp
 Number of steps
 Number of steps between printing profile
 Number of steps between performing diagnostics
 Relaxation frequency omega
 Applied force(T or F)) ?
 Initial density and velocity for the Poiseuille force
 Final velocity for the Poise force
 Linear obstacle(T or F)?
 Obstacle height?
 Obstacle id           1
 Length of the obstacle(multple of 2)           8
 File for output : 5 chars
 *****************************************
 Lattice BGK model, 2D with 9 velocities
...
</par>

To postprocess
<par>
$ python post.py  *.raw
</par>

<p align="center"><img src="img/u.gif"/></p>
