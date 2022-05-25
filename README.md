# symWannier

A collection of codes to construct Wannier functions using crystal symmetry.

## How to use

We need to modify pw2wannier90 and wannier90.

Copy src/pw2wannier90/qe6.6/pw2wannier90.f90 to your QE directory and compile it.

Copy src/wannier90-3.1.0/kmesh.f90 to your wannier90 directory and compile it. (This is just for writing bvectors in nnkp file.)

There are sample input files and a script file, test.sh, in examples/Cu.
Please modify test.sh and run.

