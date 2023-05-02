# MPI-FORTRAN-90-Nonlinear-Thermal-Conduction
Nonlinear Thermal Conduction: Parallel programming using MPI with FORTRAN 90. NOTE: Parallel storage scheme CRS, BLAS and LAPACK used for efficiency. Programmed using EMACS, and tests ran by submitting jobscript.slm on cluster (super computer based in europe)

This assignment is about practical aspects of solving sparse systems of nonlinear equations using the inexact Newton method. In particular we will use parallel Newton–CG and look at a case study on nonlinear thermal combustion in a self–heating medium in a partially insulated square domain. This is the same model as the one considered in the first assignment. This problem arises when investigating critical parameter values for underground repositories of self–heating waste which are partially covered by buildings. Beyond certain critical parameter values the steady state solutions may become very large or even unbounded leading to an explosion in the medium. 

For more information, see 

[1] Adler J., Thermal–explosion theory for a slab with partial insulation, Combustion and
Flame 50, 1983, pp. 1–7.

[2] Greenway P. and Spence A., Numerical calculation of critical points for a slab with
partial insulation, Combustion and Flame 62, 1985, pp. 141–156. 

[3] Kelley CT., Iterative Methods for Linear and Nonlinear Equations, SIAM, Philadelphia,
1995. 
