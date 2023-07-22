How to run and compile my program: 

type 'make' to compile code 

write 'sbatch jobscript.slm' to send machine code to cluster 

PROGRAMS MODIFIED: 

 - Newton: Inexact newtons method for non-lienar diffusion equation
 - func: computes F(U)
 - Jacobian: computes F'(U)
 - vecdot: computes vector dot product for distributed vectors a,b 
 - matmult: computes matrix vector product 
 - negu: determines if computed solution U is negative 
 - cg: computes CG algorithm to find solution to Ax = b
 - main: main function for running newton, returns results and saves solution
 - input: file which contains 'm' and 'lambda' (can be changed by user)
   
