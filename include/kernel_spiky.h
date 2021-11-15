#include <Particle.h>

/*
Compute the spiky kernel between two particles for gradient estimation. 

Input:
	Particle p_i: A particle in the fluid simulation.
	Particle p_j: Another particle in the fluid simulation.
	double h: Kernel Radius.

Output:
	double result: The output of the spiky kernel based on the position of p_i, p_j.
*/
double kernel_spiky(Particle &p_i, Particle &p_j, double h);