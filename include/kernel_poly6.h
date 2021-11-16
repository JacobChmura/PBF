#include <Particle.h>
#include <cmath>

#define PI 3.14159265358979323846 

/*
Compute the Poly6 kernel between two particles for density estimation. 

Input:
	Particle p_i: A particle in the fluid simulation.
	Particle p_j: Another particle in the fluid simulation.
	double h: Kernel Radius.

Output:
	double result: The output of the Poly6 kernel based on the position of p_i, p_j.
*/
double kernel_poly6(Particle &p_i, Particle &p_j, double h);
