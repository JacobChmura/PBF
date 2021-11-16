#include <Eigen/Dense>
#include <Particle.h>
#include <list>
#include <kernel_spiky.h>
/*
Apply Vorticity Confinement to the simulated fluid. 

Input:
	std::list<Particle> fluid: List of particles in the fluid simulation.
	double kernel_h: the kernel radius for the fluid simulation.

Output:
	modifies the velocity of each particle in the fluid by adding forces due to vorticity.
*/
void apply_vorticity(std::list<Particle> &fluid, double kernel_h);