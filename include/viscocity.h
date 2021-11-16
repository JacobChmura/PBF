#include <Eigen/Dense>
#include <Particle.h>
#include <list>
#include <kernel_poly6.h>

/*
Apply Viscocity Smoothing to the simulated fluid. 

Input:
	std::list<Particle> fluid: List of particles in the fluid simulation.
	double kernel_h: the kernel radius for the fluid simulation.
	double viscocity_c: the scaling factor for viscocity velocity update

Output:
	modifies the velocity of each particle in the fluid by performing viscocity smoothing.
*/
void apply_viscocity(std::list<Particle> &fluid, double kernel_h, double viscocity_c);