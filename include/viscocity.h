#include <Eigen/Dense>
#include <vector>

/*
Apply Viscocity Smoothing to the simulated fluid. 

Input:
	std::vector<Particle> fluid: List of particles in the fluid simulation.
	double kernel_h: the kernel radius for the fluid simulation.
	double viscocity_c: the scaling factor for viscocity velocity update

Output:
	modifies the velocity of each particle in the fluid by performing viscocity smoothing.
*/
//void apply_viscocity(std::vector<Particle> &fluid, double kernel_h, double viscocity_c);
