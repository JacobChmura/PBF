#include <Eigen/Dense>
#include <vector>
/*
Apply Vorticity Confinement to the simulated fluid. 

Input:
	std::vector<Particle> fluid: List of particles in the fluid simulation.
	double kernel_h: the kernel radius for the fluid simulation.
	double vorticity_epsilon: the scaling factor for vorticity force.
	doubel dt: simulation time step
Output:
	modifies the velocity of each particle in the fluid by adding forces due to vorticity.
*/
//void apply_vorticity(std::vector<Particle> &fluid, double kernel_h, double vorticity_epsilon, double dt);
