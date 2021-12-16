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
void apply_vorticity(const Eigen::Ref<const Eigen::MatrixXd> &x_new, std::vector<std::vector<int>> &neighbours, Eigen::MatrixXd &v, Eigen::MatrixXd omega, Eigen::MatrixXd eta, Eigen::MatrixXd N, Eigen::MatrixXd vorticity_f, double vorticity_epsilon, double kernel_h, double dt);
