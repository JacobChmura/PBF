#include <Eigen/Dense>
#include <vector>

/*
Apply Viscocity Smoothing to the simulated fluid. 

Input:
	const Eigen::Ref<const Eigen::MatrixXd> &x_new: Position of particles in the system on the after Jacobi Update.
	std::vector<std::vector<int>> &neighbours: Vector at index (i) stores indices of neighbouring particles to particle i.
	Eigen::MatrixXd &v_new: temporary viscocity accounted for velocity
	double viscocity_c: the scaling factor for viscocity velocity update
	double kernel_h: the kernel radius for the fluid simulation.

Modifies:
	Eigen::MatrixXd &v: updated velocity of the fluid system after XPSH viscocity.
*/
void apply_viscocity(const Eigen::Ref<const Eigen::MatrixXd> &x_new, std::vector<std::vector<int>> & neighbours, 
					Eigen::MatrixXd &v, Eigen::MatrixXd &v_new, double viscocity_c, double kernel_h);
