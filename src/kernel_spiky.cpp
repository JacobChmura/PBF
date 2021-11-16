#include <kernel_spiky.h>

void kernel_spiky(Eigen::Vector3d &result, Particle &p_i, Particle &p_j, double h){
	result.setZero();

	double r = (p_i.x_new - p_j.x_new).norm(); 

	if (r <= h and r >= 0.00000001){
		// Distance within kernel radius and non-negligble
		result = -(45.0/(PI*pow(h, 6))) * (pow(h-r, 2)/r) * (p_i.x_new - p_j.x_new);
	}
}