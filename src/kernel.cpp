#include <kernel.h>

void kernel_spiky(Eigen::Vector3d &result, Particle &p_i, Particle &p_j, double h){
	result.setZero();

	double r = (p_i.x_new - p_j.x_new).norm(); 
	if (r <= h and r >= 0.00000001){ // check p_i != p_j divide by 0
		// Distance within kernel radius and non-negligble
		result = -(45.0/(PI*pow(h, 6))) * (pow(h-r, 2)/r) * (p_i.x_new - p_j.x_new);
	}
}

double kernel_poly6(Particle &p_i, Particle &p_j, double h){
	double result = 0;
	double r = (p_i.x_new - p_j.x_new).norm(); 
	if (r <= h) result = (315.0/(64.0*PI*pow(h, 9))) * pow((pow(h,2)-pow(r,2)), 3);
	return result;
}


double kernel_poly6(double r, double h){
	double result = 0;
	if (r <= h) result = (315.0/(64.0*PI*pow(h, 9))) * pow((pow(h,2)-pow(r,2)), 3);
	return result;
}
