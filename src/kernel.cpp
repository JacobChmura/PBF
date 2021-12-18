#include <kernel.h>

void kernel_spiky(Eigen::Vector3d &result, Eigen::Ref<const Eigen::RowVector3d> p_i, Eigen::Ref<const Eigen::RowVector3d> p_j, double h){
	result.setZero();

	double r = (p_i - p_j).norm(); 
	if (r <= h and r >= NUMERICAL_EPS){ // Distance within kernel radius and non-negligble
		result = -(45.0/(PI*pow(h, 6))) * (pow(h-r, 2)/r) * (p_i - p_j);
	}
}

double kernel_poly6(Eigen::Vector3d p_i, Eigen::Vector3d p_j, double h){
	double result = 0;
	
	double r = (p_i - p_j).norm(); 
	if (r <= h) result = (315.0/(64.0*PI*pow(h, 9))) * pow((pow(h,2)-pow(r,2)), 3);
	return result;
}

double kernel_poly6(double r, double h){
	double result = 0;
	if (r <= h) result = (315.0/(64.0*PI*pow(h, 9))) * pow((pow(h,2)-pow(r,2)), 3);
	return result;
}
