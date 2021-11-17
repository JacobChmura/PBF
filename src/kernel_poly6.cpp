#include <kernel_poly6.h>

double kernel_poly6(Particle &p_i, Particle &p_j, double h){
	double result;
	double r = (p_i.x_new - p_j.x_new).norm(); 
	if (r > h){
		// not within kernel radius
		result = 0.0;
	}
	else{
		result = (315.0/(64.0*PI*pow(h, 9))) * pow((pow(h,2)-pow(r,2)), 3);
	}
	return result;
}


double kernel_poly6(double r, double h){
	double result;
	
	if (r > h){
		// not within kernel radius
		result = 0.0;
	}
	else{
		result = (315.0/(64.0*PI*pow(h, 9))) * pow((pow(h,2)-pow(r,2)), 3);
	}
	return result;
}