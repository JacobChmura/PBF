#include <Particle.h>
#include <Eigen/Dense>
#include <cmath>

#define PI 3.14159265358979323846 

/*
Compute the spiky kernel between two particles for gradient estimation. 

Input:
	Particle p_i: A particle in the fluid simulation.
	Particle p_j: Another particle in the fluid simulation.
	double h: Kernel Radius.

Output:
	Eigen::Vector3d result: The output of the spiky kernel based on the position of p_i, p_j.
*/
void kernel_spiky(Eigen::Vector3d &result, Particle &p_i, Particle &p_j, double h);


/*
Compute the Poly6 kernel between two particles for density estimation. 

Input:
	Particle p_i: A particle in the fluid simulation.
	Particle p_j: Another particle in the fluid simulation.
	double h: Kernel Radius.

Output:
	double result: The output of the Poly6 kernel based on the position of p_i, p_j.
*/
double kernel_poly6(Particle &p_i, Particle &p_j, double h);


/*
Compute the Poly6 kernel given an explicit radius r from origin. 

Input:
	double r: radius of interaction
	double h: Kernel Radius.

Output:
	double result: The output of the Poly6 kernel based on the radius r.
*/
double kernel_poly6(double r, double h);
