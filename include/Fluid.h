#ifndef FLUID
#define FLUID

#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <Particle.h>
#include <Eigen/Dense>

class Fluid{
public:
	// Basic State Parameters
	std::vector<Particle> fluid; // Particles comprising the fluid system
	Eigen::MatrixXd fluid_state;
	int num_particles;

	double particle_mass; // Mass of each particle
	double rho; // Rest density
	double gravity_f; // Force of gravity on the system
	double user_f; // Force of user applied force on the system

	// Jacobi Parameters
	int jacobi_iterations;

	// Constraint Force Mixing Parameter
	double cfm_epsilon;

	// Kernel Parameters
	double kernel_h;

	// Tensile Stability Parameters
	double tensile_k;
	double tensile_delta_q;
	int tensile_n;

	// Viscocity Parameters
	double viscocity_c;

	// Vorticity Parameters
	double vorticity_epsilon;

	// Simulation Parameters
	double t;
	double dt;

	// Bounding box
	double lower_bound;
	double upper_bound;

	// Temp variables
	Eigen::Vector3d ker_res;
	double s_corr;

public:
	
	//Initialize a blank new fluid.
	Fluid(double particle_mass, double rho, double gravity_f, double user_f, int jacobi_iterations, 
			double cfm_epsilon, double kernel_h, double tensile_k, double tensile_delta_q, int tensile_n, 
			double viscocity_c, double vorticity_epsilon, double lower_bound, double upper_bound, double dt);

	// Load a scene (random or dam break or double dam break)
	void init_state(Eigen::MatrixXd &fluid_state);
	
	//Perform a timestep update.
	void step(Eigen::MatrixXd &fluid_state);

};

#endif