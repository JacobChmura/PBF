#ifndef FLUID
#define FLUID

#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <SpatialHashGrid.h>
#include <Eigen/Dense>

class Fluid{
public:
	// Global State Parameters
        // Matrices are of size (num_particles x 3) store vector quantities for each particles
        // Vectors of size (num_particles) store scalar quantities for each particles
        // -----
        
        Eigen::MatrixXd x_new, v, dP, omega, eta, N, vorticity_f, cell_coord;
        Eigen::VectorXd density, c, lambda, c_grad_norm, gravity_f;
        std::vector<std::vector<int>> neighbours;

        // ----
	//std::vector<Particle> fluid; // Particles comprising the fluid system
	//Eigen::MatrixXd fluid_state;
        SpatialHashGrid grid; // spatial hash grid for efficient neighbourhood search
	int num_particles;

	double particle_mass; // Mass of each particle
	double rho; // Rest density
	//double gravity_f; // Force of gravity on the system
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
	
	/*
        Initialize a blank new fluid.

        Input
                int num_particles: the number of particles in the fluid system.
                double particle_mass: the mass of each fluid particle.
                double rho: the rest density of the fluid system.

                double gravity_f: the externel force of gravity.
                double user_f: the externel force supplied by user.

                int jacobi_iterations: the number of steps in the jacobi constraint solver.
                double cfm_epsilon: constraint force mixing term for regularization.
                double kernel_h: the width of the kernels.

                double tensile_k: parameter for artificial pressure. See equation [13]
                double tensile_delta_q: parameter for artificial pressure. See equation [13]
                int tensile_n: parameter for artificial pressure. See equation [13]

                double viscocity_c: weight of XSPH viscocity. See equation [17]
                double vorticity_epsilon: parameter for vorticity confinement. See equation [16]

                double lower_bound: the lower bounding box coordinate
                double upper_bound: the upper bounding box coordinate

                double dt: time step delta


        */
	Fluid(int num_particles, double particle_mass, double rho, double gravity_f, double user_f, int jacobi_iterations, 
			double cfm_epsilon, double kernel_h, double tensile_k, double tensile_delta_q, int tensile_n, 
			double viscocity_c, double vorticity_epsilon, double lower_bound, double upper_bound, double dt);

        /*
        Given the fluid state global matrix (shared by visualization), initialize the particles
        of the system, and the spatial hash grid.
        
        Input:

                Eigen::MatrixXd &fluid_state: pointer the global fluid state matrix
                                which stores the position of all the particles in the system
                                
                                ex. fluid_state(i, 0) = x coordinate of particle i. 
        */
	void init_state(Eigen::MatrixXd &fluid_state);
	
        /*
        Take step in the fluid system, and update the global fluid state particle positions.

        Input: 
                Eigen::MatrixXd &fluid_state: pointer the global fluid state matrix we modify.
                Eigen::MatrixXd &colors: FOR DEBUGGING neighbhours.
        */	
	void step(Eigen::MatrixXd &fluid_state, Eigen::MatrixXd &colors);

};

#endif
