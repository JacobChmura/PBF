#ifndef FLUID
#define FLUID

#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <SpatialHashGrid.h>
#include <Eigen/Dense>

class Fluid{
public:
	//For Density Plot UI
	double t;
	double avg_density;
	double max_density;

private:	

	/*
	All matrices are of size (num_particles, 3)
		-> first axis indexes into the particle of the fluid
		-> second axis is (x, y, z)
	
	All vectors are of size (num_particles, )
	*/

	//Basic Parameters
	int num_particles;
	double dt;
	double particle_mass; 
	double rho; // Rest density
	double kernel_h; // Kernel radius
	double gravity_f_; // Force of gravity on the system
	double user_f_; // Force of user applied force on the system

	// Vectorized copy, num particles unknown at compile time
	Eigen::VectorXd gravity_f;
	Eigen::VectorXd user_f;

	// Global State Parameters
        Eigen::MatrixXd x_new; // position prediction updated on every Jacobi iteration
        Eigen::MatrixXd v_new; // velocity prediction
        Eigen::MatrixXd v; // fluid velocity matrix
        
        // Jacobi Loop Parameters
        int jacobi_iterations;
        Eigen::MatrixXd dP; // position correction
        Eigen::VectorXd density; // fluid density
        Eigen::VectorXd c; // constraint equation
        Eigen::VectorXd lambda; // constraint lagrange
        Eigen::VectorXd c_grad_norm; // constraint grad norm

        // Neighorhood Parameters
        SpatialHashGrid grid; // spatial hash grid for efficient neighbourhood search
        std::vector<std::vector<int>> neighbours; // Lookup table for neighbours. Ex. neighbours[i] = list of particles neighouring particle i

	// Constraint Force Mixing Parameter Equation [11] https://mmacklin.com/pbf_sig_preprint.pdf
	double cfm_epsilon;

	// Tensile Stability Parameters Equation [13] https://mmacklin.com/pbf_sig_preprint.pdf (artifical pressure term)
	int tensile_n;
	double tensile_k;
	double tensile_delta_q;
	double s_corr; 
        double tensile_stability_denom; 

	// // Viscocity and Vorticity Parameters Equation [16, 17] https://mmacklin.com/pbf_sig_preprint.pdf
	double viscocity_c;
	double vorticity_epsilon;
	Eigen::MatrixXd omega;
	Eigen::MatrixXd eta;
	Eigen::MatrixXd N;
	Eigen::MatrixXd vorticity_f;

	// Bounding box for hash grid and collision
	double lower_bound;
	double upper_bound;

	// Temp variables
	Eigen::Vector3d ker_res; // stores spiky kernel 3-vector 
        Eigen::Vector3d c_grad_temp; // stores cumulative gradient magnitude
        Eigen::MatrixXd user_force_direction; // stores user force direction for external force 


public:
	
	/*
        Initialize a blank new fluid.

        Input
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
	Fluid(double particle_mass, double rho, double gravity_f, double user_f, int jacobi_iterations, 
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
                Eigen::MatrixXd &colors: for showing neighbours.
		
		Eigen::Vector3d mouse_pos in world coordinates
		bool add_user_force: if true, use the mouse_pos to compute an externel force in the direction of mouse_pos on the system
		
		bool use_viscocity: whether or not to use XPSH viscocity (can change on keydown)
		bool use_vorticity: whether or not to use Vorticity confinement (can change on keydown)

        */	
	void step(Eigen::MatrixXd &fluid_state, Eigen::MatrixXd &colors, Eigen::Vector3d mouse_pos, bool add_user_force, bool use_viscocity, bool use_vorticity);

};

#endif
