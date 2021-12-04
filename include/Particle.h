#ifndef PARTICLE
#define PARTICLE

#include <Eigen/Dense>
#include <set>
class Particle{
public:
	// Basic State Parameters
	double m; 
	double rho = 0;
	Eigen::Vector3d x; // World Space Position
	Eigen::Vector3d v; // Velocity
	Eigen::Vector3d f; // Forces

	// Neighourhood Parameters
        int global_idx; // each particle stores a global index into the fluid
	std::set<int> neighbours; // contains global indices of neighbouring particles 

	// Jacobi Parameters
	double c; // density constraint result
	double lambda; // magnitude of force solving constraint gradient
	double c_grad_neighorhood_norm; // 2-norm of constraint gradient 
	Eigen::Vector3d dP; // total position update including corrections from density constraints
	Eigen::Vector3d x_new; // updated position
	Eigen::Vector3d v_new; // updated velocity

	// Vorticity Parameters
	Eigen::Vector3d omega; // the curl at the particles position
	Eigen::Vector3d eta; // differential operator acting on omega
	Eigen::Vector3d N; // Normalized location vector
	Eigen::Vector3d vorticity_f; // corrective vorticity force

        // Spatial Grid Hashing (for removing from Spatial-Hash)
        Eigen::Vector3d cell_coord;

public:
	/*
	Initialize a new fluid particle.

        Input:
                Eigen::Vector3d x_init: the initial position of the particle in the scene.
                double m: the mass of the particle.
                int global_idx: unique global index of the particle (indexes into shared fluid_state)
	*/
	Particle(Eigen::Vector3d X_init, double m, int global_idx);
};

#endif
