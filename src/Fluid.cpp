#include <Fluid.h>
#include <kernel_poly6.h>
#include <kernel_spiky.h>
#include <vorticity.h>
#include <visocity.h>

Fluid::Fluid(double particle_mass, double rho, double gravity_f, double user_f, int jacobi_iterations, 
			double cfm_epsilon, double kernel_h, double tensile_k, double tensile_delta_q, int tensile_n, 
			double viscocity_c, double vorticity_epsilon, double dt){

	this->particle_mass = particle_mass;
	this->rho = rho;

	this->gravity_f = gravity_f;
	this->user_f = user_f;

	this->jacobi_iterations = jacobi_iterations;

	this->cfm_epsilon = cfm_epsilon;
	this->kernel_h = kernel_h;

	this->tensile_k = tensile_k;
	this->tensile_delta_q = tensile_delta_q;
	this->tensile_n = tensile_n;

	this->viscocity_c = viscocity_c;
	this->vorticity_epsilon = vorticity_epsilon;

	this->t = 0;
	this->dt = dt;
}	

Fluid::init_state(){}

Fluid::step(){

	auto accumulate_c_grad_norm = [&](Particle &p_i, Particle &p_j){
		Eigen::Vector3d c_grad_norm_ij;
		c_grad_norm_ij.setZero();

		if (p_i == p_j){
			for (Particle &p_k : p_i.neighbors){
				c_grad_norm_ij += (1 / this->rho) * kernel_spiky(p_i, p_k, this->kernel_h);
			}
		}
		else{
			c_grad_norm_ij = (-1 / this->rho) * kernel_spiky(p_i, p_j, this->kernel_h);
		}

		p_i.c_grad_neighorhood_norm += c_grad_norm_ij.norm();
	};

	auto compute_s_corr = [&](Particle &p_i, Particle &p_j){
		return - this->tensile_k * pow(kernel_poly6(p_i, p_j, this->kernel_h) / kernel_poly6(this->tensile_delta_q), this->tensile_n);
	};


	// Apply External forces
	for(Particle &p : this->fluid){
		p.f.setZero();
		
		// Gravity Force
		p.f[1] = - p.m * this->gravity_f;

		// User Force (TODO)

		p.v += this->dt * p.f;
		p.x_new += this->dt * p.v; 
	}

	// Get Neighours (Naive O(n^2))
	for(Particle &p_i : this->fluid){
		p.neighbors.clear();

		for(Particle &p_j : this->fluid){
			if ((p_i.x_new - p_j.x_new).norm() <= 2 * this->kernel_h) p_i.neighbors.push_back(p_j);
		}
	}

	// Jacobi Loop (each component should be parrellized)
	for(int i = 0; i < this->jacobi_iterations; i++){

		// Compute densities
		for(Particle &p_i : this->fluid){
			p_i.rho = 0;
			p_i.lambda = 0;
			p_i.c_grad_neighorhood_norm = 0;

			// Compute density estimate 
			for(Particle &p_j : p_i.neighbors){
				p_i.rho += p_i.m * kernel_poly6(p_i, p_j, this->kernel_h);
				accumulate_c_grad_norm(p_i, p_j);
			}

			// Compute constraint lambda
			p_i.c = (p_i.rho / this->rho) - 1;
			p_i.lambda = - p_i.c / (p_i.c_grad_neighorhood_norm + this->cfm_epsilon);
		}

		// Compute dP
		for(Particle &p_i : this->fluid){
			p_i.dP.setZero();

			for(Particle &p_j : p_j.neighbors){
				p_i.dP += ((p_i.lambda + p_j.lambda + compute_s_corr(p_i, p_j)) / this->rho) * kernel_spiky(p_i, p_j);
			}
		}

		// Collision Detection with boundary box and solid (TODO)
		for (Particle &p : this->fluid){
			p.x_new += p.dP;

			// Collision Detect correct p.x_new
		}
	}

	// Update Velocity
	for(Particle &p: this->fluid){
		p_i.dX = p_i.x_new - p_i.x;
		p_i.v += p_i.dX / this->dt;
	}

	// Vorticity (O(n^2))
	apply_vorticity(this->fluid);

	// Viscocity (O(n^2))
	apply_viscocity(this->fluid);

	// Update Position
	for (Particle &p : this->fluid){
		p_i.x = p_i.x_new;
	}

	this->t += this->dt;
}