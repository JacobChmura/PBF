#include <Fluid.h>
#include <kernel_poly6.h>
#include <kernel_spiky.h>
#include <vorticity.h>
#include <viscocity.h>

#include <iostream>

Fluid::Fluid(int num_particles, double particle_mass, double rho, double gravity_f, double user_f, int jacobi_iterations, 
			double cfm_epsilon, double kernel_h, double tensile_k, double tensile_delta_q, int tensile_n, 
			double viscocity_c, double vorticity_epsilon, double lower_bound, double upper_bound, double dt){

        this->num_particles = num_particles;
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

	this->lower_bound = lower_bound;
	this->upper_bound = upper_bound;

        this->grid = SpatialHashGrid(lower_bound, 2 * kernel_h, num_particles);
}	

void Fluid::init_state(Eigen::MatrixXd &fluid_state){
	this->fluid_state = fluid_state;
              
	for(int i = 0; i < this->num_particles; i++){
		Particle p(fluid_state.row(i), this->particle_mass, i); // create particles
		this->fluid.push_back(p); // push into global fluid state
                this->grid.insert(p); // add to hash grid 
	}
}

void Fluid::step(Eigen::MatrixXd &fluid_state){

	auto accumulate_c_grad_norm = [&](Particle &p_i, Particle &p_j){
		Eigen::Vector3d c_grad_norm_ij;
		c_grad_norm_ij.setZero();

		if (&p_i == &p_j){
			for (int particle_idx : p_i.neighbours){
                                Particle p_k = this->fluid[particle_idx];
				kernel_spiky(this->ker_res, p_i, p_k, this->kernel_h);
				c_grad_norm_ij += (1 / this->rho) * this->ker_res;
			}
		}
		else{
			kernel_spiky(this->ker_res, p_i, p_j, this->kernel_h);
			c_grad_norm_ij = (-1 / this->rho) * this->ker_res;
		}

		p_i.c_grad_neighorhood_norm += c_grad_norm_ij.norm();
	};

	auto compute_s_corr = [&](double &s_corr, Particle &p_i, Particle &p_j){
		s_corr = - this->tensile_k * pow(kernel_poly6(p_i, p_j, this->kernel_h) / kernel_poly6(this->tensile_delta_q, this->kernel_h), this->tensile_n);
	};


	// Apply External forces
	for(Particle &p : this->fluid){
		p.f.setZero();
		
		// Gravity Force
		p.f[1] = - p.m * this->gravity_f;

		// User Force (TODO)
		p.v += this->dt * p.f;
		p.x_new = p.x + this->dt * p.v; 
	}
	std::cout << "\tApplied Externel forces.\n";

        // Get neighbours using spatial hash grid
        for(Particle &p: this->fluid){
                this->grid.findNeighbours(p);
        }
	std::cout << "\tFound neighbors.\n";
        
	//Jacobi Loop (each component should be parrellized)
	for(int i = 0; i < this->jacobi_iterations; i++){
		std::cout << "\tJacobi Loop:\n";
		// Compute densities
		for(Particle &p_i : this->fluid){
			p_i.rho = 0;
			p_i.lambda = 0;
			p_i.c_grad_neighorhood_norm = 0;

			// Compute density estimate 
                        //std::cout << "num neibhours: " << p_i.neighbours.size() << "\n";
			for(int particle_idx: p_i.neighbours){
                                Particle p_j = this->fluid[particle_idx];
				p_i.rho += p_i.m * kernel_poly6(p_i, p_j, this->kernel_h);
				accumulate_c_grad_norm(p_i, p_j);
			}

			// Compute constraint lambda
			p_i.c = (p_i.rho / this->rho) - 1;
			//std::cout << p_i.c << "\n";
			p_i.lambda = - p_i.c / (p_i.c_grad_neighorhood_norm + this->cfm_epsilon);
		}

		std::cout << "\t\tComputed Densities\n";
		// Compute dP
		for(Particle &p_i : this->fluid){
			p_i.dP.setZero();

			for(int particle_idx : p_i.neighbours){
                                Particle p_j = this->fluid[particle_idx];
				kernel_spiky(this->ker_res, p_i, p_j, this->kernel_h);
				compute_s_corr(this->s_corr, p_i, p_j);
				p_i.dP += ((p_i.lambda + p_j.lambda + this->s_corr) / this->rho) * this->ker_res;
			}
		}
		std::cout << "\t\tComputed dP\n";
		// Collision Detection with boundary box and solid (TODO)
		for (Particle &p : this->fluid){
			p.x_new += p.dP;

			// Collision Detect correct p.x_new (naive)
			if (p.x_new[0] < this->lower_bound) p.x_new[0] = this->lower_bound;
			if (p.x_new[0] > this->upper_bound) p.x_new[0] = this->upper_bound;
			if (p.x_new[1] < this->lower_bound) p.x_new[1] = this->lower_bound;
			if (p.x_new[1] > this->upper_bound) p.x_new[1] = this->upper_bound;
			if (p.x_new[2] < this->lower_bound) p.x_new[2] = this->lower_bound;
			if (p.x_new[2] > this->upper_bound) p.x_new[2] = this->upper_bound;
		}
		std::cout << "\t\tCollision\n";
	}

	// Update Velocity
	for(Particle &p: this->fluid){
		p.dX = p.x_new - p.x;
		p.v += p.dX / this->dt;
	}

	// Vorticity (O(n^2))
	// apply_vorticity(this->fluid, this->kernel_h, this->vorticity_epsilon, this->dt);

	//	apply_viscocity(this->fluid, this->kernel_h, this->viscocity_c);


	// Update Position and spatial hash grid
	for (Particle &p : this->fluid){
		p.x = p.x_new;
                this->grid.update(p);
	}

	// Write to global state
	for (int i = 0; i < this->num_particles; i++){
		fluid_state.row(i) = (this->fluid)[i].x;
	}
	this->t += this->dt;

	std::cout << "\tUpdate\n";
}
