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

        this->grid = SpatialHashGrid(lower_bound, upper_bound, kernel_h, num_particles);
}	

void Fluid::init_state(Eigen::MatrixXd &fluid_state){
	this->fluid_state = fluid_state;
              
	for(int i = 0; i < this->num_particles; i++){
		Particle p(fluid_state.row(i), this->particle_mass, i); // create particles
		this->fluid.push_back(p); // push into global fluid state
                this->grid.insert(p); // add to hash grid 
	}
}

void Fluid::step(Eigen::MatrixXd &fluid_state, Eigen::MatrixXd &colors){

	auto accumulate_c_grad_norm = [&](Particle &p_i, Particle &p_j){
		Eigen::Vector3d c_grad_norm_ij;
		c_grad_norm_ij.setZero();

		if (p_i.global_idx == p_j.global_idx){
			for (int particle_idx : p_i.neighbours){
                                Particle p_k = this->fluid[particle_idx];
				kernel_spiky(this->ker_res, p_i, p_k, this->kernel_h);
				c_grad_norm_ij += this->ker_res;
			}
		}
		else{
			kernel_spiky(this->ker_res, p_i, p_j, this->kernel_h);
			c_grad_norm_ij = - this->ker_res;
		}
                c_grad_norm_ij /= this->rho;
		p_i.c_grad_neighorhood_norm += c_grad_norm_ij.norm();
	};

	auto compute_s_corr = [&](double &s_corr, Particle &p_i, Particle &p_j){
		s_corr = - this->tensile_k * pow(kernel_poly6(p_i, p_j, this->kernel_h) / kernel_poly6(this->tensile_delta_q, this->kernel_h), this->tensile_n);
	};


	// Apply External forces
        //Eigen::Vector3d vel = this->fluid[0].v;
        //Eigen::Vector3d x_pos = this->fluid[0].x;
        //Eigen::Vector3d x_pos_new = this->fluid[0].x_new;
        //std::cout << "Particle Velocity before Update: " << vel[1] << "\n";
        //std::cout << "Particle x_new before Update: " << x_pos_new[1] << "\n"; 
        //std::cout << "Particle x before Update: " << x_pos[1] << "\n"; 

	for(Particle &p : this->fluid){
		p.f.setZero();
		
		// Gravity Force
		p.f[1] = - p.m * this->gravity_f;

		// User Force (TODO)
		p.v += this->dt * p.f;
		p.x_new = p.x + this->dt * p.v; 
	}
	std::cout << "\tApplied Externel forces.\n";
        //vel = this->fluid[0].v;
        //x_pos = this->fluid[0].x;
        //x_pos_new = this->fluid[0].x_new;
        //std::cout << "\nParticle Velocity After Gravity: " << vel[1] << "\n"; 
        //std::cout << "Particle x_new After Gravity: " << x_pos_new[1] << "\n"; 
        //std::cout << "Particle x After Gravity: " << x_pos[1] << "\n"; 

        // Get neighbours using spatial hash grid
        for(Particle &p: this->fluid){
                this->grid.findNeighbours(p);
        }
	std::cout << "\tFound neighbors.\n";
        
	//Jacobi Loop (each component should be parrellized)
	for(int i = 0; i < this->jacobi_iterations; i++){
		// Compute densities
		for(Particle &p_i : this->fluid){
			p_i.rho = 0;
			p_i.lambda = 0;
			p_i.c_grad_neighorhood_norm = 0;

			// Compute density estimate 
			for(int particle_idx: p_i.neighbours){
                                Particle p_j = this->fluid[particle_idx];
				p_i.rho += p_i.m * kernel_poly6(p_i, p_j, this->kernel_h);
				accumulate_c_grad_norm(p_i, p_j);
			}
			// Compute constraint lambda
			p_i.c = (p_i.rho / this->rho) - 1;
			p_i.lambda = - p_i.c / (p_i.c_grad_neighorhood_norm + this->cfm_epsilon);
		}

                Particle p_first = this->fluid[0];
                
                //std::cout << "Rest Density: " << this->rho << " |  Particle Density: " << p_first.rho << " |  Constraint : " << p_first.c << " | Num neighbours: " << p_first.neighbours.size() << "\n"; 
                //for (int particle_idx: p_first.neighbours){
                //        Particle p_first_neighbour = this->fluid[particle_idx];
                //        std::cout << "Neighbhour " << p_first_neighbour.global_idx << " |  Particle Density: " << p_first_neighbour.rho << " |  Constraint : " << p_first_neighbour.c << "\n"; 

                //}

		std::cout << "\t\tComputed Densities\n";
		// Compute dP
		for(Particle &p_i : this->fluid){
			p_i.dP.setZero();

			for(int particle_idx : p_i.neighbours){
                                Particle p_j = this->fluid[particle_idx];
				kernel_spiky(this->ker_res, p_i, p_j, this->kernel_h);
				compute_s_corr(this->s_corr, p_i, p_j);
				p_i.dP += (p_i.lambda + p_j.lambda + this->s_corr) * this->ker_res;
			}
                        p_i.dP /= this->rho;
		}
                p_first = this->fluid[0];
                ///std::cout << "lambda: " << p_first.lambda << "\ndP norm: " << p_first.dP.norm() << " \n";
                ///for (int particle_idx: p_first.neighbours){
                ///        Particle p_first_neighbour = this->fluid[particle_idx];
                ///        std::cout << "lambda: " << p_first_neighbour.lambda << "\ndP norm: " << p_first_neighbour.dP.norm() << " \n";

                ///}
                //exit(0);
		std::cout << "\t\tComputed dP\n";
		// Collision Detection with boundary box and solid (TODO)
		for (Particle &p : this->fluid){
			p.x_new += 0.001 * p.dP;
                        
                        double damp = 1;
			//Collision Detect correct p.x_new (naive)
			if (p.x_new[0] < this->lower_bound){ 
                                p.x_new[0] = this->lower_bound;
                                if (p.v[0] < 0) p.v[0] *= -damp;

                        }
                        if (p.x_new[0] > this->upper_bound){ 
                                p.x_new[0] = this->upper_bound;
                                if (p.v[0] > 0) p.v[0] *= -damp;
                        }
                        if (p.x_new[1] < this->lower_bound){
                                p.x_new[1] = this->lower_bound;
                                if (p.v[1] < 0) p.v[1] *= -damp;
                        }
			if (p.x_new[1] > this->upper_bound){
                                p.x_new[1] = this->upper_bound;
                                if (p.v[1] > 0) p.v[1] *= -damp;
                        }
                        if (p.x_new[2] < this->lower_bound){
                                p.x_new[2] = this->lower_bound;
                                if (p.v[2] < 0) p.v[2] *= -damp;
                        }
			if (p.x_new[2] > this->upper_bound){
                                p.x_new[2] = this->upper_bound;
                                if (p.v[2] > 0) p.v[2] *= -damp;
                        }
		}
		std::cout << "\t\tCollision Detection\n";
	}

        //vel = this->fluid[0].v;
        //x_pos = this->fluid[0].x;
        //x_pos_new = this->fluid[0].x_new;
        //std::cout << "\nParticle Velocity After collision: " << vel[1] << "\n"; 
        //std::cout << "Particle x_new After collision: " << x_pos_new[1] << "\n"; 
        //std::cout << "Particle x After collision: " << x_pos[1] << "\n"; 
	//Update Velocity
	for(Particle &p: this->fluid){
		p.v = (p.x_new - p.x) / this->dt;
	}

        //vel = this->fluid[0].v;
        //x_pos = this->fluid[0].x;
        //x_pos_new = this->fluid[0].x_new;
        //std::cout << "\nParticle Velocity After update: " << vel[1] << "\n"; 
        //std::cout << "Particle x_new After update: " << x_pos_new[1] << "\n"; 
        //std::cout << "Particle x After update: " << x_pos[1] << "\n\n"; 
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

        for (int i = 0; i < num_particles; i++){
                colors.row(i) << 0, 0, 1;
        }
        colors.row(0) << 1, 0, 0; // track particle 0 in red
        for (auto n_idx : this->fluid[0].neighbours){
                if (n_idx != 0){
                        colors.row(n_idx) << 0, 1, 0; // track neighbours in green
                }
        }

	this->t += this->dt;

	//std::cout << "\tUpdate\n";
}
