#include <Fluid.h>
#include <kernel.h>
#include <vorticity.h>
#include <viscocity.h>

#include <iostream>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

#define DEBUG 1

Fluid::Fluid(int num_particles, double particle_mass, double rho, double gravity_f, double user_f, int jacobi_iterations, 
			double cfm_epsilon, double kernel_h, double tensile_k, double tensile_delta_q, int tensile_n, 
			double viscocity_c, double vorticity_epsilon, double lower_bound, double upper_bound, double dt){

        this->num_particles = num_particles;
	this->particle_mass = particle_mass;
	this->rho = rho;

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

        this->grid = SpatialHashGrid(lower_bound, upper_bound, kernel_h);


        // testing vectorization
        this->x_new.resize(num_particles, 3);
        this->v.resize(num_particles, 3);
        this->dP.resize(num_particles, 3);
        this->omega.resize(num_particles, 3);
        this->eta.resize(num_particles, 3);
        this->N.resize(num_particles, 3);
        this->vorticity_f.resize(num_particles, 3);
        this->cell_coord.resize(num_particles, 3);

        this->density.resize(num_particles);
        this->c.resize(num_particles);
        this->lambda.resize(num_particles);
        this->c_grad_norm.resize(num_particles);

        this->gravity_f = Eigen::VectorXd::Constant(num_particles, gravity_f);

        std::vector<std::vector<int>> neighbours_init(num_particles);
        this->neighbours = neighbours_init;

}	

void Fluid::init_state(Eigen::MatrixXd &fluid_state){
        this->v.setZero();
        this->grid.insert(fluid_state);
}

void Fluid::step(Eigen::MatrixXd &fluid_state, Eigen::MatrixXd &colors){
        auto t0 = Clock::now();

	// Apply External forces
        this->v.col(1) -= this->particle_mass * this->gravity_f;
        this->x_new = fluid_state + this->dt * this->v;

        auto t1 = Clock::now();
        if (DEBUG) std::cout << "\n------------------------------------------\nApplied Externel forces [" << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count() << " s]\n";

        // Get neighbours using spatial hash grid
        this->grid.findNeighbours(this->x_new, this->neighbours);

        auto t2 = Clock::now();
        if (DEBUG) std::cout << "Found Neighbours [" << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " s]\n";

	//Jacobi Loop 
	for(int i = 0; i < this->jacobi_iterations; i++){
                // reset jacobi state
                this->density.setZero();
                this->lambda.setZero();
                this->c_grad_norm.setZero();

                Eigen::Vector3d c_grad_temp;
                
                for(int p_i = 0; p_i < this->num_particles; p_i++){
                        
                        // compute densities
                        
                        for(int p_j : this->neighbours[i]){
                                this->density[p_i] += this->particle_mass * kernel_poly6(this->x_new.row(p_i), this->x_new.row(p_j), this->kernel_h);

                                // accumulate gradient norm
                                c_grad_temp.setZero();
                                if (p_i == p_j){
                                        for(int p_k : this->neighbours[p_i]){
                                               kernel_spiky(this->ker_res, this->x_new.row(p_i), this->x_new.row(p_k), this->kernel_h);
                                               c_grad_temp += this->ker_res;
                                        }
                                }
                                else{
                                        kernel_spiky(this->ker_res, this->x_new.row(p_i), this->x_new.row(p_j), this->kernel_h);
                                        c_grad_temp = this->ker_res;
                                }

                                this->c_grad_norm[i] += (c_grad_temp / this->rho).norm();
                        }

                        // Compute constraint and lambda
                        this->c[p_i] = (this->density[p_i] / this->rho) - 1;
                        this->lambda[p_i] = -this->c[p_i] / (this->c_grad_norm[p_i] + this->cfm_epsilon);
                }

                auto t3 = Clock::now();
                if (DEBUG) std::cout << "Computed Constraints [" << std::chrono::duration_cast<std::chrono::seconds>(t3 - t2).count() << " s]\n";

		// Compute dP
                this->dP.setZero();
                for(int p_i = 0; p_i < this->num_particles; p_i++){
                        for(int p_j : this->neighbours[i]){
                                kernel_spiky(this->ker_res, this->x_new.row(p_i), this->x_new.row(p_j), this->kernel_h);
                                // precompute this denominator

		                this->s_corr = - this->tensile_k * pow(kernel_poly6(this->x_new.row(p_i), this->x_new.row(p_j), this->kernel_h) / kernel_poly6(this->tensile_delta_q, this->kernel_h), this->tensile_n);
                                this->dP.row(p_i) += (this->lambda[p_i] + this->lambda[p_j] + this->s_corr) * this->ker_res; 
                        }
                        this->dP.row(p_i) /= this->rho;
                }

                auto t4 = Clock::now();
                if (DEBUG) std::cout << "Computed Position Correction [" << std::chrono::duration_cast<std::chrono::seconds>(t4 - t3).count() << " s]\n";

		// Collision Detection with boundary box and solid
                for(int p_i = 0; p_i < this->num_particles; p_i++){
			this->x_new.row(p_i) += 0.005 * this->dP.row(p_i); // TODO CALIBRATE SIMULATION
                        
			//Collision Detect correct p.x_new (naive)
			if (this->x_new.row(p_i)[0] < this->lower_bound){ 
                                this->x_new.row(p_i)[0] = this->lower_bound;
                                if (this->v(p_i, 0) < 0) this->v(p_i, 0) *= -1;
                        }
			if (this->x_new.row(p_i)[0] > this->upper_bound){ 
                                this->x_new.row(p_i)[0] = this->upper_bound;
                                if (this->v(p_i, 0) > 0) this->v(p_i, 0) *= -1;
                        }
			if (this->x_new.row(p_i)[1] < this->lower_bound){ 
                                this->x_new.row(p_i)[1] = this->lower_bound;
                                if (this->v(p_i, 1) < 0) this->v(p_i, 1) *= -1;
                        }
			if (this->x_new.row(p_i)[1] > this->upper_bound){ 
                                this->x_new.row(p_i)[1] = this->upper_bound;
                                if (this->v(p_i, 1) > 0) this->v(p_i, 1) *= -1;
                        }
			if (this->x_new.row(p_i)[2] < this->lower_bound){ 
                                this->x_new.row(p_i)[2] = this->lower_bound;
                                if (this->v(p_i, 2) < 0) this->v(p_i, 2) *= -1;
                        }
			if (this->x_new.row(p_i)[2] > this->upper_bound){ 
                                this->x_new.row(p_i)[2] = this->upper_bound;
                                if (this->v(p_i, 2) > 0) this->v(p_i, 2) *= -1;
                        }
                }
                auto t5 = Clock::now();
                if (DEBUG) std::cout << "Collision Detection [" << std::chrono::duration_cast<std::chrono::seconds>(t5 - t4).count() << " s]\n";
	}

	//Update Velocity
        this->v = (this->x_new - fluid_state) / this->dt;

	// Vorticity (O(n^2))
	// apply_vorticity(this->fluid, this->kernel_h, this->vorticity_epsilon, this->dt);
	// apply_viscocity(this->fluid, this->kernel_h, this->viscocity_c);


	// Update Position and spatial hash grid
        fluid_state = this->x_new;
        this->grid.update(fluid_state);


        auto t6 = Clock::now();
        if (DEBUG) std::cout << "Simulation Step Total Time [" << std::chrono::duration_cast<std::chrono::seconds>(t6 - t0).count() << " s]\n----------------------------------------\n";

        // Debugging using colors for now.
        for (int i = 0; i < num_particles; i++){
                colors.row(i) << 0, 0, 1; 
        }
        colors.row(0) << 1, 0, 0; // track particle 0 in red
        for (auto n_idx : this->neighbours[0]){
                if (n_idx != 0){
                        colors.row(n_idx) << 0, 1, 0; // track neighbours in green
                }
        }

	this->t += this->dt;
}
