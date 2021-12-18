#include <Fluid.h>
#include <kernel.h>
#include <vorticity.h>
#include <viscocity.h>


#include <iostream>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

#define DEBUG 0

Fluid::Fluid(double particle_mass, double rho, double gravity_f, double user_f, int jacobi_iterations, 
			double cfm_epsilon, double kernel_h, double tensile_k, double tensile_delta_q, int tensile_n, 
			double viscocity_c, double vorticity_epsilon, double lower_bound, double upper_bound, double dt){

        
	this->particle_mass = particle_mass;
	this->rho = rho;

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

        this->tensile_stability_denom = pow(kernel_poly6(tensile_delta_q, kernel_h), tensile_n);

        this->gravity_f_ = gravity_f;
        this->user_f_ = user_f;

}	

void Fluid::init_state(Eigen::MatrixXd &fluid_state){
        int num_particles = fluid_state.rows();
        this->num_particles = num_particles;
        this->x_new.resize(num_particles, 3);
        this->v.resize(num_particles, 3);
        this->v_new.resize(num_particles, 3);
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

        this->gravity_f = Eigen::VectorXd::Constant(num_particles, this->gravity_f_);
        this->user_f = Eigen::VectorXd::Constant(num_particles, this->user_f_);

        std::vector<std::vector<int>> neighbours_init(num_particles);
        this->neighbours = neighbours_init;

        v.setZero();
        grid.update(fluid_state);
}

void Fluid::step(Eigen::MatrixXd &fluid_state, Eigen::MatrixXd &colors, Eigen::Vector3d mouse_pos, bool add_user_force, bool use_viscocity, bool use_vorticity){
        auto t0 = Clock::now();

	// Apply External forces
        v.col(1) -= particle_mass * gravity_f;
                   
        // Apply User Foces if applicable
        if (add_user_force) {
                Eigen::MatrixXd test = (-fluid_state).rowwise() + mouse_pos.transpose();
                test.rowwise().normalize();
                v += particle_mass * user_f(0) * test;
        }
        // Predict next time step position
        x_new = fluid_state + dt * v;

        auto t1 = Clock::now();
        if (DEBUG) std::cout << "\n------------------------------------------\nApplied Externel forces [" << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms]\n";

        // Get neighbours using spatial hash grid
        grid.findNeighbours(x_new, neighbours);

        auto t2 = Clock::now();
        if (DEBUG) std::cout << "Found Neighbours [" << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms]\n";

	//Jacobi Loop 
	for(int i = 0; i < jacobi_iterations; i++){
                // reset jacobi state
                density.setZero();
                lambda.setZero();
                c_grad_norm.setZero();

                for(int p_i = 0; p_i < num_particles; p_i++){
                        
                        // compute densities
                        for(int p_j : neighbours[p_i]){
                                density[p_i] += particle_mass * kernel_poly6(x_new.row(p_i), x_new.row(p_j), kernel_h);

                                // accumulate gradient norm
                                c_grad_temp.setZero();
                                if (p_i == p_j){ 
                                        for(auto p_k : neighbours[p_i]){
                                               kernel_spiky(ker_res, x_new.row(p_i), x_new.row(p_k), kernel_h);
                                               c_grad_temp += ker_res;
                                        }
                                }
                                else{
                                        kernel_spiky(ker_res, x_new.row(p_i), x_new.row(p_j), kernel_h);
                                        c_grad_temp = ker_res;
                                }

                                c_grad_norm[p_i] += (c_grad_temp / rho).norm();
                        }

                        // Compute constraint and lambda
                        c[p_i] = (density[p_i] / rho) - 1;
                        lambda[p_i] = -c[p_i] / (c_grad_norm[p_i] + cfm_epsilon);
                }
                
                auto t3 = Clock::now();
                if (DEBUG) std::cout << "Computed Constraints [" << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() << " ms]\n";

                //std::cout << "Density (jacobi iteration " << i << ") : " << density.sum() / density.size() << std::endl;

		// //Compute dP
                dP.setZero();
                for(int p_i = 0; p_i < num_particles; p_i++){
                        for(int p_j : neighbours[p_i]){
                                kernel_spiky(ker_res, x_new.row(p_i), x_new.row(p_j), kernel_h);
		                s_corr = - tensile_k * pow(kernel_poly6(x_new.row(p_i), x_new.row(p_j), kernel_h), tensile_n) / tensile_stability_denom;           
                                dP.row(p_i) += (lambda[p_i] + lambda[p_j] + s_corr) * ker_res; 
                        }
                        dP.row(p_i) /= rho;
                }
                x_new += 0.005 * dP; // Calibrate!
                auto t4 = Clock::now();
                if (DEBUG) std::cout << "Computed Position Correction [" << std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count() << " ms]\n";

		// Collision Detection with boundary box
                for(int p_i = 0; p_i < num_particles; p_i++){

			//Collision Detect correct p.x_new (naive)
                        for(int axis = 0; axis < 3; axis++){
                                if (x_new.row(p_i)[axis] < lower_bound){ 
                                        x_new.row(p_i)[axis] = lower_bound;
                                        if (v(p_i, axis) < 0) v(p_i, axis) *= -3;
                                }
                                if (x_new.row(p_i)[axis] > upper_bound){ 
                                        x_new.row(p_i)[axis] = upper_bound;
                                        if (v(p_i, axis) > 0) v(p_i, axis) *= -3;
                                }
                        }
                }



                auto t5 = Clock::now();
                if (DEBUG) std::cout << "Collision Detection [" << std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4).count() << " ms]\n";
	}

	//Update Velocity
        v = (x_new - fluid_state) / dt;
        
        // Vorticity and Viscocity
        auto t6 = Clock::now();
        if (use_vorticity) apply_vorticity(x_new, neighbours, v, omega, eta, N, vorticity_f, vorticity_epsilon, kernel_h, dt);

        auto t7 = Clock::now();
        if (DEBUG) std::cout << "Vorticity [" << std::chrono::duration_cast<std::chrono::milliseconds>(t7 - t6).count() << " ms]\n";


        if (use_viscocity) apply_viscocity(x_new, neighbours, v, v_new, viscocity_c, kernel_h);
        auto t8 = Clock::now();
        if (DEBUG) std::cout << "Viscocity [" << std::chrono::duration_cast<std::chrono::milliseconds>(t8 - t7).count() << " ms]\n";


	// Update Position and spatial hash grid
        fluid_state = x_new;
        grid.update(fluid_state);

        auto t9 = Clock::now();
        if (DEBUG) std::cout << "Simulation Step Total Time [" << std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t0).count() << " ms]\n----------------------------------------\n";

        // Update Density results for plotting
        avg_density = density.sum() / density.size();
        max_density = density.maxCoeff();

        // // Debugging using colors for now.
        // colors.row(0) << 1, 0, 0; // track particle 0 in red
        // for (int i = 1; i < num_particles; i++){
        //         colors.row(i) << 0, 0, 1; 
        // }
        // for (auto n_idx : neighbours[0]){
        //         if (n_idx != 0){
        //                 colors.row(n_idx) << 0, 1, 0; // track neighbours in green
        //         }
        // }

	t += dt;
}
