#include <Fluid.h>
#include <kernel.h>
#include <vorticity.h>
#include <viscocity.h>

#include <iostream>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

#define DEBUG 0

Fluid::Fluid(int num_particles, double particle_mass, double rho, double gravity_f, double user_f, int jacobi_iterations, 
			double cfm_epsilon, double kernel_h, double tensile_k, double tensile_delta_q, int tensile_n, 
			double viscocity_c, double vorticity_epsilon, double lower_bound, double upper_bound, double dt){

        this->num_particles = num_particles;
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

        this->gravity_f = Eigen::VectorXd::Constant(num_particles, gravity_f);
        this->user_f = Eigen::VectorXd::Constant(num_particles, user_f);

        std::vector<std::vector<int>> neighbours_init(num_particles);
        this->neighbours = neighbours_init;

        this->tensile_stability_denom = pow(kernel_poly6(tensile_delta_q, kernel_h), tensile_n);

}	

void Fluid::init_state(Eigen::MatrixXd &fluid_state){
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
        if (DEBUG) std::cout << "\n------------------------------------------\nApplied Externel forces [" << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " s]\n";

        // Get neighbours using spatial hash grid
        grid.findNeighbours(x_new, neighbours);

        auto t2 = Clock::now();
        if (DEBUG) std::cout << "Found Neighbours [" << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " s]\n";

	//Jacobi Loop 
	for(int i = 0; i < jacobi_iterations; i++){
                // reset jacobi state
                density.setZero();
                lambda.setZero();
                c_grad_norm.setZero();

                //test      
                Eigen::ArrayXXd NEIGHBOURS_MATRIX = Eigen::ArrayXXd::Constant(num_particles, num_particles, 0.);

                for(int ii = 0; ii < neighbours.size(); ii++){
                        for(auto jj: neighbours[ii]){
                                NEIGHBOURS_MATRIX(ii, jj) = 1;
                        }
                }

                //std::cout << "Neighbours MATRIX\n" << NEIGHBOURS_MATRIX << std::endl;
                
                // Eigen::MatrixXd C_GRAD_PAIRWISE_SERIAL = Eigen::MatrixXd::Constant(num_particles, num_particles, 0.);;
                // Eigen::MatrixXd POLY_6_SERIAL = Eigen::MatrixXd::Constant(num_particles, num_particles, 0.);


                // auto serial_start = Clock::now();
                // for(int p_i = 0; p_i < num_particles; p_i++){
                        
                //         // compute densities
                //         for(int p_j : neighbours[p_i]){
                //                 POLY_6_SERIAL(p_i, p_j) = kernel_poly6(x_new.row(p_i), x_new.row(p_j), kernel_h);

                //                 density[p_i] += particle_mass * kernel_poly6(x_new.row(p_i), x_new.row(p_j), kernel_h);

                //                 // if (p_i == 0){
                //                 // // SUM OF NORM IS NOT NORM OF SUM
                //                 //         std::cout << "Particle : " << p_i << ", " << p_j << " density add term " << particle_mass * kernel_poly6(x_new.row(p_i), x_new.row(p_j), kernel_h) << std::endl;
                //                 // }


                //                 // accumulate gradient norm
                //                 c_grad_temp.setZero();

                //                 // DEBUGGING!!!
                //                 if (p_i == p_j){
                //                 //if (true){ 
                //                         for(int p_k : neighbours[p_i]){
                //                                kernel_spiky(ker_res, x_new.row(p_i), x_new.row(p_k), kernel_h);
                //                                c_grad_temp += ker_res;
                //                                if (p_i == 1){
                //                                 // SUM OF NORM IS NOT NORM OF SUM
                //                                       //  std::cout << "Particle : " << p_i << ", " << p_k << " ker res norm: / rho " << ker_res.norm() / rho << std::endl;
                //                                }
                //                         }
                //                         C_GRAD_PAIRWISE_SERIAL(p_i, p_j) = (c_grad_temp / rho).norm();
                //                 }
                //                 else{
                //                         kernel_spiky(ker_res, x_new.row(p_i), x_new.row(p_j), kernel_h);
                //                         c_grad_temp = ker_res;

                //                         //std::cout << "Particle : " << p_i << ", " << p_j << " ker res norm: " << ker_res.norm() << std::endl;
                //                 }

                                
                //                 c_grad_norm[p_i] += (c_grad_temp / rho).norm();
                //         }

                //         // Compute constraint and lambda
                //         c[p_i] = (density[p_i] / rho) - 1;
                //         lambda[p_i] = -c[p_i] / (c_grad_norm[p_i] + cfm_epsilon);

                // }

                // auto serial = Clock::now();
                //std::cout << "Serial -> " << std::chrono::duration_cast<std::chrono::milliseconds>(serial - serial_start).count() << " s]\n";


                // ======================= VECTORIZED OPS ========================

                // --------- Stuff we can PreCompute --------
                //auto vec_time_0 = Clock::now();

                Eigen::VectorXd ONES = Eigen::VectorXd::Constant(num_particles, 1.);
                Eigen::ArrayXXd ONES_BLOCK = Eigen::ArrayXXd::Constant(num_particles, num_particles, 1.);
                Eigen::ArrayXXd ZEROS_BLOCK = Eigen::ArrayXXd::Constant(num_particles, num_particles, 0.);

                double poly6_prefactor = 315.0 / (64.0 * PI * pow(kernel_h, 9));
                double spiky_prefactor = 45.0 / (PI * pow(kernel_h, 6));
                double inverse_rho = 1. / rho;

                Eigen::ArrayXd CFM_EPS_VECTOR = Eigen::ArrayXd::Constant(num_particles, cfm_epsilon);
                Eigen::ArrayXXd kernel_h_matrix = Eigen::ArrayXXd::Constant(num_particles, num_particles, kernel_h);
                Eigen::MatrixXd off_diag_boolean = Eigen::MatrixXd::Constant(num_particles, num_particles, 1.);
                off_diag_boolean.diagonal().setZero();
                Eigen::ArrayXXd OFF_DIAG_BOOLEAN = off_diag_boolean;

                //auto vec_time_1 = Clock::now();
                //std::cout << "Pre-Step Compute Time [" << std::chrono::duration_cast<std::chrono::milliseconds>(vec_time_1 - vec_time_0).count() << " ms]\n";

                // --------- Distance Matrix --------------
                Eigen::ArrayXXd DISTANCE_MATRIX;
                Eigen::MatrixXd XXt = x_new * x_new.transpose();
                Eigen::MatrixXd NORM = XXt.diagonal() * ONES.transpose();
                DISTANCE_MATRIX = (NORM + NORM.transpose() - 2 * XXt);
                DISTANCE_MATRIX = (DISTANCE_MATRIX.max(0.)).sqrt();

                Eigen::ArrayXXd DISTANCE_WITHIN_KERNEL_MASK;
                DISTANCE_WITHIN_KERNEL_MASK = (kernel_h - DISTANCE_MATRIX < 0).select(ZEROS_BLOCK, ONES_BLOCK);
                //auto vec_time_2 = Clock::now();
                //std::cout << "Distance Matrix Compute Time [" << std::chrono::duration_cast<std::chrono::milliseconds>(vec_time_2 - vec_time_1).count() << " ms]\n";

                // ---------------Poly 6 kernel ------------
                Eigen::ArrayXXd POLY_6 =  NEIGHBOURS_MATRIX * DISTANCE_WITHIN_KERNEL_MASK * (poly6_prefactor * (kernel_h_matrix.square() - DISTANCE_MATRIX.square()).cube());

                //auto vec_time_3 = Clock::now();
               //std::cout << "Poly6 Compute Time [" << std::chrono::duration_cast<std::chrono::milliseconds>(vec_time_3 - vec_time_2).count() << " ms]\n";


                // ------------Density / Constraint ---------
                Eigen::ArrayXd DENSITY_RESULT = (particle_mass * POLY_6).rowwise().sum();
                Eigen::ArrayXd CONSTRAINT_RESULT = inverse_rho * DENSITY_RESULT - 1.;

                //auto vec_time_4 = Clock::now();
                //std::cout << "Constraint Compute Time [" << std::chrono::duration_cast<std::chrono::milliseconds>(vec_time_4 - vec_time_3).count() << " ms]\n";
                

               // std::cout << "\n\n";
                
               // std::cout << "\t\t Poly 6 Diff: " << (POLY_6_SERIAL - POLY_6.matrix()).norm() << std::endl;
                //std::cout << "Poly 6 Serial: \n" << POLY_6_SERIAL << std::endl;
              //  std::cout << "Poly 6 Vector: \n" << POLY_6 << std::endl;



                // if (POLY_6.hasNaN()){
                //         std::cout << "Distance Matrix: " << DISTANCE_MATRIX << "\n\n";

                //         std::cout << "DISTANCE_MATRIX pre square root:\n" << (NORM + NORM.transpose() - 2 * XXt) << std::endl;
                //         exit(0);
                // }
                //std::cout << "\t\tDensity Diff: " << (density - DENSITY_RESULT.matrix()).norm() << std::endl;
               // std::cout << "Density Serial: \n" << density.transpose() << std::endl;
               // std::cout << "Density Vector: \n" << DENSITY_RESULT.transpose() << std::endl;

              //  std::cout << "\t\tConstraint Diff: " << (CONSTRAINT_RESULT.matrix() - c).norm() << std::endl;
               // std::cout << "Constraints Serial: \n" << c.transpose() << std::endl;
                //std::cout << "Constraint Vector: \n" << CONSTRAINT_RESULT.transpose() << std::endl;
              //  std::cout << "\n\n";
                //exit(0);
                // ----------- Off diag spiky -------------

                Eigen::ArrayXXd SPIKY_MATRIX_SHARED = DISTANCE_WITHIN_KERNEL_MASK * (spiky_prefactor * (kernel_h_matrix - DISTANCE_MATRIX).square());



                Eigen::ArrayXXd SPIKY_MATRIX_OFF_DIAG = inverse_rho * OFF_DIAG_BOOLEAN * SPIKY_MATRIX_SHARED;
                //SPIKY_MATRIX_OFF_DIAG *= (inverse_rho * OFF_DIAG_BOOLEAN);
  
               // auto vec_time_5 = Clock::now();
               // std::cout << "Off diag Compute Time [" << std::chrono::duration_cast<std::chrono::milliseconds>(vec_time_5 - vec_time_4).count() << " ms]\n";



                // ----------- On Diag Spiky ---------------
                Eigen::ArrayXXd NON_SELF_NEIGHBOURS = NEIGHBOURS_MATRIX * OFF_DIAG_BOOLEAN;
                
                Eigen::ArrayXXd SPIKY_MATRIX_ON_DIAG = inverse_rho * NON_SELF_NEIGHBOURS * SPIKY_MATRIX_SHARED;
                SPIKY_MATRIX_ON_DIAG = (SPIKY_MATRIX_ON_DIAG.rowwise().sum());

                // Full spiky matrix
                Eigen::MatrixXd spiky_matrix_on_diag = SPIKY_MATRIX_ON_DIAG.matrix().asDiagonal();
                Eigen::ArrayXXd SPIKY_MATRIX = (SPIKY_MATRIX_OFF_DIAG + spiky_matrix_on_diag.array());

                //auto vec_time_6 = Clock::now();
                //std::cout << "On diag Compute Time [" << std::chrono::duration_cast<std::chrono::milliseconds>(vec_time_6 - vec_time_5).count() << " ms]\n";



                // ---------- C grad norm and Lambda -------------
                Eigen::ArrayXd C_GRAD_RESULT = SPIKY_MATRIX.rowwise().sum();
                Eigen::ArrayXd DENOM = (C_GRAD_RESULT + CFM_EPS_VECTOR).inverse();
                Eigen::ArrayXd LAMBDA_VECTORIZED = - CONSTRAINT_RESULT * DENOM;
  
                //auto vec_time_7 = Clock::now();
               // std::cout << "Lambda Compute Time [" << std::chrono::duration_cast<std::chrono::milliseconds>(vec_time_7 - vec_time_6).count() << " ms]\n";

                //auto vector = Clock::now();
                //std::cout << "Vector -> " << std::chrono::duration_cast<std::chrono::milliseconds>(vector - vec_time_1).count() << " s]\n";
                
                lambda = LAMBDA_VECTORIZED;
                c = C_GRAD_RESULT;
                
                auto t3 = Clock::now();
                if (DEBUG) std::cout << "Computed Constraints [" << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() << " s]\n";



                //------- DP ---------------------
  //              std::cout << "\n\n";

  //               auto serial_dp_start = Clock::now();

  //               Eigen::ArrayXXd SERIAL_S_CORR = Eigen::ArrayXXd::Constant(num_particles, num_particles, 0.);

		// //Compute dP
  //               Eigen::ArrayXXd test_weights = Eigen::ArrayXXd::Constant(num_particles, num_particles, 0.);
  //               dP.setZero();
  //               for(int p_i = 0; p_i < num_particles; p_i++){
  //                       for(int p_j : neighbours[p_i]){
  //                               kernel_spiky(ker_res, x_new.row(p_i), x_new.row(p_j), kernel_h);

		//                 s_corr = - tensile_k * pow(kernel_poly6(x_new.row(p_i), x_new.row(p_j), kernel_h), tensile_n) / tensile_stability_denom;
  //                               SERIAL_S_CORR(p_i, p_j) = s_corr;
                                
  //                               dP.row(p_i) += (lambda[p_i] + lambda[p_j] + s_corr) * ker_res; 
  //                               test_weights(p_i, p_j) = lambda[p_i] + lambda[p_j] + s_corr;
  //                       }
  //                       dP.row(p_i) /= rho;
  //               }

                // auto serial_dp_end = Clock::now();
                // std::cout << "Serial DP Compute Time:  [" << std::chrono::duration_cast<std::chrono::milliseconds>(serial_dp_end - serial_dp_start).count() << " ms]\n";
                
               // auto vector_dp_start = Clock::now();

                Eigen::ArrayXXd S_CORR = - tensile_k * POLY_6.pow(tensile_n) / tensile_stability_denom;

                Eigen::ArrayXXd OUTER_SUM_LAMBDAS = LAMBDA_VECTORIZED.rowwise().replicate(LAMBDA_VECTORIZED.size()) + LAMBDA_VECTORIZED.transpose().colwise().replicate(LAMBDA_VECTORIZED.size());
                Eigen::ArrayXXd WEIGHT_FACTOR;
                WEIGHT_FACTOR = (DISTANCE_MATRIX < 0.00001).select(ZEROS_BLOCK, (OUTER_SUM_LAMBDAS + S_CORR) * (- SPIKY_MATRIX_SHARED / DISTANCE_MATRIX));
                WEIGHT_FACTOR *= (inverse_rho * NEIGHBOURS_MATRIX);

                Eigen::ArrayXd X_NEW_X_COORD = x_new.col(0);
                Eigen::ArrayXd X_NEW_Y_COORD = x_new.col(1);
                Eigen::ArrayXd X_NEW_Z_COORD = x_new.col(2);

                
                Eigen::ArrayXXd OUTER_SUB_X_NEW_X_COORD = X_NEW_X_COORD.rowwise().replicate(X_NEW_X_COORD.size()) - X_NEW_X_COORD.transpose().colwise().replicate(X_NEW_X_COORD.size());
                Eigen::ArrayXXd OUTER_SUB_X_NEW_Y_COORD = X_NEW_Y_COORD.rowwise().replicate(X_NEW_Y_COORD.size()) - X_NEW_Y_COORD.transpose().colwise().replicate(X_NEW_Y_COORD.size());
                Eigen::ArrayXXd OUTER_SUB_X_NEW_Z_COORD = X_NEW_Z_COORD.rowwise().replicate(X_NEW_Z_COORD.size()) - X_NEW_Z_COORD.transpose().colwise().replicate(X_NEW_Z_COORD.size());
                
                Eigen::ArrayXd DP_X = (WEIGHT_FACTOR * OUTER_SUB_X_NEW_X_COORD).rowwise().sum();
                Eigen::ArrayXd DP_Y = (WEIGHT_FACTOR * OUTER_SUB_X_NEW_Y_COORD).rowwise().sum();
                Eigen::ArrayXd DP_Z = (WEIGHT_FACTOR * OUTER_SUB_X_NEW_Z_COORD).rowwise().sum();

                Eigen::ArrayXXd DP_VECTORIZED = Eigen::ArrayXXd::Constant(num_particles, 3, 0.);
                DP_VECTORIZED << DP_X, DP_Y, DP_Z;
                //std::cout << "\t\tdP Diff: " << (DP_VECTORIZED.matrix() - dP).norm() << std::endl;
                

                //auto vector_dp_end = Clock::now();
                //std::cout << "Vector DP Compute Time: [" << std::chrono::duration_cast<std::chrono::milliseconds>(vector_dp_end - vector_dp_start).count() << " ms]\n";
                

                dP = DP_VECTORIZED;
                // ---------------------

                auto t4 = Clock::now();
                if (DEBUG) std::cout << "Computed Position Correction [" << std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count() << " s]\n";


                // --------------------- COLLISION DETECTION ----------------
                x_new += 0.005 * dP;
		// Collision Detection with boundary box and solid
                for(int p_i = 0; p_i < num_particles; p_i++){

			//Collision Detect correct p.x_new (naive)
                        for(int axis = 0; axis < 3; axis++){
                                if (x_new.row(p_i)[axis] < lower_bound){ 
                                        x_new.row(p_i)[axis] = lower_bound;
                                        if (v(p_i, axis) < 0) v(p_i, axis) *= -1;
                                }
                                if (x_new.row(p_i)[axis] > upper_bound){ 
                                        x_new.row(p_i)[axis] = upper_bound;
                                        if (v(p_i, axis) > 0) v(p_i, axis) *= -1;
                                }
                        }
                }



                auto t5 = Clock::now();
               // std::cout << "Seralized Collision Detection Compute Time: [" << std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4).count() << " ms]\n";
                if (DEBUG) std::cout << "Collision Detection [" << std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4).count() << " s]\n";

              //  std::cout << "\n\n";
	}

	//Update Velocity
        v = (x_new - fluid_state) / dt;
        
        // Vorticity and Viscocity
        auto t6 = Clock::now();

        if (use_vorticity) apply_vorticity(x_new, neighbours, v, omega, eta, N, vorticity_f, vorticity_epsilon, kernel_h, dt);
        
        auto t7 = Clock::now();
        if (DEBUG) std::cout << "Vorticity [" << std::chrono::duration_cast<std::chrono::milliseconds>(t7 - t6).count() << " s]\n";

        if (use_viscocity) apply_viscocity(x_new, neighbours, v, v_new, viscocity_c, kernel_h);
        auto t8 = Clock::now();
        if (DEBUG) std::cout << "Viscocity [" << std::chrono::duration_cast<std::chrono::milliseconds>(t8 - t7).count() << " s]\n";

	// Update Position and spatial hash grid
        fluid_state = x_new;
        grid.update(fluid_state);

        auto t9 = Clock::now();
        if (DEBUG) std::cout << "Simulation Step Total Time [" << std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t0).count() << " s]\n----------------------------------------\n";

        // Debugging using colors for now.
        colors.row(0) << 1, 0, 0; // track particle 0 in red
        for (int i = 1; i < num_particles; i++){
                colors.row(i) << 0, 0, 1; 
        }
        for (auto n_idx : neighbours[0]){
                if (n_idx != 0){
                        colors.row(n_idx) << 0, 1, 0; // track neighbours in green
                }
        }

	t += dt;
}
