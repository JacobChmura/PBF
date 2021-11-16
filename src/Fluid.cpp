#include <Fluid.h>

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
	this->dt = dt;
}	

Fluid::init_state(){}

Fluid::step(){
	/*
	for all particles i :
		apply externel forces (gravity and user input)
			v_i = v_i + dt * f_ext(i)
		predict positions 
			x_i_new = x_i + dt * v_i

	for all particles i :
		find neighooring particles using x_i_new

	while (iter < num iters){
		for all particles i:
			calculate lambda_i
				lambda_i = -C_i(p_1, ..., p_n) / sum_k|grad_pk C_i|^2 + epsilon

		for all particles i:
			calculate dP :
				calculate s_corr
					s_corr = -k * (poly6(p_i, p_j, h) / poly6(dq, h))^n

				dp = 1 / rho * (sum_j [lambda_i + lambda_j + s_corr)* spiky_kernel(p_i, p_j, h)])

		collision detection and response
			(boundary box)
			(solid object)

		for all particles i:
			update position
				x_i_new += dp

	for all particles i:
		update velocity
			v_i = 1 / dt * (x_new - x)

		apply vorticity confinement
			compute omega_i = sum_j[(v_j - v_i).cross() spiky_kernel(p_i, p_j, h)
			compute eta = nabla|omega_i|
			compute N = eta / eta.norm()

			f_i_vort = eps * (N.cross(omega_i))

		apply viscocity
			v_i_new = v_i + c * sum_j[(v_j - v_i) * poly6(p_i, p_j, h)]
		update position
			x_i = x_i_new
	}
	*/

}