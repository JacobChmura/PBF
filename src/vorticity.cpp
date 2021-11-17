#include <vorticity.h>

void apply_vorticity(std::list<Particle> &fluid, double kernel_h, double vorticity_epsilon, double dt){
	Eigen::Vector3d ker_res;

	// Compute curl
	for (Particle &p_i : fluid){
		p_i.omega.setZero();

		for(Particle &p_j : fluid){
			kernel_spiky(ker_res, p_i, p_j, kernel_h);
			p_i.omega += (p_j.v - p_i.v).cross(ker_res);
		}
	}

	// Compute differential operator norm
	for (Particle &p_i : fluid){
		p_i.eta.setZero();

		for(Particle &p_j : fluid){
			kernel_spiky(ker_res, p_i, p_j, kernel_h);
			p_i.eta += p_j.omega.norm() * ker_res;
		}
		p_i.N = p_i.eta / p_i.eta.norm();

		// Get vorticity corrective force and update velocity
		p_i.vorticity_f = vorticity_epsilon * p_i.N.cross(p_i.omega);
		p_i.v += dt * p_i.vorticity_f;
	}

}