#include <vorticity.h>

void apply_vorticity(std::list<Particle> &fluid, double kernel_h){

	// Compute curl
	for (Particle &p_i : fluid){
		p_i.omega.setZero();

		for(Particle &p_j : fluid){
			p_i.omega += (p_j.v - p_i.v) * kernel_spiky(p_i, p_j, kernel_h);
		}
	}

	// Compute differential operator norm
	for (Particle &p_i : fluid){
		p_i.eta.setZero();

		for(Particle &p_j : fluid){
			p_i.eta += p_j.omega.norm() * kernel_spiky(p_i, p_j, kernel_h);
		}
		p_i.N = p_i.eta / p_i.eta.norm();

		// Get vorticity corrective force and update velocity
		p_i.vorticity_f = p_i.N.cross(p_i.omega);
		p_i.v += dt * p_i.vorticity_f;
	}

}