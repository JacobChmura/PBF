#include <viscocity.h>

void apply_viscocity(std::vector<Particle> &fluid, double kernel_h, double viscocity_c){

	for(Particle &p_i : fluid){
		p_i.v_new = p_i.v;
		for(int particle_idx : p_i.neighbours){
                        Particle p_j = fluid[particle_idx];
			p_i.v_new += viscocity_c * kernel_poly6(p_i, p_j, kernel_h) * (p_j.v - p_i.v);
		}
	}

	for (Particle &p_i : fluid){
		p_i.v = p_i.v_new;
	}
}
