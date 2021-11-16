#include <visocity.h>

void apply_viscocity(std::list<Particle> &fluid, double kernel_h, double viscocity_c){

	for(Particle &p_i : fluid){
		p_i.v_new = p_i.v;
		for(Particle &p_j : p_i.neighbors){
			p_i.v_new += viscocity_c * kernel_poly6(p_i, p_j, kernel_h) * (p_j.v - p_i.v);
		}
	}

	for (Particle &p_i : fluid){
		p_i.v = p_i.v_new;
	}
}