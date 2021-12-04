#include <Particle.h>

Particle::Particle(Eigen::Vector3d x_init, double mass, int global_idx){
	this->x = x_init;
	this->m = mass;
        this->global_idx = global_idx;
	this->v.setZero();
	this->f.setZero();
}	

