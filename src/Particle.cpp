#include <Particle.h>

Particle::Particle(Eigen::Vector3d x_init, int global_idx){
	this->x = x_init;
        this->global_idx = global_idx;
	this->v.setZero();
	this->f.setZero();
}	

