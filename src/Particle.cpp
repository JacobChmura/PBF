#include <Particle.h>

Particle::Particle(Eigen::Vector3d X_init, double mass){
	this->m = mass;
	this->x = X_init;
	this->X = X_init;

	this->v.setZero();
	this->f.setZero();
}	

