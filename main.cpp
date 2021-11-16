#include <iostream>

#include <Particle.h>


// Particle Parameters
double PARTICLE_MASS = 1.0;
double RHO = 6000.0;

// External Force Parameters
double GRAVITY_F = 9.8;
double USER_F = 20.0;

// Jacobi Parameters
int JACOBI_ITERATIONS = 4;

// Constraint Fore Mixing Relaxation
double CFM_EPSILON = 600.0;

// Kernel Parameters
double KERNEL_h = 0.1;

// Tensile Stability Parameters
double TENSILE_k = 0.1;
double TENSILE_delta_q = 0.2 * KERNEL_h;
int TENSILE_n = 4;

// Vorticity Parameters
double VORTICITY_c = 0.01;


int main(int argc, char **argv) {

    std::cout<<"Start PBF \n";
    return 0;
}