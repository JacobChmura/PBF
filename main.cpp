#include <iostream>
#include <Particle.h>
#include <Eigen/Dense>
#include <thread>
#include <visualization.h>
#include <Fluid.h>

// Particle Parameters
double PARTICLE_MASS = 1.0;
double RHO = 6000.0;

// External Force Parameters
double GRAVITY_F = 9.8;
double USER_F = 20.0;

// Jacobi Parameters
int JACOBI_ITERATIONS = 1;

// Constraint Force Mixing Relaxation
double CFM_EPSILON = 600.0;

// Kernel Parameters
double KERNEL_h = 0.1;

// Tensile Stability Parameters
double TENSILE_k = 0.1;
double TENSILE_delta_q = 0.2 * KERNEL_h;
int TENSILE_n = 4;

// Visocity Parameters
double VISCOCITY_c = 0.01;

// Vorticity Parameters
double VORTICITY_EPSILON = 0.0005;

// Simulation Parameters
double dt = 0.01;

// Bounding Box Extrema
double LOWER_BOUND = -1;
double UPPER_BOUND = 1;

// dummy test
int num_particles = 1000;
bool simulating = true; 


// Bounding box
Eigen::Vector3d m;
Eigen::Vector3d M;

// Corners of the bounding box
Eigen::MatrixXd V_box(8,3);

// Edges of the bounding box
Eigen::MatrixXi E_box(12,2);

// Global state
Fluid fluid(num_particles, PARTICLE_MASS, RHO, GRAVITY_F, USER_F, JACOBI_ITERATIONS, 
			CFM_EPSILON, KERNEL_h, TENSILE_k, TENSILE_delta_q, TENSILE_n, 
			VISCOCITY_c, VORTICITY_EPSILON, LOWER_BOUND, UPPER_BOUND, dt);

Eigen::MatrixXd fluid_state = Eigen::MatrixXd::Random(num_particles, 3);
Eigen::MatrixXd colors = Eigen::MatrixXd::Zero(num_particles, 3);
Eigen::MatrixXd velocity = Eigen::MatrixXd::Zero(num_particles, 3);
// --------------

void simulate(){
        int flag;
	while(simulating){
		std::cout << "step.\n";
		fluid.step(fluid_state, colors);
                //flag = getchar();
	}
}

bool draw(igl::opengl::glfw::Viewer &viewer) {
    //update vertex positions using simulation
    Visualize::update_vertex_positions(fluid_state, colors);
    return false;
}


int main(int argc, char **argv) {
        for (int i = 0; i < num_particles; i++){
                colors(i, 2) = 1;
        }
	std::cout<<"Start PBF \n";
        
        // TEST
        for (int i = 0; i < num_particles; i++){
                fluid_state(i, 1) = 0;
                if (fluid_state(i, 0) < 0){
                        fluid_state(i, 0) = 0;
                }
        }
	// --- Initialize setup ----
	m << LOWER_BOUND, LOWER_BOUND, LOWER_BOUND;	
	M << UPPER_BOUND, UPPER_BOUND, UPPER_BOUND;

	V_box <<
	m(0), m(1), m(2),
	M(0), m(1), m(2),
	M(0), M(1), m(2),
	m(0), M(1), m(2),
	m(0), m(1), M(2),
	M(0), m(1), M(2),
	M(0), M(1), M(2),
	m(0), M(1), M(2);
	E_box <<
	0, 1,
	1, 2,
	2, 3,
	3, 0,
	4, 5,
	5, 6,
	6, 7,
	7, 4,
	0, 4,
	1, 5,
	2, 6,
	7 ,3;
	// --------------------

	fluid.init_state(fluid_state); // random init
        
        /// for debugging
	std::thread simulation_thread(simulate);
        simulation_thread.detach();

	// Visualization Schema
        Visualize::setup(fluid_state, velocity);
	Visualize::add_object_to_scene(V_box, E_box, Eigen::RowVector3d(1, 0, 0));

	Visualize::viewer().callback_post_draw = &draw;
	Visualize::viewer().launch_init(true,false,"Position Based Fluids",0,0);
	Visualize::viewer().launch_rendering(true);
	simulating = false;
	Visualize::viewer().launch_shut();
}
