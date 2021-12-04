#include <iostream>
#include <thread>
#include <visualization.h>
#include <Fluid.h>
#include "setup.h"

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
int num_particles = 500;
bool simulating = true; 
int SIMULATION_SCENE = 0;

// Bounding box vertices and edges 
Eigen::MatrixXd V_box(8,3);
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

bool key_down_callback(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifiers){
        if (key == ' '){ // pause / resume
                if (simulating){
                        std::cout << "Pausing Simulation.\n";
                        simulating = false; 
                }
                else{
                        std::cout << "Resuming Simulation.\n";
                        simulating = true;
                        std::thread simulation_thread(simulate);
                        simulation_thread.detach();
                }
        }
        else if (key == '0') { // restart a dam fall
        }
        else if (key == '1') { // restart a double dam fall
        }
        else if (key == '2') { // restart a dam break
        }
        else if (key == '3') { // restart a double dam break
        }
        else if (key == 'v') { // toggle vorticity confinement
        }
        else if (key == 'x') { // toggle XSPH viscocity
        }
        return false;
}

int main(int argc, char **argv) {
        for (int i = 0; i < num_particles; i++){
                colors(i, 2) = 1;
        }
	std::cout<<"Start PBF \n";

	// --- Initialize setup ----
        setup(SIMULATION_SCENE, LOWER_BOUND, UPPER_BOUND, fluid_state, V_box, E_box);
	fluid.init_state(fluid_state); 
	// --------------------
        
        // Launch new thread for simulation
	std::thread simulation_thread(simulate);
        simulation_thread.detach();

	// Visualization Schema
        Visualize::setup(fluid_state, velocity);
	Visualize::add_object_to_scene(V_box, E_box, Eigen::RowVector3d(1, 0, 0));

	Visualize::viewer().callback_post_draw = &draw;
        Visualize::viewer().callback_key_down = &key_down_callback;
	Visualize::viewer().launch_init(true,false,"Position Based Fluids",0,0);
	Visualize::viewer().launch_rendering(true);
	simulating = false;
	Visualize::viewer().launch_shut();
}
