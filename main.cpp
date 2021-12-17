#include <iostream>
#include <thread>
#include <visualization.h>
#include <Fluid.h>
#include "setup.h"
#include <unistd.h>
#include <map>
#include <string>

// Particle Parameters
double PARTICLE_MASS = 1.0;
double RHO = 6000.0;
int num_particles = 2000;

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
double VISCOCITY_c = 0.0001;

// Vorticity Parameters
double VORTICITY_EPSILON = 0.0001;

// Simulation Parameters
double dt = 0.001;
bool simulating = true; 
int SIMULATION_SCENE = 0;
std::map<int, std::string> SIMULATION_MODE = {{0, "Dam Fall"}, {1, "Dam Break"}, {2, "Double Dam Fall"}, {3, "Double Dam Break"}};

// Bounding Box Extrema
double LOWER_BOUND = -1;
double UPPER_BOUND = 1;

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


// test mouse
Eigen::Vector3d mouse_pos;
bool add_user_force;
bool user_force_mode;
bool use_viscocity;
bool use_vorticity;

void simulate(){
        int flag;
	while(simulating){
		//std::cout << "step.\n";
		fluid.step(fluid_state, colors, mouse_pos, add_user_force, use_viscocity, use_vorticity);
                //flag = getchar(); 
        }
}

bool draw(igl::opengl::glfw::Viewer &viewer) {
    //update vertex positions using simulation
    Visualize::update_vertex_positions(fluid_state, colors);
    return false;
}

bool mouse_down(igl::opengl::glfw::Viewer &viewer, int x, int y){
        Visualize::get_mouse_down_pos(viewer, mouse_pos);
        if (user_force_mode) add_user_force  = true;
        return false;
}

bool mouse_up(igl::opengl::glfw::Viewer &viewer, int x, int y){
        add_user_force = false;
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
        else if (key == '0' || key == '1' || key == '2' || key == '3') { // restart a simulation
                simulating = false;

                sleep(1); // not sure about this
                std::cout << "Restarting simulation " << SIMULATION_MODE[key-'0'] << std::endl;
                fluid_state = Eigen::MatrixXd::Random(num_particles, 3);
                setup(key - '0', LOWER_BOUND, UPPER_BOUND, fluid_state, V_box, E_box, SIMULATION_MODE[key-'0']);
                fluid.init_state(fluid_state); 

                simulating = true;
                std::thread simulation_thread(simulate);
                simulation_thread.detach();
        }
        else if (key == '4'){ // toggle user force mode
                user_force_mode = !user_force_mode;
                if(user_force_mode){
                        std::cout << "User Force Mode Enabled.\n";
                }
                else{
                        std::cout << "User Force Mode Disabled.\n";
                }
        }
        else if (key == '5') { // toggle vorticity confinement
                use_vorticity = !use_vorticity;
                if(use_vorticity){
                        std::cout << "Vorticity Enabled.\n";
                }
                else{
                        std::cout << "Vorticity Disabled.\n";
                }
        }
        else if (key == '6') { // toggle XSPH viscocity
                use_viscocity = !use_viscocity;
                if(use_viscocity){
                        std::cout << "Viscocity Enabled.\n";
                }
                else{
                        std::cout << "Viscocity Disabled.\n";
                }
        }
        return false;
}

int main(int argc, char **argv) {

        // Parse command line arguments
        if (argc > 1){
                if (argc > 2){
                        std::cerr << "Usage: " << argv[0] << " SIMULATION SCENE" << std::endl;
                        return 1;
                }

                SIMULATION_SCENE = std::stoi(argv[1]);
        }
        for (int i = 0; i < num_particles; i++){
                colors(i, 2) = 1;
        }

	// Initialize setup
        setup(SIMULATION_SCENE, LOWER_BOUND, UPPER_BOUND, fluid_state, V_box, E_box, SIMULATION_MODE[SIMULATION_SCENE]);
	fluid.init_state(fluid_state); 
        
        // Launch new thread for simulation
	std::thread simulation_thread(simulate);
        simulation_thread.detach();

	// Visualization Schema
        Visualize::setup(fluid_state, velocity);
	Visualize::add_object_to_scene(V_box, E_box, Eigen::RowVector3d(1, 0, 0));

	Visualize::viewer().callback_post_draw = &draw;
        Visualize::viewer().callback_key_down = &key_down_callback;

        // user force
        Visualize::viewer().callback_mouse_down = &mouse_down;
        Visualize::viewer().callback_mouse_up = &mouse_up;

	Visualize::viewer().launch_init(true,false,"Position Based Fluids",0,0);
	Visualize::viewer().launch_rendering(true);
	simulating = false;
	Visualize::viewer().launch_shut();
}

