#include <iostream>
#include <unistd.h>
#include <thread>

#include <map>
#include <string>

#include <visualization.h>
#include <Fluid.h>
#include "setup.h"

#define CAPTURE 0 // if you want to save frames


// ---------------- Constants ----------------
double dt = 0.0001;
double particle_mass = 1.0;
double kernel_h = 0.1;

// External forces
double gravity_f = 9.8; 
double user_f = 20.0; 

// Bounding Box
double lower_bound = -1;
double upper_bound = 1;
Eigen::MatrixXd V_box(8,3);
Eigen::MatrixXi E_box(12,2);

int num_particles = 2000;
int jacobi_iterations = 4;
double rho = 10000.0; // rest density
double cfm_epsilon = 60.0; // Constraint Force Mixing Relaxation

// Tensile Stability 
int tensile_n = 4;
double tensile_k = 0.1;
double tensile_delta_q = 0.2 * kernel_h;

// Visocity / Vorticity 
double viscocity_c = 0.0001;
double vorticity_epsilon = 0.0001;


// ----------- Simulation Parameters ------------
bool simulating = true; 
int simulation_scene = 0;
std::map<int, std::string> simulation_modes = {{0, "Dam Fall"}, {1, "Dam Break"}, {2, "Double Dam Fall"}, {3, "Double Dam Break"}, {4, "Floor"}};



// ----------- Global state -------------
Fluid fluid(particle_mass, rho, gravity_f, user_f, jacobi_iterations, 
			cfm_epsilon, kernel_h, tensile_k, tensile_delta_q, tensile_n, 
			viscocity_c, vorticity_epsilon, lower_bound, upper_bound, dt);
Eigen::MatrixXd fluid_state;
Eigen::MatrixXd colors;
Eigen::MatrixXd velocity;

// ----------- Input / Output ------------
bool add_user_force;
bool user_force_flag;
bool use_viscocity = true;
bool use_vorticity = true;
Eigen::Vector3d mouse_pos;

bool write_frame_flag;
int frame_number = 0;
std::string experiment_id = "debug";


// ----------- Simulation Loop -----------
void simulate(){
        while(simulating){
                fluid.step(fluid_state, colors, mouse_pos, add_user_force, use_viscocity, use_vorticity);
                Visualize::add_density(fluid.t, fluid.avg_density, fluid.max_density);
                write_frame_flag = true;               
        }    
}

// ----------- IGL Callbacks -------------
bool draw(igl::opengl::glfw::Viewer &viewer) {
        //update vertex positions using simulation
        Visualize::update_vertex_positions(fluid_state, colors);

        // if we are intending to capture frames, and we have received an update from the fluid, use openGL to read all pixels
        if (CAPTURE && write_frame_flag){
                std::string path = "../data/" + experiment_id + "/frames/" + std::to_string(frame_number++) + ".png";
                Visualize::write_frame(path);
                write_frame_flag = false;
        }
        return false;
}

bool mouse_down(igl::opengl::glfw::Viewer &viewer, int x, int y){
        Visualize::get_mouse_down_pos(viewer, mouse_pos);
        if (user_force_flag) add_user_force  = true;
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
        else if (key == '0' || key == '1' || key == '2' || key == '3' || key == '4') { // restart a simulation
                simulating = false;
                sleep(1);
                std::cout << "Restarting simulation " << simulation_modes[key-'0'] << std::endl;
                setup(num_particles, key - '0', lower_bound, upper_bound, fluid_state, V_box, E_box, colors, simulation_modes[key-'0']);
                fluid.init_state(fluid_state); 

                simulating = true;
                std::thread simulation_thread(simulate);
                simulation_thread.detach();
        }
        else if (key == 'F'){ // toggle user force mode
                user_force_flag = !user_force_flag;
                if(user_force_flag){
                        std::cout << "User Force Mode Enabled.\n";
                }
                else{
                        std::cout << "User Force Mode Disabled.\n";
                }
        }
        else if (key == 'V') { // toggle vorticity confinement
                use_vorticity = !use_vorticity;
                if(use_vorticity){
                        std::cout << "Vorticity Enabled.\n";
                }
                else{
                        std::cout << "Vorticity Disabled.\n";
                }
        }
        else if (key == 'X') { // toggle XSPH viscocity
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
        if (parse_args(argc, argv, simulation_scene, num_particles, experiment_id)) return 1;

	// Initialize setup
        setup(num_particles, simulation_scene, lower_bound, upper_bound, fluid_state, V_box, E_box, colors, simulation_modes[simulation_scene]);
	fluid.init_state(fluid_state); 
        
        // Launch new thread for simulation
	std::thread simulation_thread(simulate);
        simulation_thread.detach();

	// Visualization Schema
        Visualize::setup(fluid_state);
	Visualize::add_object_to_scene(V_box, E_box, Eigen::RowVector3d(1, 0, 0));

	Visualize::viewer().callback_post_draw = &draw;
        Visualize::viewer().callback_key_down = &key_down_callback;

        // user force
        Visualize::viewer().callback_mouse_down = &mouse_down;
        Visualize::viewer().callback_mouse_up = &mouse_up;

        /* 
        Density Chart
        Reference: https://github.com/dilevin/CSC417-a3-finite-elements-3d
        */ 
        Visualize::viewer_menu().callback_draw_custom_window = [&](){
                // Define next window position + size
                ImGui::SetNextWindowPos(ImVec2(180.f * Visualize::viewer_menu().menu_scaling(), 0), ImGuiCond_FirstUseEver);
                ImGui::SetNextWindowSize(ImVec2(100, 10), ImGuiCond_FirstUseEver);
                ImGui::Begin("Density Plot", nullptr, ImGuiWindowFlags_NoSavedSettings);

                ImVec2 min = ImGui::GetWindowContentRegionMin();
                ImVec2 max = ImGui::GetWindowContentRegionMax();
                max.x = ( max.x - min.x ) / 2;
                max.y -= min.y + ImGui::GetTextLineHeightWithSpacing() * 3;

                Visualize::plot_density("Average Density", 1, ImVec2(0,0.1), ImVec2(0,2 * rho), ImGui::GetColorU32(ImGuiCol_PlotLines));
                Visualize::plot_density("Max Density", 2, ImVec2(0,0.1), ImVec2(0,2 * rho), ImGui::GetColorU32(ImGuiCol_HeaderActive));

                ImGui::End();
        };

        Visualize::viewer().launch_init(true,false,"Position Based Fluids",0,0);
        Visualize::viewer().launch_rendering(true);
}

