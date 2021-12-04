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
int num_particles = 500;
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
        elif (key == '0') { // restart a dam fall
                pass; 
        }
        elif (key == '1') { // restart a double dam fall
                pass; 
        }
        elif (key == '2') { // restart a dam break
                pass; 
        }
        elif (key == '3') { // restart a double dam break
                pass; 
        }
        elif (key == 'v') { // toggle vorticity confinement
                pass; 
        }
        elif (key == 'x') { // toggle XSPH viscocity
                pass; 
        }
        return false;
}

void build_scene(Eigen::MatrixXd &fluid_state, int simulation_type){ 
       
        switch(simulation_type){
                case 0: { // dam fall
                        double HI = 0.25;
                        double LO = -0.25;
                        double range = HI - LO;       
                        
                        fluid_state = (fluid_state + Eigen::MatrixXd::Constant(fluid_state.rows(), fluid_state.cols(), 1.))*range/2.;
                        fluid_state = (fluid_state + Eigen::MatrixXd::Constant(fluid_state.rows(), fluid_state.cols(), LO));

                        break;
                }
                case 1: {// double dam fall
                        int num_particles_block_one = int(num_particles / 2);
                        int num_particles_block_two = fluid_state.rows() - num_particles_block_one;

                        Eigen::MatrixXd block_one = Eigen::MatrixXd::Random(num_particles_block_one, fluid_state.cols()); 
                        Eigen::MatrixXd block_two = Eigen::MatrixXd::Random(num_particles_block_two, fluid_state.cols()); 

                        double HI_block_one = 0.125;
                        double LO_block_one = -0.125;
                        double range_block_one = HI_block_one - LO_block_one;       
                        
                        double HI_block_two = 0.4;
                        double LO_block_two = 0.15;
                        double range_block_two = HI_block_one - LO_block_one;       

                        block_one = (block_one + Eigen::MatrixXd::Constant(block_one.rows(), block_one.cols(), 1.))*range_block_one/2.;
                        block_one = (block_one + Eigen::MatrixXd::Constant(block_one.rows(), block_one.cols(), LO_block_one));

                        block_two = (block_two + Eigen::MatrixXd::Constant(block_two.rows(), block_two.cols(), 1.))*range_block_two/2.;
                        block_two = (block_two + Eigen::MatrixXd::Constant(block_two.rows(), block_two.cols(), LO_block_two));

                        fluid_state.block(0, 0, block_one.rows(), block_one.cols()) = block_one;
                        fluid_state.block(num_particles_block_one, 0, block_two.rows(), block_two.cols()) = block_two;
                        break;
                }
                case 2: {// dam break

                        int num_particles_block_one = int(num_particles / 4);
                        int num_particles_floor = fluid_state.rows() - num_particles_block_one;
                       
                        double HI_block_one = 0.125;
                        double LO_block_one = -0.125;
                        double range_block_one = HI_block_one - LO_block_one;       

                        double HI_surface_x = 0.5;
                        double LO_surface_x = -0.5;
                        double range_surface_x = HI_surface_x - LO_surface_x;
                        double HI_surface_y = -0.95;
                        double LO_surface_y = -1;
                        double range_surface_y = HI_surface_y - LO_surface_y;
                        double HI_surface_z = 0.5;
                        double LO_surface_z = -0.5;
                        double range_surface_z = HI_surface_z - LO_surface_y;


                        Eigen::MatrixXd block_one = Eigen::MatrixXd::Random(num_particles_block_one, fluid_state.cols()); 
                        block_one = (block_one + Eigen::MatrixXd::Constant(block_one.rows(), block_one.cols(), 1.))*range_block_one/2.;
                        block_one = (block_one + Eigen::MatrixXd::Constant(block_one.rows(), block_one.cols(), LO_block_one));


                        Eigen::VectorXd surface_x = Eigen::VectorXd::Random(num_particles_floor);
                        Eigen::VectorXd surface_y = Eigen::VectorXd::Random(num_particles_floor);
                        Eigen::VectorXd surface_z = Eigen::VectorXd::Random(num_particles_floor);

                        surface_x = (surface_x + Eigen::VectorXd::Constant(num_particles_floor, 1.))*range_surface_x/2;
                        surface_x = (surface_x + Eigen::VectorXd::Constant(num_particles_floor, LO_surface_x));
                        surface_y = (surface_y + Eigen::VectorXd::Constant(num_particles_floor, 1.))*range_surface_y/2;
                        surface_y = (surface_y + Eigen::VectorXd::Constant(num_particles_floor, LO_surface_y));
                        surface_z = (surface_z + Eigen::VectorXd::Constant(num_particles_floor, 1.))*range_surface_z/2;
                        surface_z = (surface_z + Eigen::VectorXd::Constant(num_particles_floor, LO_surface_z));


                        fluid_state.block(0, 0, block_one.rows(), block_one.cols()) = block_one;
                        fluid_state.block(num_particles_block_one, 0, num_particles_floor, 1) = surface_x;
                        fluid_state.block(num_particles_block_one, 1, num_particles_floor, 1) = surface_y;
                        fluid_state.block(num_particles_block_one, 2, num_particles_floor, 1) = surface_z;
                        break;
                }
                case 3: {// double dam break
                        std::cout << "Not efficient enough for this right now.\n";
                        break;
                }
                default: {
                        std::cout << "Random Simulation State.\n";
                }
        } 
}

/* 
 
   scene types
        dam fall: cube drops into nothing
        dam break: cube drops into surface of water
        double dam fall: two cubes drop into nothing
        double dam fall: two cubes drop into surface of water
*/
int main(int argc, char **argv) {
        for (int i = 0; i < num_particles; i++){
                colors(i, 2) = 1;
        }
	std::cout<<"Start PBF \n";
        
        build_scene(fluid_state, 0);

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
        Visualize::viewer().callback_key_down = &key_down_callback;
	Visualize::viewer().launch_init(true,false,"Position Based Fluids",0,0);
	Visualize::viewer().launch_rendering(true);
	simulating = false;
	Visualize::viewer().launch_shut();
}
