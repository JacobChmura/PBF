#ifndef SETUP_H
#define SETUP_H

#include <Eigen/Dense>
#include <iostream>
#include <map>
#include <string>

/*
Main setup function. Called when program is launched.

        Builds bounding box vertices and edges.
        Writes into fluid state global matrix based on simulation scene.

Input: 
        int simulation_scene: predefined the type of scene we are simulating
                Simulation Scenes:
                        0 = dam fall: cube drops into nothing
                        1 = dam break: cube drops into surface of water
                        2 = double dam fall: two cubes drop into nothing
                        3 = double dam fall: two cubes drop into surface of water
        
        double lower_bound, upper_bound: contain the bounding box extrema world-space coordinates.


Modifies:
        Eigen::MatrixXd &fluid_state (num_particles x 3): pointer the global fluid state matrix. In this function we 
                                        initialize it based on the simulation scene.

        Eigen::MatrixXd &V_box, E_box: containg the geometry of the bounding box.
*/
void setup(int simulation_scene, double lower_bound, double upper_bound, 
                Eigen::MatrixXd &fluid_state, Eigen::MatrixXd &V_box, Eigen::MatrixXi &E_box, std::string simulation_scene_str){

        std::cout << "\n=========== POSITION BASED FLUIDS:  " << simulation_scene_str << " =============\n\tPause Simulation\t\t [SpaceBar]\n\tRestart Dam Fall\t\t [0]\n\tRestart Dam Break\t\t [1]\n\tRestart Double Dam Fall\t\t [2]\n\tRestart Double Dam Break\t [3]\n\n\tToggle User Force Mode\t\t [f]\n\tToggle Vorticity Confinement\t [v]\n\tToggle XPSH Viscocity\t\t [x]\n========================================================\n";

        int num_particles = fluid_state.rows();
        // Boundary box geometry
	V_box <<
        lower_bound, lower_bound, lower_bound,
        upper_bound, lower_bound, lower_bound,
        upper_bound, upper_bound, lower_bound,
        lower_bound, upper_bound, lower_bound,
        lower_bound, lower_bound, upper_bound,
        upper_bound, lower_bound, upper_bound,
        upper_bound, upper_bound, upper_bound,
        lower_bound, upper_bound, upper_bound,

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

        // Write into fluid_state
        switch(simulation_scene){
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
                        int num_particles_block_one = int(num_particles / 8);
                        int num_particles_block_two = int(num_particles / 8);
                        int num_particles_floor = fluid_state.rows() - num_particles_block_one - num_particles_block_two;

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

                        double HI_surface_x = 0.5;
                        double LO_surface_x = -0.5;
                        double range_surface_x = HI_surface_x - LO_surface_x;
                        double HI_surface_y = -0.95;
                        double LO_surface_y = -1;
                        double range_surface_y = HI_surface_y - LO_surface_y;
                        double HI_surface_z = 0.5;
                        double LO_surface_z = -0.5;
                        double range_surface_z = HI_surface_z - LO_surface_y;

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
                        fluid_state.block(num_particles_block_one, 0, block_two.rows(), block_two.cols()) = block_two;
                        fluid_state.block(num_particles_block_one + num_particles_block_two, 0, num_particles_floor, 1) = surface_x;
                        fluid_state.block(num_particles_block_one + num_particles_block_two, 1, num_particles_floor, 1) = surface_y;
                        fluid_state.block(num_particles_block_one + num_particles_block_two, 2, num_particles_floor, 1) = surface_z;
                        break;
                }
                default: {
                        std::cout << "Random Simulation State.\n";
                }
        } 
} 

#endif
