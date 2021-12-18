#include <vorticity.h>
#include <kernel.h>

void apply_vorticity(const Eigen::Ref<const Eigen::MatrixXd> &x_new, std::vector<std::vector<int>> &neighbours, 
        Eigen::MatrixXd &v, Eigen::MatrixXd omega, Eigen::MatrixXd eta, Eigen::MatrixXd N, 
        Eigen::MatrixXd vorticity_f, double vorticity_epsilon, double kernel_h, double dt){

        Eigen::Vector3d ker_res, v_ji, N_i, omega_i; // for cross product

        int num_particles = x_new.rows();
        omega.setZero();
        eta.setZero();
        N.setZero();

        // Compute curl
        for(int p_i = 0; p_i < num_particles; p_i++){
                for(int p_j: neighbours[p_i]){
                        kernel_spiky(ker_res, x_new.row(p_i), x_new.row(p_j), kernel_h);
                        v_ji = v.row(p_j) - v.row(p_i);
                        omega.row(p_i) += (v_ji).cross(ker_res); 
                }
        }

        // Compute differential operator norm
        for(int p_i = 0; p_i < num_particles; p_i++){
                for(int p_j: neighbours[p_i]){
                        kernel_spiky(ker_res, x_new.row(p_i), x_new.row(p_j), kernel_h);
                        eta.row(p_i) += omega.row(p_j).norm() * ker_res;
                }
                if (eta.row(p_i).norm() > 0){
                        N.row(p_i) = eta.row(p_i) / eta.row(p_i).norm();
                }
                N_i = N.row(p_i);
                omega_i = omega.row(p_i);
                vorticity_f.row(p_i) = vorticity_epsilon * (N_i).cross(omega_i); 
        }

        // Update with force due to vorticity
        v += dt * vorticity_f;
}
