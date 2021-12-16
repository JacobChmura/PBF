#include <vorticity.h>
#include <kernel.h>

void apply_vorticity(const Eigen::Ref<const Eigen::MatrixXd> &x_new, std::vector<std::vector<int>> &neighbours, Eigen::MatrixXd &v, Eigen::MatrixXd omega, Eigen::MatrixXd eta, Eigen::MatrixXd N, Eigen::MatrixXd vorticity_f, double vorticity_epsilon, double kernel_h, double dt){
        Eigen::Vector3d ker_res;
        Eigen::Vector3d tmp_vec1, tmp_vec2; // for cross product ... should be refactored 
        int num_particles = x_new.rows();
        
        omega.setZero();
        eta.setZero();
        N.setZero();
        // Compute curl
        for(int p_i = 0; p_i < num_particles; p_i++){
                for(int p_j: neighbours[p_i]){
                        kernel_spiky(ker_res, x_new.row(p_i), x_new.row(p_j), kernel_h);
                        tmp_vec1 = v.row(p_j) - v.row(p_i);
                        tmp_vec2 = ker_res;
                        omega.row(p_i) += tmp_vec1.cross(tmp_vec2); 
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
                tmp_vec1 = N.row(p_i);
                tmp_vec2 = omega.row(p_i);
                vorticity_f.row(p_i) = vorticity_epsilon * tmp_vec1.cross(tmp_vec2); 
        }

        v += dt * vorticity_f;
}
