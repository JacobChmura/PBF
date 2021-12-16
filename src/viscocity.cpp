#include <viscocity.h>
#include <kernel.h>

void apply_viscocity(const Eigen::Ref<const Eigen::MatrixXd> &x_new, std::vector<std::vector<int>> &neighbours, Eigen::MatrixXd &v, Eigen::MatrixXd &v_new, double viscocity_c, double kernel_h){
        int num_particles = x_new.rows();

        v_new = v;
        for(int p_i = 0; p_i < num_particles; p_i++){
                for(int p_j : neighbours[p_i]){
                        v_new.row(p_i) += viscocity_c * kernel_poly6(x_new.row(p_i), x_new.row(p_j), kernel_h) * (v.row(p_i) - v.row(p_j));
                } 
        }

        v = v_new;

}
