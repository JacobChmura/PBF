#ifndef  VISUALIZATION_H
#define  VISUALIZATION_H

#define IMGUI_DEFINE_MATH_OPERATORS

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>

#include <deque>

#include <Eigen/Dense>


namespace Visualize {

	igl::opengl::glfw::Viewer & viewer();
    igl::opengl::glfw::imgui::ImGuiMenu & viewer_menu();
	
    void setup(const Eigen::MatrixXd &fluid_state);

    void add_object_to_scene(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, Eigen::RowVector3d color);
    void update_vertex_positions(Eigen::Ref<const Eigen::MatrixXd> pos, Eigen::Ref<const Eigen::MatrixXd> colors);

    void get_mouse_down_pos(igl::opengl::glfw::Viewer &viewer, Eigen::Vector3d& mouse_pos);

    void add_energy(float t, float average_density, float max_density);
    bool plot_energy(const char *label, unsigned int type, ImVec2 T_bounds, ImVec2 V_bounds, ImU32 plot_col);
}


#endif
