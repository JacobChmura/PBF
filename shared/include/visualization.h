#ifndef  VISUALIZATION_H
#define  VISUALIZATION_H

#define IMGUI_DEFINE_MATH_OPERATORS

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>

#include <Eigen/Dense>


namespace Visualize {

	igl::opengl::glfw::Viewer & viewer();

	
    void setup(const Eigen::MatrixXd &fluid_state, const Eigen::MatrixXd &velocity);

    void add_object_to_scene(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, Eigen::RowVector3d color);
    
    void update_vertex_positions(Eigen::Ref<const Eigen::MatrixXd> pos, Eigen::Ref<const Eigen::MatrixXd> colors);
}


#endif
