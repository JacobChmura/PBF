#include <visualization.h> 
#include <iostream>

//libigl viewer
namespace Visualize {

    igl::opengl::glfw::Viewer g_viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    
    Eigen::VectorXd const *g_fluid_state;

    Eigen::Vector3d g_mouse_win;

}

igl::opengl::glfw::Viewer & Visualize::viewer() { return g_viewer; }


// setup
void Visualize::setup(const Eigen::MatrixXd &fluid_state){
	g_viewer.core().background_color.setConstant(1.0);
	g_viewer.core().is_animating = true;
	g_viewer.data().point_size = 5.;
}	

// add_object to scene
void Visualize::add_object_to_scene(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, Eigen::RowVector3d color){

	// Plot the edges of the bounding box
	for (unsigned i=0;i<E.rows(); ++i)
	 g_viewer.data().add_edges
	 (
		V.row(E(i,0)),
		V.row(E(i,1)),
		color
	 );
}

// update vertex positions
void Visualize::update_vertex_positions(Eigen::Ref<const Eigen::MatrixXd> pos, Eigen::Ref<const Eigen::MatrixXd> colors){
	g_viewer.data().set_points(pos,colors);

}


// mouse events 
void Visualize::get_mouse_down_pos(igl::opengl::glfw::Viewer &viewer, Eigen::Vector3d& mouse_pos){
        g_mouse_win = Eigen::Vector3d(g_viewer.current_mouse_x, viewer.core().viewport(3) - g_viewer.current_mouse_y, 0.);
        igl::unproject(g_mouse_win, g_viewer.core().view, g_viewer.core().proj, g_viewer.core().viewport, mouse_pos);
}
