#include <visualization.h> 
#include <iostream>

//libigl viewer
namespace Visualize {

	igl::opengl::glfw::Viewer g_viewer;
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	
	Eigen::VectorXd const *g_fluid_state;
	Eigen::Vector3d g_mouse_win;

	std::deque<std::array<float, 3> > g_density; // time, average fluid density, max fluid density
}

igl::opengl::glfw::Viewer & Visualize::viewer() { return g_viewer; }
igl::opengl::glfw::imgui::ImGuiMenu & Visualize::viewer_menu() { return menu; }


/* 
Reference: https://github.com/dilevin/CSC417-a3-finite-elements-3d
*/ 
void Visualize::setup(const Eigen::MatrixXd &fluid_state){

	Visualize::g_viewer.plugins.push_back(&menu);
	menu.callback_draw_viewer_menu = [&]()
		{
			ImGuiStyle& style = ImGui::GetStyle();
			style.WindowRounding = 5.3f;
			style.FrameRounding = 2.3f;
			style.ScrollbarRounding = 0;

			style.Colors[ImGuiCol_Text]                  = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);
			style.Colors[ImGuiCol_TextDisabled]          = ImVec4(0.70f, 0.70f, 0.70f, 1.00f);
			style.Colors[ImGuiCol_WindowBg]              = ImVec4(0.8f, 0.8f, 0.8f, 1.00f);
			style.Colors[ImGuiCol_PopupBg]               = ImVec4(0.05f, 0.05f, 0.10f, 0.85f);
			style.Colors[ImGuiCol_Border]                = ImVec4(0.70f, 0.70f, 0.70f, 0.65f);
			style.Colors[ImGuiCol_BorderShadow]          = ImVec4(1.00f, 0.00f, 0.00f, 0.00f);
			style.Colors[ImGuiCol_FrameBg]               = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
			style.Colors[ImGuiCol_FrameBgHovered]        = ImVec4(0.90f, 0.80f, 0.80f, 0.40f);
			style.Colors[ImGuiCol_FrameBgActive]         = ImVec4(0.90f, 0.65f, 0.65f, 0.45f);
			style.Colors[ImGuiCol_TitleBg]               = ImVec4(0.960f, 0.960f, 0.960f, 1.0f);
			style.Colors[ImGuiCol_TitleBgCollapsed]      = ImVec4(0.960f, 0.960f, 0.960f, 1.0f);
			style.Colors[ImGuiCol_TitleBgActive]         = ImVec4(0.960f, 0.960f, 0.960f, 1.0f);
			style.Colors[ImGuiCol_MenuBarBg]             = ImVec4(0.01f, 0.01f, 0.02f, 0.80f);
			style.Colors[ImGuiCol_ScrollbarBg]           = ImVec4(0.20f, 0.25f, 0.30f, 0.60f);
			style.Colors[ImGuiCol_ScrollbarGrab]         = ImVec4(0.55f, 0.53f, 0.55f, 0.51f);
			style.Colors[ImGuiCol_ScrollbarGrabHovered]  = ImVec4(0.56f, 0.56f, 0.56f, 1.00f);
			style.Colors[ImGuiCol_ScrollbarGrabActive]   = ImVec4(0.56f, 0.56f, 0.56f, 0.91f);
			style.Colors[ImGuiCol_CheckMark]             = ImVec4(0.90f, 0.90f, 0.90f, 0.83f);
			style.Colors[ImGuiCol_SliderGrab]            = ImVec4(0.70f, 0.70f, 0.70f, 0.62f);
			style.Colors[ImGuiCol_SliderGrabActive]      = ImVec4(0.30f, 0.30f, 0.30f, 0.84f);
			style.Colors[ImGuiCol_Button]                = ImVec4(0.48f, 0.72f, 0.89f, 1.00f);
			style.Colors[ImGuiCol_ButtonHovered]         = ImVec4(0.50f, 0.69f, 0.99f, 1.00f);
			style.Colors[ImGuiCol_ButtonActive]          = ImVec4(0.80f, 0.50f, 0.50f, 1.00f);
			style.Colors[ImGuiCol_Header]                = ImVec4(0.44f, 0.61f, 0.86f, 1.00f);
			style.Colors[ImGuiCol_HeaderHovered]         = ImVec4(0.44f, 0.61f, 0.86f, 1.00f);
			style.Colors[ImGuiCol_HeaderActive]          = ImVec4(0.44f, 0.61f, 0.86f, 1.00f);
			style.Colors[ImGuiCol_ResizeGrip]            = ImVec4(0.70f, 0.70f, 0.70f, 1.00f);
			style.Colors[ImGuiCol_ResizeGripHovered]     = ImVec4(0.70f, 0.70f, 0.70f, 1.00f);
			style.Colors[ImGuiCol_ResizeGripActive]      = ImVec4(0.70f, 0.70f, 0.70f, 1.00f);
			style.Colors[ImGuiCol_PlotLines]             = ImVec4(0.00f, 1.00f, 0.00f, 1.00f);
			style.Colors[ImGuiCol_PlotLinesHovered]      = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
			style.Colors[ImGuiCol_PlotHistogram]         = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
			style.Colors[ImGuiCol_PlotHistogramHovered]  = ImVec4(1.00f, 0.60f, 0.00f, 1.00f);
			style.Colors[ImGuiCol_TextSelectedBg]        = ImVec4(0.00f, 0.00f, 1.00f, 0.35f);

			// Draw parent menu content
			menu.draw_viewer_menu();
		};

	g_viewer.core().background_color.setConstant(1.);
	g_viewer.core().is_animating = true;
	g_viewer.data().point_size = 3.; // determines visual appeal
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

/* 
Density Plot

Added from: https://github.com/dilevin/CSC417-a3-finite-elements-3d
*/ 
bool Visualize::plot_density(const char *label, unsigned int type, ImVec2 t_bounds, ImVec2 E_bounds, ImU32 plot_col) {
		using namespace ImGui; 

		unsigned int num_lines = 3; //always odd number because I want lines to cross at 0,0

		ImGuiContext& g = *GImGui;
		const ImGuiStyle& style = g.Style;

		//ImGUI stuff that I don't understand (taken from code example here: https://github.com/ocornut/imgui/issues/786)
		const ImGuiStyle& Style = GetStyle();
		const ImGuiIO& IO = GetIO();
		ImDrawList* DrawList = GetWindowDrawList();
		ImGuiWindow* Window = GetCurrentWindow();
	
		if (Window->SkipItems)
			return false;

		// header and spacing
		int hovered = IsItemActive() || IsItemHovered(); // IsItemDragged() ?
		Dummy(ImVec2(0,3));

		// prepare canvas
		ImVec2 Canvas = GetContentRegionAvail();
		Canvas.y = std::min(Canvas.y, 150.f);
		Canvas = CalcItemSize(Canvas, style.FramePadding.x * 2.0f, style.FramePadding.y * 2.0f);
		ImRect bb(Window->DC.CursorPos, Window->DC.CursorPos + Canvas);
		ItemSize(bb);
		const ImGuiID id = Window->GetID(label);
		
		RenderFrame(bb.Min, bb.Max, GetColorU32(ImGuiCol_FrameBg, 1), true, Style.FrameRounding);
		
		//local grid coordinates are -1,1 in both directions
		auto pix_to_normalized =  [&bb, &Canvas](ImVec2 pix) {  return ImVec2((pix.x-bb.Min.x)/Canvas.x,(pix.y-bb.Min.y)/Canvas.y); };
		auto normalized_to_pix =  [&bb, &Canvas] (ImVec2 norm) {  return ImVec2(norm.x*Canvas.x + bb.Min.x, Canvas.y - norm.y*Canvas.y + bb.Min.y); };
		auto data_to_normalized =  [&] (ImVec2 state) {  return ImVec2((state.x - g_density[0][0])/(t_bounds.y), (state.y - E_bounds.x)/(E_bounds.y - E_bounds.x)); };
		
		DrawList->AddText(bb.Min, GetColorU32(ImGuiCol_Text), label);

		//background grid centered on origin        
		for (float i = 0.f; i <= 1.f; i+= 1.f/static_cast<float>(num_lines-1)) {
			DrawList->AddLine(
				normalized_to_pix(ImVec2(0.f, i)),
				normalized_to_pix(ImVec2(1.f, i)),
				GetColorU32(ImGuiCol_TextDisabled), 1.2);
		}  

		//plot density
		if(g_density.size() == 0)
			return false;

		while ((g_density[g_density.size()-1][0] - g_density[0][0]) > t_bounds.y) {
			g_density.pop_front();
		}

		for(unsigned int ei=type; ei <= type; ++ei) {
			
			bool clip_p1;
			bool clip_p2;
			for(unsigned int i=0; i<g_density.size()-1; ++i) {
				clip_p1 = false;
				clip_p2 = false;

				ImVec2 p1 = data_to_normalized(ImVec2(g_density[i][0], g_density[i][ei]));
				ImVec2 p2 = data_to_normalized(ImVec2(g_density[i+1][0], g_density[i+1][ei]));

				if(p1.x < 0.f || p1.x > 1.f || p1.y < 0.f || p1.y > 1.f) {
					clip_p1 = true;
				}

				if(p2.x < 0.f || p2.x > 1.f || p2.y < 0.f || p2.y > 1.f) {
					clip_p2 = true;
				}

				p1.x = ImMin(ImMax(p1.x,0.f),1.f);
				p1.y = ImMin(ImMax(p1.y,0.f),1.f);
				p2.x = ImMin(ImMax(p2.x,0.f),1.f);
				p2.y = ImMin(ImMax(p2.y,0.f),1.f);

				if(!clip_p1 || !clip_p2) {
					DrawList->AddLine(
						normalized_to_pix(p1),
						normalized_to_pix(p2),
						plot_col, 2);
				}
			}
		}
		return true;
	}

void Visualize::add_density(float t, float average_density, float max_density) {
	//update plotting cache 
	std::array<float,3> tmp;
	tmp[0] = t;
	tmp[1] = average_density;
	tmp[2] = max_density;
	g_density.push_back(tmp);
}

void Visualize::write_frame(std::string &path){
	const int width  = g_viewer.core().viewport(2);
    const int height = g_viewer.core().viewport(3);

    std::unique_ptr<GLubyte[]> pixels(new GLubyte[width * height * 4]);
    glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels.get());
    
    igl::stbi_write_png(path.c_str(), width, height, 4, pixels.get() + width * (height - 1) * 4, - width * 4);
}