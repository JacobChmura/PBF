#include <SpatialHashGrid.h>
#include <iostream>

SpatialHashGrid::SpatialHashGrid(){
}

SpatialHashGrid::SpatialHashGrid(double l_bound, double u_bound, double cell_size){
        this->lower_bound << l_bound, l_bound, l_bound;
        this->upper_bound << u_bound, u_bound, u_bound; 
        this->cell_size = cell_size;
        this->range = u_bound - l_bound;
}

void SpatialHashGrid::insert(const Eigen::Ref<const Eigen::MatrixXd> &fluid_state){
        Eigen::Vector3d cell_coord;
        std::tuple<int, int, int> hashed_coord;
        for (int i = 0; i < fluid_state.rows(); i++){
                cell_coord = this->WorldToCell(fluid_state.row(i));
                hashed_coord = this->hash(cell_coord);

                // this key does not exist yet, create an empty set at that location
                if (this->cells.find(hashed_coord) == this->cells.end()) this->cells[hashed_coord] = std::set<int>(); 
                // Add the particle global index
                this->cells[hashed_coord].insert(i);
        }
}


void SpatialHashGrid::update(const Eigen::Ref<const Eigen::MatrixXd> &fluid_state){
        for(auto &element: this->cells){
                element.second.clear();
        }
        this->insert(fluid_state); 
}

void SpatialHashGrid::findNeighbours(const Eigen::Ref<const Eigen::MatrixXd> &x_new, std::vector<std::vector<int>> &neighbours){
        // Store offset cell coordinate 
        Eigen::Vector3d neighbourhood_coord, cell_coord;
        std::tuple<int, int, int> hashed_coord;
        
        for (int i = 0; i < x_new.rows(); i++){
                // Clear neighbours from previous iteratoin
                neighbours[i].clear(); 

                // We use a cell granularity of (kernel_width) and loop through the neighbouring grid cells (conversative estimate on neighbours) 
                cell_coord = this->WorldToCell(x_new.row(i)); 

                // loop through 8 neighbourhing cells to find neighbours
                for (int dx = -1; dx <= 1; dx ++){
                        for (int dy = -1; dy <= 1; dy++){
                                for (int dz = -1; dz <= 1; dz ++){
                                        // Construct cell coordinate of the current neighbourhood 
                                        neighbourhood_coord << dx, dy, dz;
                                        neighbourhood_coord += cell_coord;
                                        hashed_coord = this->hash(neighbourhood_coord);
                                        
                                        // Check if there are particles in this hash
                                        if (this->cells.find(hashed_coord) != this->cells.end()){

                                                // Add each one (including yourself) as a neighbour
                                                for (auto particle_idx: cells[hashed_coord]){
                                                        if (! std::count(neighbours[i].begin(), neighbours[i].end(), particle_idx)) neighbours[i].push_back(particle_idx);
                                                        
                                                }
                                        }
                                }
                        } 
                }
        }


}

Eigen::Vector3d SpatialHashGrid::WorldToCell(Eigen::Vector3d world_coord){
        // Coordinates between 0 and 1 within grid
        Eigen::Vector3d cell_coord = (world_coord - this->lower_bound) / this->range;

        // Cell number in each axis
        cell_coord *= ((this->range / this->cell_size) - 1);
        return cell_coord;
}

std::tuple<int, int, int> SpatialHashGrid::hash(Eigen::Vector3d cell_coord){
        return std::make_tuple(cell_coord[0], cell_coord[1], cell_coord[2]); 

}
