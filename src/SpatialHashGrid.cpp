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

void SpatialHashGrid::insert(Particle &p){
        // Convert from world coordiantes to cell coordinates and hash
        Eigen::Vector3d cell_coord = this->WorldToCell(p.x);
        p.cell_coord = cell_coord;
        std::tuple<int, int, int> hashed_coord = this->hash(cell_coord);

        if (this->cells.find(hashed_coord) == this->cells.end()){
                // this key does not exist yet
                this->cells[hashed_coord] = std::set<int>(); 
        }
        // Add the particle global index
        this->cells[hashed_coord].insert(p.global_idx);
              
        //std::cout << "Particle: " << p.global_idx << "\nWorld Coord: " << p.x << "\nCell Coord: " << cell_coord << "\nHash: " << std::get<0>(hashed_coord) << ", " << std::get<1>(hashed_coord) << ", " << std::get<2>(hashed_coord) << "\n\n";
}

void SpatialHashGrid::remove(Particle &p){
        std::tuple<int,int,int> hashed_coord = this->hash(p.cell_coord);
        this->cells.erase(hashed_coord);
}

void SpatialHashGrid::update(Particle &p){
        this->remove(p);
        this->insert(p);
}

void SpatialHashGrid::findNeighbours(Particle &p){
        p.neighbours.clear(); // clear neighbours from previous iteration
       
        // Store offset cell coordinate 
        Eigen::Vector3d neighbourhood_coord;
        std::tuple<int, int, int> hashed_coord;

        // We use a cell granularity of (kernel_width) and loop through the neighbouring grid cells (conversative estimate on neighbours) 
        Eigen::Vector3d cell_coord = this->WorldToCell(p.x_new); 

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

                                        // Add each one (not including the particle itself) as a neighbour
                                        for (auto particle_idx: cells[hashed_coord]){
                                                if (particle_idx != p.global_idx) p.neighbours.insert(particle_idx);
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
