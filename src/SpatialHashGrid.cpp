#include <SpatialHashGrid.h>
#include <iostream>

SpatialHashGrid::SpatialHashGrid(){
}

SpatialHashGrid::SpatialHashGrid(double l_bound, double u_bound, double cell_size, int num_particles){
        this->lower_bound << l_bound, l_bound, l_bound;
        this->upper_bound << u_bound, u_bound, u_bound;
        this->cell_size = cell_size;
        this->table_size = 2 * num_particles;
}

void SpatialHashGrid::insert(Particle &p){
        // Convert from world coordiantes to cell coordinates and hash
        Eigen::Vector3d cell_coord = this->WorldToCell(p.x);

        p.cell_coord = cell_coord;
        std::tuple<int, int, int> hashed_coord = this->hash(cell_coord);

        //std::cout << "Particle: " << p.global_idx << "\nWorld Coord: " << p.x << "\nCell Coord: " << cell_coord << "\nHash: " << std::get<0>(hashed_coord) << ", " << std::get<1>(hashed_coord) << ", " << std::get<2>(hashed_coord) << "\n\n";
        if (this->cells.find(hashed_coord) == this->cells.end()){
                // this key does not exist yet
                this->cells[hashed_coord] = std::set<int>(); 
        }
        // add
        this->cells[hashed_coord].insert(p.global_idx);
              
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

        Eigen::Vector3d cell_coord = this->WorldToCell(p.x_new); 
        Eigen::Vector3d neighbourhood_coord;
        std::tuple<int, int, int> hashed_coord;

        // loop through 8 neighbourhing cells to find neighbours
        for (int dx = -1; dx <= 1; dx ++){
                for (int dy = -1; dy <= 1; dy++){
                        for (int dz = -1; dz <= 1; dz ++){
                                
                                neighbourhood_coord << dx, dy, dz;
                                neighbourhood_coord += cell_coord;
                                neighbourhood_coord = neighbourhood_coord.cwiseMax(0).cwiseMin(this->table_size);
       
                                hashed_coord = this->hash(neighbourhood_coord);
                                if (this->cells.find(hashed_coord) != this->cells.end()){
                                        for (auto particle_idx: cells[hashed_coord]){
                                                p.neighbours.insert(particle_idx);
                                        }
                                }
                        }
                } 
        }
}

Eigen::Vector3d SpatialHashGrid::WorldToCell(Eigen::Vector3d world_coord){
        Eigen::Vector3d cell_coord;
        // Get value between 0 and 1
        cell_coord[0] = (world_coord[0] - this->lower_bound[0]) / (this->upper_bound[0] - this->lower_bound[0]);
        cell_coord[1] = (world_coord[1] - this->lower_bound[1]) / (this->upper_bound[1] - this->lower_bound[1]);
        cell_coord[2] = (world_coord[2] - this->lower_bound[2]) / (this->upper_bound[2] - this->lower_bound[2]);

        cell_coord[0] = floor(cell_coord[0] * ((this->upper_bound[0] - this->lower_bound[0] / this->cell_size) - 1));
        cell_coord[1] = floor(cell_coord[1] * ((this->upper_bound[1] - this->lower_bound[1] / this->cell_size) - 1));
        cell_coord[2] = floor(cell_coord[2] * ((this->upper_bound[2] - this->lower_bound[2] / this->cell_size) - 1));
        return cell_coord;
}

std::tuple<int, int, int> SpatialHashGrid::hash(Eigen::Vector3d cell_coord){
        return std::make_tuple(cell_coord[0], cell_coord[1], cell_coord[2]); 

}
