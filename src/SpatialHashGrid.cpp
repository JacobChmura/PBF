#include <SpatialHashGrid.h>

SpatialHashGrid::SpatialHashGrid(){
}

SpatialHashGrid::SpatialHashGrid(double l_bound, double cell_size, int num_particles){
        this->lower_bound << l_bound, l_bound, l_bound;
        this->cell_size = cell_size;

        this->table_size = 2 * num_particles;
        // this->P1 = 73856093;
        // this->P2 = 19349663;
        // this->P3 = 83492791;
}

void SpatialHashGrid::insert(Particle &p){
        // Convert from world coordiantes to cell coordinates and hash
        Eigen::Vector3d cell_coord = this->WorldToCell(p.x);
        p.cell_coord = cell_coord;
        int hashed_coord = this->hash(cell_coord);

        if (this->cells.find(hashed_coord) == this->cells.end()){
                // this key does not exist yet
                this->cells[hashed_coord] = std::set<int>(); 
        }
        // add
        this->cells[hashed_coord].insert(p.global_idx);
              
}

void SpatialHashGrid::remove(Particle &p){
        int hashed_coord = this->hash(p.cell_coord);
        this->cells.erase(hashed_coord);
}

void SpatialHashGrid::update(Particle &p){
        this->remove(p);
        this->insert(p);
}

void SpatialHashGrid::findNeighbours(Particle &p){
       p.neighbours.clear(); // clear neighbours from previous iteration

       Eigen::Vector3d cell_coord = this->WorldToCell(p.x); 
       int hashed_coord = this->hash(cell_coord);

        if (this->cells.find(hashed_coord) != this->cells.end()){
                for (auto particle_idx: cells[hashed_coord]){
                        p.neighbours.insert(particle_idx);
                }
        }
}

Eigen::Vector3d SpatialHashGrid::WorldToCell(Eigen::Vector3d world_coord){
        Eigen::Vector3d cell_coord;
        cell_coord[0] = floor(world_coord[0] - this->lower_bound[0] / this->cell_size);        
        cell_coord[1] = floor(world_coord[1] - this->lower_bound[1] / this->cell_size);        
        cell_coord[2] = floor(world_coord[2] - this->lower_bound[2] / this->cell_size);        
        
        return cell_coord;
}


int SpatialHashGrid::hash(Eigen::Vector3d cell_coord){
        return (int(cell_coord[0]) * this->P1 ^ int(cell_coord[1]) * this->P2 ^ int(cell_coord[2]) * this->P2) % this->table_size;
}
