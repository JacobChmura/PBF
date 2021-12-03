#ifndef SPATIALHASHGRID
#define SPATIALHASHGRID

#include <Eigen/Dense>
#include <math.h>
#include <vector>
#include <map>
#include <set>
#include <Particle.h>
#include <utility> // pair
class SpatialHashGrid{

public: 
        // Initialize a Spacial Hash Grid with a given boundary and granularity
        SpatialHashGrid(double l_bound, double u_bound, double cell_size, int num_particles);
        SpatialHashGrid(); // default constructor

        // Insert a particle into our grid
        void insert(Particle &p);

        // Remove a particle from our grid 
        void remove(Particle &p);

        // Update the grid map
        void update(Particle &p);

        // Find Neighours
        void findNeighbours(Particle &p);
        
private:
        
        // Maps a cell hash to a list of particle global indices
        std::map<std::tuple<int, int, int>, std::set<int>> cells;
        Eigen::Vector3d lower_bound, upper_bound;
        double cell_size;
        int table_size;

        // Convert World Coordiantes to Cell Coordinates
        Eigen::Vector3d WorldToCell(Eigen::Vector3d world_coord);

        // Hash map cell coordinates to a bucket in our hash map
        std::tuple<int, int, int> hash(Eigen::Vector3d cell_coord);
};



#endif
