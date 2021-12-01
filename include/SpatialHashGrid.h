#ifndef SPATIALHASHGRID
#define SPATIALHASHGRID

#include <Eigen/Dense>
#include <math.h>
#include <vector>
#include <map>
#include <set>
#include <Particle.h>

class SpatialHashGrid{

public: 
        // Initialize a Spacial Hash Grid with a given boundary and granularity
        SpatialHashGrid(double l_bound, double cell_size, int num_particles);
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
        std::map<int, std::set<int>> cells;
        Eigen::Vector3d lower_bound;
        double cell_size;

        // Convert World Coordiantes to Cell Coordinates
        Eigen::Vector3d WorldToCell(Eigen::Vector3d world_coord);

        // Hash map cell coordinates to a bucket in our hash map
        int hash(Eigen::Vector3d cell_coord);

        // Hash functions large primes
        int table_size;
        int P1 = 73856093;
        int P2 = 19349663;
        int P3 = 83492791;

};



#endif
