#ifndef SPATIALHASHGRID
#define SPATIALHASHGRID

#include <Eigen/Dense>
#include <math.h>
#include <map>
#include <set>
#include <tuple>
#include <vector>
#include <Particle.h>

class SpatialHashGrid{

public: 
        /*
        Initialize a new Spatial Hash Grid for efficient neighbourhood queries.

        Input:
                double l_bound: the minimum world space coord. (assumed the same for x, y, z)
                double u_bound: the maximum world space coord. (assumed the same for x, y, z)
                double cell_size: the distance in world-space that each cell in our grid occupies.
                                  (in our case this is equal to the kernel size)
        */
        SpatialHashGrid(); // default constructor
        SpatialHashGrid(double l_bound, double u_bound, double cell_size);

        /* Insert the particle p into the grid. Requires finding the cell that it is contained in, 
         * then adding the particle global index into this hash.

        Input: 
                Particle p: the particle to insert.
        */
        void insert(const Eigen::Ref<const Eigen::MatrixXd> &fluid_state);

        /* After each simulation step, we need to update our cell grid. Currently this is done by 
         * removing all particles from the grid and inserting them again. This should be optimized.

        Input: 
                Particle p: the particle to update.
        */
        void update(const Eigen::Ref<const Eigen::MatrixXd> &fluid_state);

        /* Find Neighours for particle p using the predicted position for the jacobi update, 
         * and update the p.neighbours set to reflect the (possibly) new neighbours.
        
        Input:
                Particle p: the particle to find neighbours for.
        */  
        void findNeighbours(const Eigen::Ref<const Eigen::MatrixXd> &x_new, std::vector<std::vector<int>> &neighbours);
        
private:
        
        // Grid Dimensionality and Size Properties
        Eigen::Vector3d lower_bound, upper_bound;
        double cell_size; 
        int range; // upper bound - lower bound 

        /*
        This is the hash map. 
        <key> A three-tuple hash which uniquely identifies each 3d cell in our grid.
        <value> A set of global indices for all the particls currently in that cell.
        */
        std::map<std::tuple<int, int, int>, std::set<int>> cells;

        /*
        Convert world-space coordinates into a cell index in the x, y, z axis.

        Input:
                Eigen::Vector3d world_coord: the 3d position in world space.

        Output:
                Eigen::Vector3d the positional indices in our 3d cell grid.
        */ 
        Eigen::Vector3d WorldToCell(Eigen::Vector3d world_coord);

        /* 
        Hash the positional indices in our 3d cell grid to a key in our hash map.
        Right now, our keys are just the tuple of positional indices.

        Input:
                Eigen::Vector3d cell_coord the positional indices in our 3d cell grid.

        Output:
                std::tuple<int, int, int> hash: contains key into our cell hash map by
                                                simply converting the input to a tuple. 
        */
        std::tuple<int, int, int> hash(Eigen::Vector3d cell_coord);
};

#endif
