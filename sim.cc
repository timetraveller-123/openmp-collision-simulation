#include <omp.h>

#include <cmath>
#include <fstream>
#include <vector>
#include <utility>
#include <iostream> 
#include "collision.h"
#include "io.h"
#include "sim_validator.h"
#include <algorithm>
#include <chrono>

void update_positions(std::vector<Particle> &particles) {
    for(int i = 0; i < (int)particles.size();i++) {
        particles[i].loc.x += particles[i].vel.x;
        particles[i].loc.y += particles[i].vel.y;
    }

}



void assign_grid(std::vector<Particle> particles, std::vector<int> grid[], std::vector<int> &particle_map, int number_of_cells, int square_size) {
    double cell_length = (double)square_size/number_of_cells;
    for(int i = 0; i < (int)particles.size(); i++) {
        int x = std::min(std::max(0,(int)((particles[i].loc.x)/cell_length)),number_of_cells - 1);
        int y = std::min(std::max(0,(int)((particles[i].loc.y)/cell_length)), number_of_cells - 1);
        grid[x + y*number_of_cells].push_back(i);
        particle_map[i] = x + y*number_of_cells;
    }

}

void clear_grid(std::vector<int> grid[], int number_of_cells) {
    for(int i = 0; i < number_of_cells*number_of_cells; i++) {
        grid[i].clear();
    }
}


void print(std::vector<Particle> particles) {
    for(auto i:particles) {
        std::cout<<i.i<<" "<<i.loc.x<<" "<<i.loc.y<<" "<<i.vel.x<<" "<<i.vel.y<<std::endl;
    }
    std::cout<<std::endl<<std::endl;
}

void print_grid(std::vector<int> grid[], int number_of_cells) {
    for(int i = 0; i<number_of_cells*number_of_cells;i++) {
        std::cout<<i<<":   ";
        for(auto j:grid[i]) {
            std::cout<<j<< " ";
        }
        std::cout<<std::endl;
        
    }
}

void detect_overlaps(std::vector<Particle> &particles, std::vector<int> grid[], std::vector<std::vector<int>> &overlaps, std::vector<int> &particle_map, int number_of_cells, int radius) {
    for(int i = 0; i < (int)particles.size(); i++) {
        overlaps[i].clear();

        for(auto j : grid[particle_map[i]]) {
            if(is_particle_overlap(particles[i].loc, particles[j].loc, radius)) {
                overlaps[i].push_back(j);
            }
        }

        int dirx[] = {1,1,0,-1,-1,-1,0,1};
        int diry[] = {0,number_of_cells,number_of_cells,number_of_cells,0,-number_of_cells,-number_of_cells,-number_of_cells};

        for(int d = 0; d < 8; d++) {
            int h = particle_map[i] +dirx[d] + diry[d];
            if(h < number_of_cells*number_of_cells && h >= 0 && (h%number_of_cells + dirx[d]) >= 0 && (h%number_of_cells + dirx[d]) < number_of_cells) {
                for(auto j : grid[h]) {
                    if(is_particle_overlap(particles[i].loc, particles[j].loc, radius)) {
                        overlaps[i].push_back(j);
                    }
                }
            }
        }


    }
}



bool resolve_collisions(std::vector<Particle> &particles, std::vector<std::vector<int>> &overlaps, int square_size, int radius) {

    bool no_collisions = true;

    for(int i = 0; i < (int)particles.size(); i++) {
        std::vector<int> collisions;


        if(is_wall_collision(particles[i].loc, particles[i].vel, square_size, radius)) {
            no_collisions = false;
            collisions.push_back(-1);
        }


        for(auto j: overlaps[i]) {
            if(is_particle_collision(particles[i].loc, particles[i].vel, particles[j].loc, particles[j].vel, radius)) {
                 no_collisions = false;
                 collisions.push_back(j);
            }
        }

        
        for(auto j: collisions) {
            if(j == -1) {
                resolve_wall_collision(particles[i].loc, particles[i].vel, square_size, radius);
            } else {
                resolve_particle_collision(particles[i].loc, particles[i].vel, particles[j].loc, particles[j].vel);
            }
        }
        
    }
    return no_collisions;
}



int main(int argc, char* argv[]) {
    // Read arguments and input file
    Params params{};
    std::vector<Particle> particles;
    read_args(argc, argv, params, particles);

    int time = 0, time2 = 0;

    // Set number of threads
    omp_set_num_threads(params.param_threads);

#if CHECK == 1
    // Initialize collision checker
    SimulationValidator validator(params.param_particles, params.square_size, params.param_radius);
    // Initialize with starting positions
    validator.initialize(particles);
    // Uncomment the line below to enable visualization (makes program much slower)
    validator.enable_viz_output("test.out");
#endif

    // TODO: this is the part where you simulate particle behavior.

    /*
    After simulating each timestep, you must call:

    #if CHECK == 1
    validator.validate_step(particles);
    #endif
    */
    
    std::vector<std::vector<int>> overlaps;
    std::vector<int> particle_map(params.param_particles, 0);
    const int number_of_cells = 4;
    std::vector<int> grid[number_of_cells*number_of_cells];

    for(int i = 0; i < params.param_particles; i++) {
        std::vector<int> v;
        overlaps.push_back(v);
    }


    for(int i = 1; i <= params.param_steps; i ++) {
        //std::cout<<i<<std::endl;
        auto begin = std::chrono::high_resolution_clock::now();        
        update_positions(particles);
        assign_grid(particles, grid, particle_map, number_of_cells, params.square_size);
        detect_overlaps(particles, grid, overlaps, particle_map, number_of_cells, params.param_radius);
        const auto last = std::chrono::high_resolution_clock::now();
        time2 += (std::chrono::duration_cast<std::chrono::microseconds>(last - begin )).count();




        auto start = std::chrono::high_resolution_clock::now();        
        while(!resolve_collisions(particles, overlaps, params.square_size, params.param_radius)) {
        }
        const auto end = std::chrono::high_resolution_clock::now();
        time += (std::chrono::duration_cast<std::chrono::microseconds>(end - start )).count();

        
        #if CHECK == 1
            // Check final positions
            validator.validate_step(particles);
        #endif
        clear_grid(grid, number_of_cells);
        
    }

    std::cout<<time<<" " <<time2<<std::endl;

}
