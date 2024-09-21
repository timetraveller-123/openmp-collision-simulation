#include <omp.h>

#include <cmath>
#include <fstream>
#include <vector>
#include <utility>
#include <iostream> 
#include "collision.h"
#include "io.h"
#include "sim_validator.h"

void update_positions(std::vector<Particle> &particles) {
    for(int i = 0; i < (int)particles.size();i++) {
        particles[i].loc.x += particles[i].vel.x;
        particles[i].loc.y += particles[i].vel.y;
    }

}

void print(std::vector<Particle> particles) {
    for(auto i:particles) {
        std::cout<<i.i<<" "<<i.loc.x<<" "<<i.loc.y<<" "<<i.vel.x<<" "<<i.vel.y<<std::endl;
    }
    std::cout<<std::endl<<std::endl;
}

void detect_overlaps(std::vector<Particle> &particles,
                    std::vector<std::pair<int,int>> &overlaps, int square_size, int radius) {

    for(int i = 0; i < (int)particles.size(); i++) {
        if(is_wall_overlap(particles[i].loc, square_size, radius)){
            overlaps.push_back(std::make_pair(i,-1));
        }
        for(int j = i + 1; j < (int)particles.size(); j++) {
            if(is_particle_overlap(particles[i].loc, particles[j].loc, radius)) {
                overlaps.push_back(std::make_pair(i,j));
            }
        }
    }

}


bool detect_collisions(std::vector<Particle> &particles, 
                        std::vector<std::pair<int,int>> &overlaps, int square_size, int radius) {
    bool no_collision = true;
    for(int i = 0; i < (int)overlaps.size();i++ ) {
        int f = overlaps[i].first;
        int s = overlaps[i].second;
        if(s == -1 && is_wall_collision(particles[f].loc, particles[f].vel, square_size, radius)) {
            resolve_wall_collision(particles[f].loc, particles[f].vel, square_size, radius);        
            no_collision = false;
        } else {
            if(is_particle_collision(particles[f].loc, particles[f].vel, particles[s].loc, particles[s].vel, radius)) {
                resolve_particle_collision(particles[f].loc, particles[f].vel, particles[s].loc, particles[s].vel);
                no_collision = false;
            }
        }

    }
    return no_collision;
}

int main(int argc, char* argv[]) {
    // Read arguments and input file
    Params params{};
    std::vector<Particle> particles;
    read_args(argc, argv, params, particles);

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
    
    std::vector<std::pair<int,int>> overlaps;
    for(int i = 1; i <= params.param_steps; i ++) {
        update_positions(particles);
        detect_overlaps(particles, overlaps, params.square_size, params.param_radius);
        while(!detect_collisions(particles, overlaps, params.square_size, params.param_radius)) {
            
        }
        #if CHECK == 1
            // Check final positions
            validator.validate_step(particles);
        #endif
    }

}
