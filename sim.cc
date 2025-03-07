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
#include <queue>
#include <random>

omp_lock_t locks[1000000];

std::vector<std::pair<int,int>> overlaps[9][1000000], overlap[9][1000000];

std::vector<int> grid[1000000];
std::vector<Particle> particles;

int num;

int number_of_cells = 0, n=0, group = 0;
int square_size = 0;
int radius = 0;
int gen[1000000];

long long assign_grid_time = 0,total_time=0,  overlap_time=0, collision_time =0, start_time = 0;

 

void update_positions() {
    int i; int s = (int)particles.size();
    #pragma omp parallel for shared(particles) private(i) num_threads(num)
    for(i = 0; i < s;i++) {
        Particle p = particles[i];
        particles[i].loc.x += p.vel.x;
        particles[i].loc.y += p.vel.y;
    }

}

void clear_grid() {
    int i;
    #pragma omp parallel for shared(grid) private(i) num_threads(num)
    for(i = 0; i < number_of_cells*number_of_cells;i++){
        grid[i].clear();
    }
}


void assign_grid() {
    clear_grid();
    int i;
    double cell_length = (double)square_size/number_of_cells;
    #pragma omp parallel for private(i)
    for(i = 0; i < (int)particles.size(); i++) {
        int x = std::min(std::max(0,(int)((particles[i].loc.x)/cell_length)),number_of_cells - 1);
        int y = std::min(std::max(0,(int)((particles[i].loc.y)/cell_length)), number_of_cells - 1);
        
        omp_set_lock(&locks[x + y*number_of_cells]);
        grid[x + number_of_cells*y].push_back(i);
        omp_unset_lock(&locks[x + y*number_of_cells]);

        
    }
}

void clear_overlap() {
    for(int d = 0; d < 9; d++){
        for(int i = 0; i < number_of_cells*number_of_cells; i++){
            overlap[d][i].clear();
        }
    }
}


void detect_overlaps() {
    
    int dirx[] = {0,1,1,0,-1,-1,-1,0,1};
    int diry[] = {0,0,1,1,1,0,-1,-1,-1};
    
    
    int k;
    clear_overlap();
    //gridwise overlap
    #pragma omp parallel for shared(overlap) private(k)  num_threads(num)
    for(k = 0; k < n*n; k++){
        for(int j = 0; j < group*group; j++){
            int i = k/n*number_of_cells*group + k%n*group + gen[j];
            if(i >= number_of_cells*number_of_cells || (i%number_of_cells)/group + (i/number_of_cells)/group*n != k)continue;
            for(int d = 0; d < 9; d++){
                int ni = i + dirx[d] + diry[d]*number_of_cells;
                if(ni  < 0 || ni >= number_of_cells*number_of_cells || (i%number_of_cells + dirx[d])  < 0 || (i%number_of_cells + dirx[d])  >= number_of_cells)continue;
                int dx = (ni%number_of_cells)/group -(i%number_of_cells)/group;
                int dy = (ni/number_of_cells)/group - (i/number_of_cells)/group;
                int dir;
                for( dir = 0; dir < 9;dir++){
                    if(dirx[dir] == dx && diry[dir] == dy)break;
                }
                
                for(auto m:grid[ni]){
                    for(auto l:grid[i]){

                        if(is_particle_overlap(particles[l].loc, particles[m].loc, radius)) {
                            overlap[dir][k].push_back(std::make_pair(l,m));
                          
                        }

                    }
                }
            }
        }
    }

     
}

bool resolve_collisions() {
    

    
    bool no_collisions = true;
    int i;
    #pragma omp parallel for  private(i) num_threads(num)
    for(i = 0; i < (int)particles.size(); i++) {
        if(is_wall_collision(particles[i].loc, particles[i].vel, square_size, radius)) {
            no_collisions = false;
            resolve_wall_collision(particles[i].loc, particles[i].vel, square_size, radius);
        }
    }

    int dirx[] = {0,1,1,0,-1,-1,-1,0,1};
    int diry[] = {0,0,1,1,1,0,-1,-1,-1};
    for(int d = 0; d < 9; d++) {
        
        int i;
        #pragma omp parallel for shared(overlap) private(i) num_threads(num)
        for(i = 0; i < n*n; i++){
            if((i*dirx[d] + (i/n)*diry[d])%2 == 1)continue;
            
                for(auto p: overlap[d][i]){
                    if(is_particle_collision(particles[p.first].loc, particles[p.first].vel, particles[p.second].loc, particles[p.second].vel, radius)) {
                        no_collisions = false;
                        resolve_particle_collision(particles[p.first].loc, particles[p.first].vel, particles[p.second].loc, particles[p.second].vel);
                    }    
                }

            
        }

        #pragma omp parallel for shared(overlap) private(i) num_threads(num)
        for(i = 0; i < n*n; i++){
            if((i*dirx[d] + (i/n)*diry[d])%2 == 0)continue;
            
                for(auto p: overlap[d][i]){
                    if(is_particle_collision(particles[p.first].loc, particles[p.first].vel, particles[p.second].loc, particles[p.second].vel, radius)) {
                        no_collisions = false;
                        resolve_particle_collision(particles[p.first].loc, particles[p.first].vel, particles[p.second].loc, particles[p.second].vel);
                    }    
                }

            
        }
    }
    return no_collisions;
}







int main(int argc, char* argv[]) {
    auto start_start = std::chrono::high_resolution_clock::now();
    // Read arguments and input file
    Params params{};
    
    read_args(argc, argv, params, particles);
    
    // Set number of threads
    omp_set_num_threads(params.param_threads);
    num = 8;

#if CHECK == 1
    // Initialize collision checker
    SimulationValidator validator(params.param_particles, params.square_size, params.param_radius);
    // Initialize with starting positions
    validator.initialize(particles);
    // Uncomment the line below to enable visualization (makes program much slower)
    //validator.enable_viz_output("test.out");
#endif

    
    
    number_of_cells = std::min(1000, (int)((params.square_size)/(2*params.param_radius)));
    group = std::max(1,std::min(10, (int)number_of_cells/5));
    n = (number_of_cells - 1)/group + 1;
    square_size = params.square_size;
    radius = params.param_radius;
    
    
    for(int i = 0; i<group*group;i++){
        gen[i] = i%group + i/group*number_of_cells;
    }

    for(int i = 0; i < n*n;i++){
        omp_init_lock(&locks[i]);
    }

    auto start_end = std::chrono::high_resolution_clock::now();
    start_time += std::chrono::duration_cast<std::chrono::nanoseconds>(start_end - start_start).count();


    for(int i = 1; i <= params.param_steps; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        int iter = 0;

        

        auto assign_start = std::chrono::high_resolution_clock::now();
        update_positions();

        assign_grid();
        auto assign_end = std::chrono::high_resolution_clock::now();
        assign_grid_time += std::chrono::duration_cast<std::chrono::nanoseconds>(assign_end - assign_start).count();

        
        auto overlap_start = std::chrono::high_resolution_clock::now();
        detect_overlaps();
        auto overlap_end = std::chrono::high_resolution_clock::now();
        overlap_time += std::chrono::duration_cast<std::chrono::nanoseconds>(overlap_end - overlap_start).count();
        
        

        auto collision_start = std::chrono::high_resolution_clock::now();
        while(!resolve_collisions()) {iter++;
        }
        auto collision_end = std::chrono::high_resolution_clock::now();
        collision_time += std::chrono::duration_cast<std::chrono::nanoseconds>(collision_end - collision_start).count();


        auto end = std::chrono::high_resolution_clock::now();
        long long t = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        total_time += t;

        //std::cout<<i<<" "<<t/1000000<<" "<<iter<<std::endl;

        
        #if CHECK == 1
            // Check final positions
            validator.validate_step(particles);
        #endif
        
    }
    std::cout<<total_time/1000000<<" "<<assign_grid_time/1000000<<" "<<overlap_time/1000000<<" "<<collision_time/1000000<< " " <<start_time/1000000<<" "<<std::endl;

}