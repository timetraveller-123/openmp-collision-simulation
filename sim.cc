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

omp_lock_t locks[1000000];
std::vector<std::pair<int,int>> connected_overlap[1000000];
std::vector<std::pair<int,int>> overlap;
//std::vector<int> comp[1000000];
bool visit[1000000];
int mi = 10000000;
long long ta = 0,ma = 0;

long long assign_grid_time = 0,total_time=0, bfs_time=0, overlap_time=0, collision_time =0, start_time = 0;

void update_positions(std::vector<Particle> &particles) {
    int i; int s = (int)particles.size();
    #pragma omp parallel for shared(particles) private(i)
    for(i = 0; i < s;i++) {
        Particle p = particles[i];
        particles[i].loc.x += p.vel.x;
        particles[i].loc.y += p.vel.y;
    }

}

void clear_grid(std::vector<std::vector<int>> &grid, int number_of_cells) {
    int i;
    #pragma omp parallel for shared(grid) private(i)
    for(i = 0; i < number_of_cells*number_of_cells; i++) {
        grid[i].clear();
    }
}


void assign_grid(std::vector<Particle> particles, std::vector<std::vector<int>> &grid, std::vector<int> &particle_map, int number_of_cells, int square_size) {
    clear_grid(grid, number_of_cells);

    double cell_length = (double)square_size/number_of_cells;
    int i;
    #pragma omp parallel for shared(grid) private(i)
    for(i = 0; i < (int)particles.size(); i++) {
        int x = std::min(std::max(0,(int)((particles[i].loc.x)/cell_length)),number_of_cells - 1);
        int y = std::min(std::max(0,(int)((particles[i].loc.y)/cell_length)), number_of_cells - 1);
        
        omp_set_lock(&locks[x + y*number_of_cells]);
        grid[x + y*number_of_cells].push_back(i);
        omp_unset_lock(&locks[x + y*number_of_cells]);

        particle_map[i] = x + y*number_of_cells;
    }
}


void detect_overlaps(std::vector<Particle> &particles, std::vector<std::vector<int>> &grid, std::vector<std::vector<int>> &overlaps, std::vector<int> &particle_map, int number_of_cells, int radius) {
    

    int i;
    grid = grid;
    particle_map = particle_map;
    number_of_cells = number_of_cells;
    #pragma omp parallel for shared(overlaps) private(i) 
    for(i = 0; i < (int)particles.size(); i++) {
        overlaps[i].clear();

        int dirx[] = {0, 1,1,0,-1,-1,-1,0,1};
        int diry[] = {0, 0,number_of_cells,number_of_cells,number_of_cells,0,-number_of_cells,-number_of_cells,-number_of_cells};

        for(int d = 0; d < 9; d++) {
            int h = particle_map[i] +dirx[d] + diry[d];
            if(h < number_of_cells*number_of_cells && h >= 0 && (particle_map[i]%number_of_cells + dirx[d]) >= 0 && (particle_map[i]%number_of_cells + dirx[d]) < number_of_cells) {
                for(auto j : grid[h]) {
                    
                    if(j != i && is_particle_overlap(particles[i].loc, particles[j].loc, radius)) {
                        overlaps[i].push_back(j);
                    }
                }
            }
        }

    }
    
}

void bfs(std::vector<Particle> &particles, std::vector<std::vector<int>> &overlaps, int s, int index, int square_size, int radius) {

    std::queue<int> q;
    connected_overlap[index].clear();
    
    s =s;
    particles = particles;
    square_size = square_size;
    radius = radius; 

    q.push(s);
    visit[s] = true;
    while(!q.empty()) {
        int n = q.front();
        
        q.pop();
        for(auto j:overlaps[n]) {
            connected_overlap[index].push_back(std::make_pair(n,j));
            if(!visit[j]) {
                visit[j] = true;
                q.push(j);
            }
        }
    }
    ta+= (long long)connected_overlap[index].size();
    
    ma = std::max(ma, (long long)connected_overlap[index].size());
}




bool resolve_collisions(std::vector<Particle> &particles, int index , int square_size, int radius) {
    bool no_collisions = true; 
    square_size = square_size;
    int i;
    #pragma omp parallel for shared(particles) private(i) 
    for(i = 0; i < (int)particles.size(); i++) {
        if(is_wall_collision(particles[i].loc, particles[i].vel, square_size, radius)) {
            no_collisions = false;
            resolve_wall_collision(particles[i].loc, particles[i].vel, square_size, radius);
        }
    }

    #pragma omp parallel for shared(particles) private(i) schedule(dynamic)
    for(i = 0; i < index; i++) {
        //if((long long)connected_overlap[i].size()*2 > ta)continue;
        int j;
        for(j = 0; j < (int)connected_overlap[i].size(); j++) {
            int f = connected_overlap[i][j].first, s = connected_overlap[i][j].second;
            
        
            if(is_particle_collision(particles[f].loc, particles[f].vel, particles[s].loc, particles[s].vel, radius)) {
                no_collisions = false;
                resolve_particle_collision(particles[f].loc, particles[f].vel, particles[s].loc, particles[s].vel);
            }
           
            
        }
    }
    
    
    return no_collisions;
}

void detect_overlap2(std::vector<Particle> &particles, std::vector<std::vector<int>> &grid,  std::vector<int> &particle_map, int number_of_cells, int radius){
    overlap.clear();
    for(int i = 0; i < (int)particles.size(); i++) {

        int dirx[] = {0, 1,1,0,-1,-1,-1,0,1};
        int diry[] = {0, 0,number_of_cells,number_of_cells,number_of_cells,0,-number_of_cells,-number_of_cells,-number_of_cells};

        for(int d = 0; d < 9; d++) {
            int h = particle_map[i] +dirx[d] + diry[d];
            if(h < number_of_cells*number_of_cells && h >= 0 && (particle_map[i]%number_of_cells + dirx[d]) >= 0 && (particle_map[i]%number_of_cells + dirx[d]) < number_of_cells) {
                for(auto j : grid[h]) {
                    
                    if(j != i && is_particle_overlap(particles[i].loc, particles[j].loc, radius)) {
                        overlap.push_back(std::make_pair(i,j));
                    }
                }
            }
        }

    }
}

bool resolve_collisions2(std::vector<Particle> &particles, int square_size, int radius) {
    bool no_collisions = true; 
    square_size = square_size;
    int i;
    #pragma omp parallel for shared(particles) private(i) 
    for(i = 0; i < (int)particles.size(); i++) {
        if(is_wall_collision(particles[i].loc, particles[i].vel, square_size, radius)) {
            no_collisions = false;
            resolve_wall_collision(particles[i].loc, particles[i].vel, square_size, radius);
        }
    }

    #pragma omp parallel for shared(overlap, particles) private(i)
    for(i = 0; i < (int)overlap.size(); i++) {
        int f = overlap[i].first, s=overlap[i].second;
        omp_set_lock(&locks[std::min(f,s)]);
        omp_set_lock(&locks[std::max(f,s)]);
        if(is_particle_collision(particles[f].loc, particles[f].vel, particles[s].loc, particles[s].vel, radius)) {
            no_collisions = false;
            resolve_particle_collision(particles[f].loc, particles[f].vel, particles[s].loc, particles[s].vel);
        }
        omp_unset_lock(&locks[std::max(f,s)]);
        omp_unset_lock(&locks[std::min(f,s)]);
    }

    return no_collisions;
}


int main(int argc, char* argv[]) {
    auto start_start = std::chrono::high_resolution_clock::now();
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

    

    std::vector<std::vector<int>> overlaps;
    std::vector<int> particle_map(params.param_particles, 0);
    const int number_of_cells = std::min(1000, (int)((params.square_size)/(2*params.param_radius)));
    std::vector<std::vector<int>> grid(number_of_cells*number_of_cells, std::vector<int>(1,0));
    //std::cout<<number_of_cells<<std::endl;
    for(int i = 0; i < params.param_particles; i++) {
        std::vector<int> v;
        overlaps.push_back(v);
        omp_init_lock(&locks[i]);

    }
    auto start_end = std::chrono::high_resolution_clock::now();
    start_time += std::chrono::duration_cast<std::chrono::nanoseconds>(start_end - start_start).count();
    for(int i = 1; i <= params.param_steps; i++) {
        auto start = std::chrono::high_resolution_clock::now();


        

        auto assign_start = std::chrono::high_resolution_clock::now();
        update_positions(particles);

        assign_grid(particles, grid, particle_map, number_of_cells, params.square_size);
        auto assign_end = std::chrono::high_resolution_clock::now();
        assign_grid_time += std::chrono::duration_cast<std::chrono::nanoseconds>(assign_end - assign_start).count();

        
        auto overlap_start = std::chrono::high_resolution_clock::now();
        detect_overlaps(particles, grid, overlaps, particle_map, number_of_cells, params.param_radius);
        //detect_overlap2(particles,grid, particle_map, number_of_cells,params.param_radius);
        auto overlap_end = std::chrono::high_resolution_clock::now();
        overlap_time += std::chrono::duration_cast<std::chrono::nanoseconds>(overlap_end - overlap_start).count();
        ma = ta = 0;
        


        auto bfs_start = std::chrono::high_resolution_clock::now();

        for(int j = 0; j<params.param_particles; j++)visit[j]=false;
        int index = 0;
        
        for(int j = 0; j < params.param_particles; j++) {
            if(visit[j])continue;
            bfs(particles, overlaps, j, index, params.square_size, params.param_radius);
            index++;
        }
        auto bfs_end = std::chrono::high_resolution_clock::now();
        bfs_time += std::chrono::duration_cast<std::chrono::nanoseconds>(bfs_end - bfs_start).count();
        mi = std::min(mi, index);
        
        
        

        auto collision_start = std::chrono::high_resolution_clock::now();
        while(!resolve_collisions(particles, index, params.square_size, params.param_radius)) {
        }
        auto collision_end = std::chrono::high_resolution_clock::now();
        collision_time += std::chrono::duration_cast<std::chrono::nanoseconds>(collision_end - collision_start).count();


        auto end = std::chrono::high_resolution_clock::now();
        long long t = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        total_time += t;

        std::cout<<i<<" "<<index<<" "<<t/1000000<<" "<<ma<<" "<<ta<<std::endl;

        
        #if CHECK == 1
            // Check final positions
            validator.validate_step(particles);
        #endif
        
    }
    std::cout<<total_time/1000000<<" "<<assign_grid_time/1000000<<" "<<overlap_time/1000000<<" "<<bfs_time/1000000<<" "<<collision_time/1000000<< " " <<start_time/1000000<<" "<<mi<<std::endl;

}