#include<physics_engine_uniform.h>
#include<thread>

PhysEngine::PhysEngine(float width, float height, int nmesh_hrz, int nmesh_vrt, int nparticles,glm::vec2 center,int threads):
// width(width), height(height), nmesh_hrz(nmesh_hrz), nmesh_vrt(nmesh_vrt), nparticles(nparticles), mesh(particles, nmesh_hrz, nmesh_vrt, width, height, center)
width(width), height(height), nmesh_hrz(nmesh_hrz), nmesh_vrt(nmesh_vrt), nparticles(nparticles),threads(threads)
{
    
    this->particles.resize(nparticles, nullptr); 
    this->collision_ind.resize(nparticles, -1); 
    

    // // initialize bodies
    // for (int i=0; i<nparticles; i++){
    //     glm::vec2 pos(rand()%1000/500.0f-1.f,rand()%1000/500.0f-1.f);
    //     // glm::vec2 vel(rand()%1000/500.0f-1.f,rand()%1000/500.0f-1.f);
    //     // glm::vec2 vel = -pos;
    //     glm::vec2 vel(0,0);
    //     // vel *= 1e3;
    //     this->particles[i] = new Particle(i,pos,vel,glm::vec2(0,0),10.,0.,0.);

    //     if (i==0){
    //         this->particles[i] = new Particle(i,glm::vec2(0,0),glm::vec2(0.014142,0),glm::vec2(0,0),1.,0.,0.);
    //     }
    //     else if (i==1){
    //         this->particles[i] = new Particle(i,glm::vec2(1.,0),glm::vec2(0,1.),glm::vec2(-1.,0),1e-8,0.,0.);
    //     }
    //     else if (i==2){
    //         this->particles[i] = new Particle(i,glm::vec2(-1.,0),glm::vec2(0,-1.),glm::vec2(-1.,0),-1e-8,0.,0.);
    //     }
    //     else if (i==3){
    //         this->particles[i] = new Particle(i,glm::vec2(0,0.5),glm::vec2(-1.4142,0),glm::vec2(0,-4),1e-2,0.,0.);
    //     }
    //     else if (i==4){
    //         // std::cout << "init" << std::endl;
    //         this->particles[i] = new Particle(i,glm::vec2(0.01,0.5),glm::vec2(-1.4142,1.),glm::vec2(-10.,-4),1e-16,0.,0.);
    //     }
    //     else if (i==5){
    //         // std::cout << "init" << std::endl;
    //         this->particles[i] = new Particle(i,glm::vec2(0,0.4292625434084361),glm::vec2(-1.4142*0.4292625434084361/0.5,0),glm::vec2(0.,-2.357),1e-16,0.,0.);
    //     }
    //     else if (i==6){
    //         // std::cout << "init" << std::endl;
    //         this->particles[i] = new Particle(i,glm::vec2(0.25,0.4330127018922193),glm::vec2(-1.23,0.71),glm::vec2(-2,-3.4641),1e-16,0.,0.);
    //     }

    //     // if (i==0){
    //     //     this->particles[i] = new Particle(i,glm::vec2(-1,-1),glm::vec2(1,1),glm::vec2(0,0),1e-2,0.,0.);
    //     // }
    //     // else if (i==1){
    //     //     this->particles[i] = new Particle(i,glm::vec2(1.,1),glm::vec2(-1,-1),glm::vec2(0,0),1e-2,0.,0.);
    //     // }
    //     this -> particles[i] -> updateFromBuffer();
    // }


    // // Spiral parameters
    // float spiral_spacing = 0.01f;      // Controls how far apart each turn of the spiral is
    // float angular_increment = 0.2f;  // Controls the angle between consecutive particles in radians
    // float base_mass = 10.0f;         // Default mass for particles
    // float base_energy = 0.0f;        // Default energy
    // float base_ang_momentum = 1000.0f;  // Default angular momentum

    // for (int i = 0; i < nparticles; i++) {
    //     float angle = i * angular_increment;
    //     float radius = i * spiral_spacing;

    //     // Position in Cartesian coordinates
    //     glm::vec2 pos = radius * glm::vec2(cos(angle), sin(angle));

    //     // Velocity tangent to the spiral (optional: scale it based on radius)
    //     glm::vec2 vel = radius * angular_increment * glm::vec2(-sin(angle), cos(angle));

    //     // Zero initial acceleration
    //     glm::vec2 acc(0.0f, 0.0f);

    //     // Initialize the particle
    //     this->particles[i] = new Particle(i, pos, vel, acc, base_mass, base_energy, base_ang_momentum);
    // }

    float center_mass = 1e5;

    // Disk parameters
    float radius_max = 3.f;          // Maximum radius of the disk
    float radius_min = 0.01f;
    float base_mass = 1e-1f;          // Default mass for particles
    // float disk_ang_mom = 3.f;       // Angular momentum scaling factor
    float disk_ang_mom = sqrt(center_mass);       // Angular momentum scaling factor
    float base_energy = 0.0f;         // Default energy
    float base_ang_momentum = 0.0f;   // Default angular momentum

    for (int i = 0; i < nparticles; i++) {
        // Generate random radius and angle within the disk
        float radius = radius_max * std::pow((static_cast<float>(rand()) / RAND_MAX),0.5)+radius_min; // Uniform distribution in area
        float angle = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;       // Random angle in [0, 2π]

        // Position in Cartesian coordinates
        glm::vec2 pos = radius * glm::vec2(cos(angle), sin(angle));

        // Random tangential velocity magnitude
        float random_factor = 0.5f + static_cast<float>(rand()) / RAND_MAX; // Random scaling factor in [0.5, 1.5]
        // float velocity_magnitude = radius * disk_ang_mom * random_factor;  // Tangential velocity scaled by radius and randomness
        float velocity_magnitude = disk_ang_mom * random_factor/radius;  // Tangential velocity scaled by radius and randomness

        // Tangential velocity (perpendicular to radial vector)
        glm::vec2 vel = velocity_magnitude * glm::vec2(-sin(angle), cos(angle));

        // Zero initial acceleration
        glm::vec2 acc(0.0f, 0.0f);

        float mass_random_factor = 0.5f + static_cast<float>(rand()) / RAND_MAX; // Random mass scaling factor in [0.5, 1.5]
        float mass = base_mass * std::pow(2,mass_random_factor);

        // Initialize the particle
        this->particles[i] = new Particle(i, pos, vel, acc, mass, base_energy, base_ang_momentum);
        this->particles[i]->calcRenderParams();
    }
    
    this->particles[0] = new Particle(0, glm::vec2(0, 0), glm::vec2(0, 0), glm::vec2(0, 0), center_mass, 10000.0f, 0.0f);

    // this -> mesh = new Mesh(this->particles,nmesh_hrz,nmesh_vrt,width,height,center);

    // initializes Barnes hut

    std::vector<int> active_ind;
    active_ind.resize(particles.size(),0);
    for (int i=0;i<int(particles.size());i++){
        active_ind[i] = i;
    }
    this -> mesh = new Mesh(this->particles,nmesh_hrz,nmesh_vrt,width,height,glm::vec2(0,0));
    // this -> mesh -> particles = particles;

    this -> mesh -> update_mesh();



    // std::cout << this->mesh->sub_mesh.size() << std::endl;
    this -> dt = 5e-3/(this -> mesh -> tot_mass)*width/3;
    this -> dt = 2e-3;
    // this -> dt = 1e-5;
    this -> dt = 5e-5;
    std::cout << dt << std::endl;

}

void PhysEngine::update(){
    // update mesh
    // std::cout << "update mesh" << std::endl;
    this -> mesh -> update_mesh();

    // for (Particle* particle: this->particles){
    //     if (particle!=nullptr){
    //         if (particle->mesh->particles.size() != 10000){
    //             std::cout << "SDFIEG" << std::endl;
    //         }
    //     }
    // }


    // for (Particle* particle: this->particles){
    //     if (particle->mesh->particles.size() == 0){
    //         std::cout << "wtf" << std::endl;
    //     }
    // }
    // update state
    // std::cout << "kdk" << std::endl;
    // this -> kick_drift_kick();

    // for (int i=0; i<this->nparticles; i++){
    //     if (particles[i] != nullptr){

    //         this->kick_drift_kick(this->particles[i]);
    //     }
    // }

    // std::vector<std::future<void>> futures;
    // futures.reserve(this->nparticles);

    // for (int i = 0; i < this->nparticles; i++) {
    //     if (particles[i] != nullptr) {
    // // for (int i: this->mesh->active_ind) {
    // //     if (particles[i] != nullptr && mesh->particles[i] != nullptr) {
    //         // Launch asynchronous tasks instead of threads
    //         futures.push_back(std::async(std::launch::async, &GravEngine::kick_drift_kick, this, particles[i]));
    //     }
    // }

    // // Wait for all tasks to finish
    // for (auto& future : futures) {
    //     future.get();
    // }

    // Create a thread pool with a size of, for example, 4 threads (adjust as needed)
    ThreadPool pool(this->threads);

    // this->pool->restart(this->threads);
    // Enqueue tasks for each active particle
    for (int i=0; i<this->nparticles; i++){
    // for (int i : this->mesh->active_ind) {
        if (particles[i] != nullptr && mesh->particles[i] != nullptr) {
            // Enqueue the task to the thread pool
            // this->pool->enqueue([this, i]() {
            pool.enqueue([this, i]() {
                // Call the original method in the task
                this->kick_drift_kick(particles[i]);
            });
        }
    }

    // No need to explicitly wait for futures, as the thread pool will manage the synchronization
    // The thread pool will block until all tasks have been completed when it is stopped
    // pool->stopPool();
    pool.stopPool();

    // // int sum_part = 0;
    // for (int i=0; i<this->nparticles; i++){
    //     this->postUpdateProcess(i);
    // }

    // std::vector<std::future<void>> futures1;
    // futures.reserve(this->nparticles);

    // for (int i = 0; i < this->nparticles; i++) {
    //     if (particles[i] != nullptr) {
    //         // Launch asynchronous tasks instead of threads
    //         futures1.push_back(std::async(std::launch::async, &GravEngine::postUpdateProcess, this, i));
    //     }
    // }

    // // Wait for all tasks to finish
    // for (auto& future : futures1) {
    //     future.get();
    // }

    pool.restart(this->threads);
    for (int i = 0; i < this->nparticles; i++) {
        if (particles[i] != nullptr) {
            // Launch asynchronous tasks instead of threads
            // this->pool->enqueue([this, i]() {
            pool.enqueue([this, i]() {
                // Call the original method in the task
                this->postUpdateProcess(i);
            });        
        }
    }
    // this->pool->stopPool();
    pool.stopPool();


    // int sum=0;
    // for (NbodyMesh* mesh: this->mesh->sub_mesh){
    //     for (Particle* particle: mesh->particles){
    //         if (particle!=nullptr){
    //             sum+=1;
    //         }
    //     }
    // }
    // // std::cout << sum << std::endl;

    // if (sum!=sum_part){
    //     std::cout << "FUCKSLDJ" << std::endl;

    // }
}

void PhysEngine::postUpdateProcess(int i){
    if (particles[i] != nullptr){
    // sum_part += 1;
    // this -> particles[i] -> calcRenderParams();

    // this -> particles[i] -> updateFromBuffer();
    
    if (this->collision_ind[i]!=-1){
        // collision handler
        // std::cout << "collision" << std::endl;
        // ensure second particle is update wrt buffer
        // this -> particles[this->collision_ind[i]] -> updateFromBuffer();

        // limit to 2 way collision
        if (this->particles[i] != nullptr && this->particles[this->collision_ind[i]]!=nullptr){
            // this -> collision_handler(this->particles[i],this->particles[this->collision_ind[i]]);
            this -> collision_ind[i] = -1;
        }
        // this -> particles[i] -> calcRenderParams();
    }
    if (particles[i] != nullptr){
        this -> particles[i] -> updateFromBuffer();
        // bc handler
        // std::cout << "bc" << std::endl;
        this -> bc_handler(this->particles[i]);
    }
}
}


void PhysEngine::kick_drift_kick(Particle* particle){
    if (particle != nullptr){
        particle->buffer_acc = this->calc_eom(particle);

        // for (Particle* particle1: this->particles){
        //     if (particle1->mesh->particles.size() == 0){
        //         std::cout << "wtf" << std::endl;
        //     }
        // }

        // std::cout << "calc" << std::endl;
        // std::cout << i << ": " << particle->buffer_pos.x << std::endl;
        // std::cout << particle->buffer_acc.x  << ", " << particle->buffer_acc.y << std::endl;
        // particle->buffer_vel += 0.5f*particle->buffer_acc*dt;

        // particle->buffer_vel += 0.5f*particle->buffer_acc*dt;

        particle->buffer_vel += particle->buffer_acc*dt;

        // vel at 0.5 step forward
        particle->buffer_pos += particle->buffer_vel*dt;
        // acc = this->calc_eom();
        // particle.acc = acc;
        particle->calcRenderParams();
    }
}


PhysEngine::~PhysEngine() {
    for (Particle* particle : this->particles){
        delete particle;
    }
}

// define gravity engine specific functions
GravEngine::GravEngine(float width, float height, int nmesh_hrz, int nmesh_vrt, int nparticles,glm::vec2 center,int threads):
PhysEngine(width,height,nmesh_hrz,nmesh_vrt,nparticles,center,threads)
{}



glm::vec2 GravEngine::calc_eom(Particle* particle){
    // calculate acceleration for each particle
    float force = 0.;
    float dist_squared;
    glm::vec2 acc(0,0);
    Mesh* mesh;
    
    if (particle->mesh != nullptr){
        for (Particle* body_j:particle->mesh->particles){
            // std::cout << "here" << std::endl;
            if (body_j!=nullptr){
                if (body_j->ind != particle->ind){
                    glm::vec2 rel_vec = body_j->pos - particle->pos;

                    dist_squared = dot(rel_vec,rel_vec);
                    if (dist_squared == 0){
                        // std::cout << std::endl <<std::endl << std::endl <<  dist_squared << std::endl;
                        dist_squared = 1e-6;
                    }
                    this -> checkCollision(dist_squared,particle->ind,body_j->ind);


                    acc += this->G*(body_j->mass)/(dist_squared*sqrt(dist_squared))*rel_vec;
                    // std::cout << dot(acc,acc) << std::endl;
                }
            }
        }

        // get higest mesh

        this->mesh_tree_calc(&acc,particle,this->mesh);
    }

    // for other meshes calc theta (width/dist)
    // move from highest to lowest. if in highest check if in submesh. if not in submesh cacl width/dist. if above threshold, calc indiv contributions

    return acc;

    // return glm::vec2(0,0);
}


// void GravEngine::mesh_tree_calc(glm::vec2* acc, Particle* particle, NbodyMesh* mesh){
void GravEngine::mesh_tree_calc(glm::vec2* acc, Particle* particle, NbodyMesh* mesh){
    // base_mesh->parent_mesh;
    // get highest mesh
    if (mesh->particles[particle->ind] != nullptr){
        // if in this mesh, find lower sub meshes
        if (!mesh->sub_mesh.empty()){
            // for (NbodyMesh* sub_mesh: mesh->sub_mesh){
            for (NbodyMesh* sub_mesh: mesh->sub_mesh){
                this->mesh_tree_calc(acc,particle,sub_mesh);
            }
        }
    }
    else{
        bool mesh_check=false;
        // calculate width/dist
        if (mesh->tot_mass!=0){
            // glm::vec2 rel_vec = mesh->com - particle->buffer_pos;
            glm::vec2 rel_vec = mesh->center - particle->pos;
            // if (mesh->char_len/dist<1){

            if (mesh->char_len/(std::abs(rel_vec.x)+std::abs(rel_vec.y))<1.2){
                float dist = dot(rel_vec,rel_vec);
                if (mesh->char_len*mesh->char_len/dist<1.2){
                    glm::vec2 rel_vec = mesh->com - particle->pos;
                    dist = sqrt(dot(rel_vec,rel_vec));
                    // if within threshold use com and tot mass of mesh
                    *acc += this->G * (mesh->tot_mass)/(dist*dist*dist)*rel_vec;
                    mesh_check = true;
                // std::cout << "mesh center calc" << std::endl;
                }
            }
            if (!mesh_check) {
                // continue down meshes
                if (!mesh->sub_mesh.empty()){
                    // for (NbodyMesh* sub_mesh: mesh->sub_mesh){
                    for (NbodyMesh* sub_mesh: mesh->sub_mesh){
                        this->mesh_tree_calc(acc,particle,sub_mesh);
                    }
                }
                else {
                    float dist_squared;
                    float dist;

                    // if already at lowest, calculate individual particle contributions
                    for (Particle* body_j:mesh->particles){
                        if (body_j != nullptr){
                            // glm::vec2 rel_vec = body_j->buffer_pos - particle->buffer_pos;
                            glm::vec2 rel_vec = body_j->pos - particle->pos;
                            dist_squared = dot(rel_vec,rel_vec);
                            dist = sqrt(dist_squared);

                            if (dist_squared == 0){
                                // std::cout << std::endl <<std::endl << std::endl <<  dist_squared << std::endl;
                                dist_squared = 1e-6;
                            }
                            // this -> checkCollision(dist,particle->ind,body_j->ind);

                            *acc += this->G*(body_j->mass)/(dist_squared*dist)*rel_vec;
                            // std::cout << "body calc" << std::endl;
                        }
                    }
                }
            }
        }
    }
}


void GravEngine::checkCollision(float dist,int particle1_ind,int particle2_ind){
    float dist_threshold = std::max(this->particles[particle1_ind]->mass,this->particles[particle2_ind]->mass);
    // dist_threshold = 0.00005*cbrt(dist_threshold);
    dist_threshold = 0.0002*cbrt(dist_threshold);
    // dist_threshold = 0.02*cbrt(dist_threshold); 
    // std::cout << dist_threshold << ", " << dist << std::endl;
    if (dist < dist_threshold){
        // this->collision_handler(particle,body_j);
        // continue;
        if (collision_ind[particle2_ind]!=particle1_ind){
            this->collision_ind[particle1_ind] = particle2_ind;
        }
    }
}

void GravEngine::collision_handler(Particle* particle1, Particle* particle2) {

    // std::cout<<"collision"<<std::endl;
    // std::cout << particle1->buffer_vel.x << ", " << particle1->buffer_vel.y << std::endl;
    // std::cout << particle2->buffer_vel.x << ", " << particle2->buffer_vel.y << std::endl;
    // Check if both particles are valid
    if (!particle1 || !particle2) return;

    // Calculate new combined properties
    float new_mass = particle1->mass + particle2->mass;
    glm::vec2 new_pos = (particle1->pos * particle1->mass + particle2->pos * particle2->mass) / new_mass;
    glm::vec2 new_vel = (particle1->vel * particle1->mass + particle2->vel * particle2->mass) / new_mass;

    // Calculate internal energy
    float kinetic_energy_1 = 0.5f * particle1->mass * dot(particle1->vel, particle1->vel);
    float kinetic_energy_2 = 0.5f * particle2->mass * dot(particle2->vel, particle2->vel);
    float combined_kinetic_energy = 0.5f * new_mass * dot(new_vel, new_vel);
    float internal_energy_add = kinetic_energy_1 + kinetic_energy_2 - combined_kinetic_energy;

    // Update particle1
    particle1->mass = new_mass;
    particle1->buffer_pos = new_pos;
    particle1->buffer_vel = new_vel;
    particle1->internal_energy += internal_energy_add + particle2->internal_energy;

    // Mark particle2 for deletion
    particle2->mesh->remove_particle(particle2->ind);
    this -> mesh -> particles[particle2->ind] = nullptr;
    this -> particles[particle2->ind] = nullptr;

    // particle2->active = false;  // Use a flag instead of deleting immediately
}

void GravEngine::bc_handler(Particle* particle){
    // apply periodic bc
    if (particle->buffer_pos.x<this->mesh->corners[0].x*0.95){
        particle->buffer_pos.x = this->mesh->corners[2].x*0.9;
        particle->buffer_vel.x = 0;
    }
    else if(particle->buffer_pos.x>this->mesh->corners[2].x*0.95){
        particle->buffer_pos.x = this->mesh->corners[0].x*0.9;
        particle->buffer_vel.x = 0;
    }

    if (particle->buffer_pos.y<this->mesh->corners[0].y*0.95){
        particle->buffer_pos.y = this->mesh->corners[2].y*0.9;
        particle->buffer_vel.y = 0;
    }
    else if(particle->buffer_pos.y>this->mesh->corners[2].y*0.95){
        particle->buffer_pos.y = this->mesh->corners[0].y*0.9;
        particle->buffer_vel.y = 0;
    }

}
