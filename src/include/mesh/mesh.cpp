#include <mesh.h>

BarnesHutMesh::BarnesHutMesh(std::vector<Particle*> particles,float width, float height, glm::vec2 center, std::vector<int> active_ind, int active)
: NbodyMesh(particles, width, height, center, nullptr),active_ind(active_ind),npart_active(active)  // Call NbodyMesh constructor
{
}

void BarnesHutMesh::modify_mesh(){
    // this->deleteTree();
    // remove references in particles to old mesh
    for (Particle* particle: this->particles){
        if (particle!=nullptr){
            particle->mesh = nullptr;
        }
    }

    // create new mesh
    // ThreadPool pool(8);
    this->checkActive();
    this->splitMesh();
}

void BarnesHutMesh::update_mesh(){
    // update which mesh particles are in
    // generates or modifies mesh
    this -> modify_mesh();
}

void BarnesHutMesh::deleteTree(){
    if (!this->sub_mesh.empty()){
        for (BarnesHutMesh* mesh: this->sub_mesh){
            mesh->deleteTree();
            delete mesh;
        }
    }
    this->sub_mesh.clear();  // Clear the vector after deletion
}


void BarnesHutMesh::reinitialize(std::vector<Particle*> particles, float width, float height, glm::vec2 center, std::vector<int> active_ind, int active) {
    this->width = width;
    this->height = height;
    this->center = center;
    this->active_ind = active_ind;
    this->tot_mass = 0;
    this->com = glm::vec2(0, 0);
    this->particles = particles;
    this->npart_active = active;
}

void BarnesHutMesh::splitMesh(){
    this -> tot_mass = 0;
    this -> com = glm::vec2(0,0);
    // if (this->active>100){
    if (this->npart_active>10){
        if (this->sub_mesh.empty()) {
            this->sub_mesh.resize(4, nullptr);  // Ensure submesh is resized only once
        }        
        
        // split mesh
        float submesh_width = this->width/2;
        float submesh_height = this->height/2;    
        
        // ----------------- PARALLEL SEEMS SLOWER ----------------
        // const int num_threads = 4;
        // const int n = this->active_ind.size(); // Total number of active particles
        // const int chunk_size = (n + num_threads - 1) / num_threads; // Chunk size for each thread

        // std::vector<std::thread> threads;  // Store threads
        // threads.reserve(32);

        // // reassign active_ind to quadrants using this->center
        // // code stuff
        // std::vector<int> quadrants[4];
        // std::vector<int> counters(4, 0);  // Array to track the number of particles per quadrant
        // std::vector<Particle*> particles_sub[4];  // Array to track the number of particles per quadrant

        // // First pass: Reserve space and populate the quadrants
        // for (int i = 0; i < 4; ++i) {
        //     quadrants[i].resize(this->particles.size(),-1);  // Reserve space for each quadrant
        //     particles_sub[i].resize(this->particles.size(),nullptr);
        // }

        // // second pass: sort
        // for (int i = 0; i < this->active_ind.size(); ++i) {
        //     // if (this->particles[ind]!=nullptr){
        //         int ind = this->active_ind[i];
        //         glm::vec2 pos = this->particles[ind]->pos;

        //         if (pos.y > this->center.y) {
        //             if (pos.x < this->center.x) {
        //                 quadrants[3][counters[3]] = ind;  // Upper-left
        //                 counters[3]++;
        //                 particles_sub[3][ind] = this->particles[ind];
        //             }
        //             else {
        //                 quadrants[2][counters[2]] = ind;  // Upper-right
        //                 counters[2]++;
        //                 particles_sub[2][ind] = this->particles[ind];
        //             }
        //         } else {
        //             if (pos.x > this->center.x) {
        //                 quadrants[1][counters[1]] = ind; // Bottom-right
        //                 counters[1]++;
        //                 particles_sub[1][ind] = this->particles[ind];
        //             }
        //             else {
        //                 quadrants[0][counters[0]] = ind;  // Bottom-left
        //                 counters[0]++;
        //                 particles_sub[0][ind] = this->particles[ind];
        //             }
        //         }
        //     // }
        // }

        // // // remove excess
        // // for (int i = 0; i < 4; ++i) {
        // //     quadrants[i].resize(counters[i]);  // Resize to remove excess space
        // // }
        // // -----------------------


        // // Lambda function for sorting a range of particles
        // auto sortChunk = [&](int start, int end) {
        //     for (int i = start; i < end; ++i) {
        //         int ind = this->active_ind[i];
        //         glm::vec2 pos = this->particles[ind]->pos;

        //         if (pos.y > this->center.y) {
        //             if (pos.x > this->center.x) {
        //                 quadrants[3][counters[3]++] = ind; // Upper-left
        //                 particles_sub[3][ind] = this->particles[ind];
        //             } else {
        //                 quadrants[2][counters[2]++] = ind; // Upper-right
        //                 particles_sub[2][ind] = this->particles[ind];
        //             }
        //         } else {
        //             if (pos.x > this->center.x) {
        //                 quadrants[1][counters[1]++] = ind; // Bottom-right
        //                 particles_sub[1][ind] = this->particles[ind];
        //             } else {
        //                 quadrants[0][counters[0]++] = ind; // Bottom-left
        //                 particles_sub[0][ind] = this->particles[ind];
        //             }
        //         }
        //     }
        // };

        // // Create and launch threads
        // // std::vector<std::thread> threads;
        // for (int t = 0; t < num_threads; ++t) {
        //     int start = t * chunk_size;
        //     int end = std::min(start + chunk_size, n); // Make sure end does not exceed n
        //     threads.emplace_back(sortChunk, start, end);
        // }

        // // Join threads
        // for (std::thread& t : threads) {
        //     if (t.joinable()) {
        //         t.join();
        //     }
        // }


        // // remove excess
        // for (int i = 0; i < 4; ++i) {
        //     quadrants[i].resize(counters[i]);  // Resize to remove excess space
        // }


        // for (int i=0; i<4; i++){
        //     glm::vec2 submesh_center = 0.5f*(this->center+this->corners[i]);

        //     // this->sub_mesh[i] = new BarnesHutMesh(this->particles,submesh_width,submesh_height,submesh_center,this->active_ind);
        //     // this->sub_mesh[i] -> checkActive();
        //     // this->sub_mesh[i] = new BarnesHutMesh(particles_sub[i],submesh_width,submesh_height,submesh_center,quadrants[i]);
        //     if (!this->sub_mesh[i]) {
        //         this->sub_mesh[i] = new BarnesHutMesh(particles_sub[i],submesh_width,submesh_height,submesh_center,quadrants[i],counters[i]);
        //     } else {
        //         this->sub_mesh[i]->reinitialize(particles_sub[i],submesh_width,submesh_height,submesh_center,quadrants[i],counters[i]);
        //     }


        //     // continue to check if need to split
        //     this->sub_mesh[i]->splitMesh();


        //     // calc com and tot mass
        //     if (this->sub_mesh[i]->tot_mass != 0){
        //         this->tot_mass += this->sub_mesh[i]->tot_mass;
        //         this->com += this->sub_mesh[i]->tot_mass*this->sub_mesh[i]->com;
        //     }
        // }


        std::vector<std::thread> threads;  // Store threads

        for (int i = 0; i < 4; i++) {
            glm::vec2 submesh_center = 0.5f * (this->center + this->corners[i]);
            // glm::vec2 submesh_center = 0.5f * (this->com + this->corners[i]);

            // threads.emplace_back([this, i, submesh_width, submesh_height, submesh_center,particles_sub,quadrants,counters]() {
            //     if (!this->sub_mesh[i]) {
            //         this->sub_mesh[i] = new BarnesHutMesh(particles_sub[i],submesh_width,submesh_height,submesh_center,quadrants[i],counters[i]);
            //     } else {
            //         this->sub_mesh[i]->reinitialize(particles_sub[i],submesh_width,submesh_height,submesh_center,quadrants[i],counters[i]);
            //     }

            threads.emplace_back([this, i, submesh_width, submesh_height, submesh_center]() {
                if (!this->sub_mesh[i]) {
                    this->sub_mesh[i] = new BarnesHutMesh(this->particles,submesh_width,submesh_height,submesh_center,this->active_ind,this->npart_active);
                } else {
                    this->sub_mesh[i]->reinitialize(this->particles,submesh_width,submesh_height,submesh_center,this->active_ind,this->npart_active);
                }

                // this->sub_mesh[i] = new BarnesHutMesh(this->particles,submesh_width,submesh_height,submesh_center,this->active_ind,this->npart_active);

                this->sub_mesh[i]->checkActive();


                // if (counters[i]!=this->sub_mesh[i]->npart_active){
                //     std::cout<<"Sdf"<<std::endl;
                // }

                // Continue to check if need to split
                this->sub_mesh[i]->splitMesh();

                // Update mass and center of mass
                if (this->sub_mesh[i]->tot_mass != 0) {
                    // std::lock_guard<std::mutex> lock(this->mutex);
                    this->tot_mass += this->sub_mesh[i]->tot_mass;
                    this->com += this->sub_mesh[i]->tot_mass * this->sub_mesh[i]->com;
                }
            });
        }

        // Wait for all threads to complete
        for (auto& t : threads) {
            if (t.joinable()) {
                t.join();
            }
        }


        // std::vector<std::future<void>> futures;  // Store futures for async tasks

        // for (int i = 0; i < 4; i++) {
        //     glm::vec2 submesh_center = 0.5f * (this->center + this->corners[i]);

        //     // futures.push_back(std::async(std::launch::async, [this, i, submesh_width, submesh_height, submesh_center]() {
        //     futures.push_back(std::async(std::launch::async, [this, i, submesh_width, submesh_height, submesh_center]() {
        //         if (!this->sub_mesh[i]) {
        //             this->sub_mesh[i] = new BarnesHutMesh(this->particles, submesh_width, submesh_height, submesh_center, this->active_ind);
        //         } else {
        //             this->sub_mesh[i]->reinitialize(this->particles,submesh_width, submesh_height, submesh_center, this->active_ind);
        //             // this->sub_mesh[i]->particles = this->particles;
        //         }
        //         // this->sub_mesh[i] = new BarnesHutMesh(this->particles, submesh_width, submesh_height, submesh_center, this->active_ind);

        //         // Continue to check if need to split
        //         this->sub_mesh[i]->splitMesh();

        //         // Update mass and center of mass
        //         if (this->sub_mesh[i]->tot_mass != 0) {
        //             // std::lock_guard<std::mutex> lock(this->mutex);
        //             this->tot_mass += this->sub_mesh[i]->tot_mass;
        //             this->com += this->sub_mesh[i]->tot_mass * this->sub_mesh[i]->com;
        //         }
        //     }));
        // }


        // // Wait for all tasks to complete
        // for (auto& fut : futures) {
        //     fut.get();  // Block until the future completes
        // }


        // Create a thread pool with 4 threads (or adjust as needed)
        // if (*pool == nullptr){
        //     ThreadPool pool(8);
        // }
        // std::mutex mutex;  // Mutex for updating shared variables safely (tot_mass and com)

        // // Start tasks for each submesh
        // for (int i = 0; i < 4; ++i) {
        //     glm::vec2 submesh_center = 0.5f * (this->center + this->corners[i]);

        //     pool.enqueue([this, i, submesh_width, submesh_height, submesh_center, &pool]() {
        //         // Create the submesh
        //         this->sub_mesh[i] = new BarnesHutMesh(this->particles, submesh_width, submesh_height, submesh_center, this->active_ind);

        //         // Continue to check if the mesh needs further splitting
        //         this->sub_mesh[i]->splitMesh(pool);

        //         // Update the total mass and center of mass for the parent mesh
        //         if (this->sub_mesh[i]->tot_mass != 0) {
        //             // std::lock_guard<std::mutex> lock(mutex);  // Lock for thread-safe updates
        //             this->tot_mass += this->sub_mesh[i]->tot_mass;
        //             this->com += this->sub_mesh[i]->tot_mass * this->sub_mesh[i]->com;
        //         }
        //     });
        // }

        // // Stop the thread pool and wait for all tasks to finish
        // pool.stopPool();

        
        if (this->tot_mass!=0){
            this->com /= tot_mass;
        }

    }
    else {
        // assign values

        if (npart_active!=0){

            if (!this->sub_mesh.empty()){
                for (BarnesHutMesh* mesh: this->sub_mesh){
                    mesh->deleteTree();
                    delete mesh;
                }
                this->sub_mesh.clear();
            }
            
            this -> com = glm::vec2(0,0);
            this -> tot_mass = 0;
            for (int ind: this->active_ind){
                this -> com += this->particles[ind]->pos * this->particles[ind]->mass;
                this -> tot_mass += this->particles[ind]->mass;
                this->particles[ind]->mesh = this;
            }
            this -> com /= tot_mass;
        }
    }
    // if (this->particles.size() != 10000){
    //     std::cout << "FDK" <<std::endl;
    // }
    
}

void BarnesHutMesh::return_corners(std::vector<std::vector<glm::vec2>>* vertices,std::vector<glm::vec2>* com){
    vertices->push_back(this->corners);
    com->push_back(this->com);
    if (!this->sub_mesh.empty()){
        for (BarnesHutMesh* mesh: this->sub_mesh){
            mesh->return_corners(vertices,com);
        }
    }
}

void BarnesHutMesh::checkActive(){

    this->npart_active = 0;

    auto isInsideBoundary = [this](int ind) {
        return this->particles[ind] != nullptr &&
               this->particles[ind]->pos.x > this->corners[0].x &&
               this->particles[ind]->pos.x <= this->corners[2].x &&
               this->particles[ind]->pos.y > this->corners[0].y &&
               this->particles[ind]->pos.y <= this->corners[2].y;
    };

    // Remove indices for particles that are outside the boundary
    std::vector<int>::iterator it = std::remove_if(
        this->active_ind.begin(), this->active_ind.end(),
        [this, &isInsideBoundary](int ind) -> bool {
            if (!isInsideBoundary(ind)) {
                // Set particles outside the boundary to nullptr
                this->particles[ind] = nullptr;
                return true; // Mark for removal
            }
            return false; // Keep if inside the boundary
        });

    // Erase the removed elements
    this->active_ind.erase(it, this->active_ind.end());

    // Count remaining active particles
    this->npart_active = int(this->active_ind.size());
}

void BarnesHutMesh::remove_particle(int ind){
    // removes particle from mesh starting at root
    this->particles[ind] = nullptr;
    if (this->parent_mesh != nullptr){
        this->parent_mesh->remove_particle(ind);
    }
}

Mesh::Mesh(std::vector<Particle*> particles, int nmesh_hrz, int nmesh_vrt, float width, float height, glm::vec2 center)
: NbodyMesh(particles, width, height, center, nullptr)  // Call NbodyMesh constructor
{
    // this->particles.resize(particles.size(),nullptr);

    // this->width = width;
    // this->height = height;

    // Declare mesh as highest tree
    // this->parent_mesh = nullptr;

    this->sub_mesh.resize(nmesh_hrz * nmesh_vrt, nullptr);  // Properly initialize sub_mesh
    float dx = width / nmesh_hrz;
    float dy = height / nmesh_vrt;

    glm::vec2 start_pos(center.x-width/2,center.y-height/2);

    // Create sub-meshes with proper initialization
    for (int i = 0; i < nmesh_hrz; i++) {
        for (int j = 0; j < nmesh_vrt; j++) {
            // Calculate the center of each sub-mesh
            glm::vec2 center(i * dx + dx / 2+start_pos.x, j * dy + dy / 2+start_pos.y);
            // Create sub-meshes and assign them
            this->sub_mesh[i * nmesh_vrt + j] = new NbodyMesh(std::vector<Particle*>(this->particles.size(),nullptr), dx, dy, center, this);
        }
    }
}



NbodyMesh::NbodyMesh(std::vector<Particle*> particles,float width, float height, glm::vec2 center, NbodyMesh* parent_mesh):
particles(particles),center(center),com(glm::vec2(0.,0.)),tot_mass(0.),tot_energy(0.),tot_ang_momentum(0.),tot_momentum(0.),parent_mesh(parent_mesh),
width(width),height(height)
// center(center),com(glm::vec2(0.,0.)),tot_mass(0.),tot_energy(0.),tot_ang_momentum(0.),tot_momentum(0.),parent_mesh(parent_mesh)
{
    // std::cout << "initialized" << std::endl;
    // this->particles.resize(particles.size(),nullptr);

    // std::cout << this->particles[0] << std::endl;
    // if (this->particles[0] != nullptr) {
    //     std::cout << this->particles[0] << std::endl;
    // } else {
    //     std::cout << "Null pointer" << std::endl;
    // }
    // define corner positions
    this->corners = {
        glm::vec2(center[0] - width / 2, center[1] - height / 2),  // lower left
        glm::vec2(center[0] + width / 2, center[1] - height / 2),  // lower right
        glm::vec2(center[0] + width / 2, center[1] + height / 2),  // upper right
        glm::vec2(center[0] - width / 2, center[1] + height / 2)   // upper left
    };

    this->char_len = std::max(width,height);

}

    // uninitialized submesh should be empty vector
    // this->sub_mesh = nullptr;


void NbodyMesh::update_mesh(){
    // update which mesh particles are in

    // current does nothing
    // std::cout << "modify mesh" << std::endl;
    
    // this -> sub_mesh.resize(4,nullptr);
    this -> modify_mesh();

    // update mesh and particle assignments
    // std::cout << "update particles" << std::endl;
    this -> update_particles();
    // std::cout << "update particles" << std::endl;
 
}

void NbodyMesh::update_particles(){
    // update which mesh particles are in and mesh states

    // reset tree states
    this->reset_mesh_tree_states();

    for (size_t i=0; i < this->particles.size();i++){
        // check if particle is active
        if (this->particles[i]!=nullptr){
            // set direction initially upward (using pos 1)
            // particle stores last lowest mesh hierarchy, use mesh_particle_updater to update particles and meshes
            if (this->particles[i]->mesh!=nullptr){
                this->particles[i]->mesh->mesh_particle_updater(particles[i], 1);  // Start by moving upward
            }
            else {
                // std::cout << "here" << std::endl;
                this->mesh_particle_updater(particles[i], 1);  // Start by moving upward
            }
            // check if particle was deactivated by mesh_particle updater
            if (this->particles[i]!=nullptr){
                if (this->particles[i]->mesh!=nullptr){
                    // update lowest mesh parameters
                    // all particles have shifted 

                    this->particles[i]->mesh->tot_mass += this->particles[i]->mass;
                    this->particles[i]->mesh->com += this->particles[i]->mass*this->particles[i]->pos;
                }
            }
        }
    }


    // propagate mass and com upward through tree
    this->propagate_mesh_tree_states();
}

void NbodyMesh::reset_mesh_tree_states(){
    // reset all particle defined states to 0
    this -> com = glm::vec2(0,0);
    this -> tot_mass = 0;
    this -> tot_energy = 0;
    this -> tot_ang_momentum = 0;
    this -> tot_momentum = 0;

    if (!this->sub_mesh.empty()){
        for (size_t i=0; i<this->sub_mesh.size(); i++){
            this->sub_mesh[i]->reset_mesh_tree_states();
        }
    }
}

void NbodyMesh::remove_particle(int ind){
    // removes particle from mesh starting at root
    this->particles[ind] = nullptr;
    if (this->parent_mesh != nullptr){
        this->parent_mesh->remove_particle(ind);
    }
}

void NbodyMesh::propagate_mesh_tree_states(){
    if (!this->sub_mesh.empty()){
        for (size_t i=0; i<this->sub_mesh.size(); i++){
            this->sub_mesh[i]->propagate_mesh_tree_states();
            if (this->sub_mesh[i]->tot_mass!=0){
                this->tot_mass += this->sub_mesh[i]->tot_mass;
                this->com += this->sub_mesh[i]->com*this->sub_mesh[i]->tot_mass;
            }
        }

        if (this->tot_mass != 0){
            this->com /= this->tot_mass;
        }
    }
    else {
        // lowest state must calculate com and total states
        // all total states should be summed when iterating through particles into lowest mesh
        // divide to get com
        if (this->tot_mass != 0){
            this->com /= this->tot_mass;
        }
        // std::cout << tot_mass << std::endl;
    }
}

void NbodyMesh::mesh_particle_updater(Particle* particle, int direction){
    // Updates particle and mesh assignments based off mesh tree
    // starts at last lowest mesh and search mesh tree upwards until within bounds
    // then search tree downward to find lowest mesh hierarchy and assigns pointer 
    // to particle
    // returns nullptr in particle if outside of all meshes

    // mesh is rectangular so only need to check against 2 points
    // Check if particle is outside the current mesh
    bool isOutside = particle->pos.x < this->corners[0].x ||
                     particle->pos.x > this->corners[2].x ||
                     particle->pos.y < this->corners[0].y ||
                     particle->pos.y > this->corners[2].y;

    if (isOutside)
        {
            // temporarily set particle mesh to null pointer
            particle->mesh = nullptr;
            // if updater moving upwards look into parent cases else break and return nothing
            if (direction == 1){
                // remove particle from current mesh if mesh has parent
                // highest mesh keeps track of all meshes to be able to reassign into a submesh

                // check if has parent
                if (this->parent_mesh != nullptr){

                    this->particles[particle->ind] = nullptr;

                    // update mesh particles
                    // update particle meshes
                    // move upwards in mesh tree using 1
                    this->parent_mesh->mesh_particle_updater(particle,1);
                }
                else{
                    // if particle not in any mesh, assigned null ptr 
                    particle->mesh = nullptr; 
                }
            }
        }
    else {
            // std::cout << "here1" << std::endl;
            // std::cout << "a" << std::endl;
            // std::cout << particle << " : " << particle->ind << std::endl;
            // std::cout << this->particles[0] << std::endl;

            // std::cout << this->particles[particle->ind] << std::endl;
            // std::cout << "a" << std::endl;

            // assign particle to current mesh if in mesh
            this->particles[particle->ind] = particle;
            // std::cout << "c" << std::endl;
            // std::cout << sub_mesh.size() << std::endl;
            // std::cout << "c" << std::endl;


            // check if lowest and assigns mesh to particle if lowest heirachy mesh
            if (this->sub_mesh.empty()){
                // std::cout << "here2" << std::endl;

                this->particles[particle->ind]->mesh=this;
            }
            else {
                // else iterate through lower meshes
                for (size_t i=0; i<this->sub_mesh.size(); i++){
                    // move downward through submeshes
                    this->sub_mesh[i]->mesh_particle_updater(particle,-1);
                    if (particle->mesh != nullptr){
                        // particle assigned to mesh
                        return;
                    }
                }
            }
        }
}