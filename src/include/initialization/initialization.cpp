#include <initialization.h>

void Initialization::init(std::vector <Particle*> *particles){}

DiskInitialization::DiskInitialization(float base_mass, float radius_max, float radius_min, float disk_ang_mom, float base_energy, float base_ang_momentum):
base_mass(base_mass),radius_max(radius_max),radius_min(radius_min),disk_ang_mom(disk_ang_mom),base_energy(base_energy),base_ang_momentum(base_ang_momentum)
{
}

void DiskInitialization::init(std::vector <Particle*> *particles){
    for (size_t i = 0; i < particles->size(); i++) {
        // Generate random radius and angle within the disk
        // float radius = radius_max * std::pow((static_cast<float>(rand()) / RAND_MAX),1.)+radius_min; // Uniform distribution in area
        float radius = radius_max * std::pow((static_cast<float>(rand()) / RAND_MAX),0.5)+radius_min; // Uniform distribution in area
        float angle = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;       // Random angle in [0, 2π]

        // Position in Cartesian coordinates
        glm::vec2 pos = radius * glm::vec2(cos(angle), sin(angle));

        // Random tangential velocity magnitude
        float random_factor = 0.5f + static_cast<float>(rand()) / RAND_MAX; // Random scaling factor in [0.5, 1.5]
        // float velocity_magnitude = radius * disk_ang_mom * random_factor;  // Tangential velocity scaled by radius and randomness
        float velocity_magnitude = 0.8*disk_ang_mom * random_factor/radius;  // Tangential velocity scaled by radius and randomness

        // Tangential velocity (perpendicular to radial vector)
        glm::vec2 vel = velocity_magnitude * glm::vec2(-sin(angle), cos(angle));

        // Zero initial acceleration
        glm::vec2 acc(0.0f, 0.0f);

        float mass_random_factor = 0.5f + static_cast<float>(rand()) / RAND_MAX; // Random mass scaling factor in [0.5, 1.5]
        float mass = base_mass * std::pow(4,(mass_random_factor-0.5)*4+0.5);

        // Initialize the particle
        (*particles)[i] = new Particle(i, pos, vel, acc, mass, base_energy, base_ang_momentum);
        (*particles)[i]->calcRenderParams();
    }
}

void DiskInitialization::reInitParticle(Particle* particle){
    // Generate random radius and angle within the disk
    float radius = 2.* std::pow((static_cast<float>(rand()) / RAND_MAX),1.)+radius_max; // Uniform distribution in area
    // float radius = this->radius_max*1.3;
    float angle = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;       // Random angle in [0, 2π]

    // Position in Cartesian coordinates
    glm::vec2 pos = radius * glm::vec2(cos(angle), sin(angle));

    // Random tangential velocity magnitude
    float random_factor = 0.5f + static_cast<float>(rand()) / RAND_MAX; // Random scaling factor in [0.5, 1.5]
    // float velocity_magnitude = radius * disk_ang_mom * random_factor;  // Tangential velocity scaled by radius and randomness
    float velocity_magnitude = 0.8*this->disk_ang_mom * random_factor/radius;  // Tangential velocity scaled by radius and randomness

    // Tangential velocity (perpendicular to radial vector)
    glm::vec2 vel = velocity_magnitude * glm::vec2(-sin(angle), cos(angle));

    // // Zero initial acceleration
    // glm::vec2 acc(0.0f, 0.0f);

    float mass_random_factor = 0.5f + static_cast<float>(rand()) / RAND_MAX; // Random mass scaling factor in [0.5, 1.5]
    float mass = this->base_mass * std::pow(2,mass_random_factor);

    // reinitialize
    particle -> pos = pos;
    particle -> vel = vel;
    particle -> acc = glm::vec2(0,0);    
    particle -> buffer_pos = pos;
    particle -> buffer_vel = vel;
    particle -> buffer_acc = glm::vec2(0,0);
    particle -> mass = mass;
    particle -> internal_energy = this->base_energy;
    particle -> ang_mom = base_ang_momentum;
}