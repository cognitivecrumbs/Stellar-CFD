#include <particles.h>

Particle::Particle(int ind, glm::vec2 pos, glm::vec2 vel, glm::vec2 acc, float mass, float internal_energy, float ang_mom):
ind(ind), pos(pos), vel(vel), acc(acc),buffer_pos(pos), buffer_vel(vel), buffer_acc(acc), mass(mass), internal_energy(internal_energy), ang_mom(ang_mom)
{
    this->calcRenderParams();
}

void Particle::updateFromBuffer(){
    this->pos = buffer_pos;
    this->vel = buffer_vel;
    this->acc = buffer_acc;
}

void Particle::calcRenderParams(){
    // calculate render point size
    // this -> render_size = this->mass;
    // this -> render_size = this->mass/(cbrt(this->mass)*cbrt(this->mass));
    this -> render_size = std::log(this->mass);
    if (this -> render_size<3){
        this -> render_size = 3;
    }
    // this -> render_size = 10;

    // calculate color
    // if (this->mesh==nullptr){
    //     this -> render_color = glm::vec3(1.,1.,1.);
    // }
    // else {
    //     // Hash the address to get a "unique" value
    //     uintptr_t intAddress = reinterpret_cast<uintptr_t>(this->mesh);

    //     intAddress = intAddress ^ (intAddress >> 33);

    //     // Extract three components from the integer address
    //     float r = static_cast<float>((intAddress & 0xFF)) / 255.0f;           // Red from lowest byte
    //     float g = static_cast<float>((intAddress >> 8) & 0xFF) / 255.0f;     // Green from next byte
    //     float b = static_cast<float>((intAddress >> 16) & 0xFF) / 255.0f;    // Blue from next byte

    //     this -> render_color = glm::vec3(r, g, b);
    // }

    // Clamp t to [0, 1]

    glm::vec3 color;
    if (this->render_size<4){
        float lum = std::pow(this->internal_energy/this->mass,0.2);
        float t = glm::clamp(lum/30, 0.0f, 1.0f);
        // if (t <= 0.6f) {
        //     float t_prime = t * 2.0f; // Scale t to [0, 1]
        //     this -> render_color = glm::mix(glm::vec3(0.2, 0, 0), glm::vec3(1, 1, 1), t_prime);
        // } else {
        //     float t_prime = (t - 0.6f) * 2.0f; // Scale t to [0, 1]
        //     this -> render_color = glm::mix(glm::vec3(1, 1, 1), glm::vec3(0.8, 0.95, 1), t_prime);
        // }
        this -> render_color = glm::mix(glm::vec3(0.3, 0.3, 0.4), glm::vec3(1, 0.8, 0.2), t);

    }
    else {
        // float r = cbrt(this->mass);
        // float t = glm::clamp(r*100+this->internal_energy, 0.0f, 1.0f);
        float t = glm::clamp(this->render_size/10, 0.0f, 1.0f);

        if (t<=0.6f){
            float t_prime = t * 2.0f; // Scale t to [0, 1]
            this -> render_color = glm::mix(glm::vec3(1., 1., 0), glm::vec3(1, 1, 1), t_prime);
        } else {
            float t_prime = (t - 0.6f) * 2.0f; // Scale t to [0, 1]
            this -> render_color = glm::mix(glm::vec3(1, 1, 1), glm::vec3(0.65, 0.95, 1), t_prime);
        }
    }

    // glm::vec3 color;
    // float vel = sqrt(dot(this->vel,this->vel));
    // float t = glm::clamp(vel/3000, 0.0f, 1.0f);

    // this -> render_color = glm::mix(glm::vec3(0.6, 0.1, 0.6), glm::vec3(1, 1., 0.2), t);


    // this->render_color = glm::vec3(1,1,1);
}