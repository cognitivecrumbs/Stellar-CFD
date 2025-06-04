#include <mesh.h>
#include <initialization.h>

// ---------------------- EIGEN USE -----------------------------------------

// #include <mesh.h>
// #include <initialization.h>

CFDMeshWrapper::CFDMeshWrapper(int nx, int ny, int ghost_cell, float dx) {
    // float dx = l / n;
    vector <float> U(4,0);

    // Resize mesh and initialize cells
    this->mesh.resize((nx+ghost_cell*2) * (ny+ghost_cell*2), nullptr);
    for (int i = 0; i < ny + ghost_cell*2; i++) {
        for (int j = 0; j < nx + ghost_cell*2; j++) {
            glm::vec2 center(-ghost_cell*dx + dx / 2 + dx * j, -ghost_cell*dx + dx / 2 + dx * i);
            this->mesh[i * (nx+2*ghost_cell) + j] = new CFDMesh(i * (nx+2*ghost_cell) + j, center, dx, U);
        }
    }

    // Create connections
    for (int i = 0; i < ny + ghost_cell*2; i++) {
        for (int j = 0; j < nx + ghost_cell*2; j++) {
            for (int k = 0; k < 4; k++) {
                // Setup boundary conditions (adiabatic walls)
                if ((i == 0 && k == 3) || (i == ny + 2*ghost_cell - 1 && k == 1) ||
                    (j == 0 && k == 2) || (j == nx + 2*ghost_cell - 1 && k == 0)) {
                    this->mesh[i * (nx+2*ghost_cell) + j]->neighbor[k] = nullptr;
                    continue;
                }

                // Internal connectivity
                if (k == 0) { // Right neighbor
                    this->mesh[i * (nx+2*ghost_cell) + j]->neighbor[k] = this->mesh[i * (nx+2*ghost_cell) + j + 1];
                } else if (k == 1) { // Bottom neighbor
                    this->mesh[i * (nx+2*ghost_cell) + j]->neighbor[k] = this->mesh[(i + 1) * (nx+2*ghost_cell) + j];
                } else if (k == 2) { // Left neighbor
                    this->mesh[i * (nx+2*ghost_cell) + j]->neighbor[k] = this->mesh[i * (nx+2*ghost_cell) + j - 1];
                } else if (k == 3) { // Top neighbor
                    this->mesh[i * (nx+2*ghost_cell) + j]->neighbor[k] = this->mesh[(i - 1) * (nx+2*ghost_cell) + j];
                }
            }
        }
    }

    // Initialize the mesh using a specific initializer
    this->initializer = new ObliqueShock();
    this->initializer = new KelvinHelholtz();
    this->initializer = new GravPressure();

    this->initializer->init(this->mesh, nx, ny, ghost_cell);
    this->ghost_cell = ghost_cell;
    this->nx = nx; this->ny = ny;

    // delete initializer; // Cleanup to prevent memory leaks
    std::vector<int> active_ind;
    active_ind.resize(this->mesh.size(),0);
    for (int i=0;i<int(this->mesh.size());i++){
        active_ind[i] = i;
    }
    this->grav_mesh = new BarnesHutCFDMesh(this->mesh,nx*dx, ny*dx, glm::vec2(float(nx)/2*dx,float(ny)/2*dx), active_ind, (nx+ghost_cell*2)*(nx+ghost_cell*2), nullptr);
    this->grav_mesh->initializeMesh();
}

void CFDMeshWrapper::update_mesh() {
    float dt = 0.0011; // Use a fixed timestep
    
    const int nx_total = this->nx + 2 * this->ghost_cell;
    const int x_boundary_min = this->ghost_cell;
    const int x_boundary_max = this->nx + this->ghost_cell;
    const int y_boundary_min = this->ghost_cell;
    const int y_boundary_max = this->ny + this->ghost_cell;

    float p,a,u,v,V; 

    ThreadPool pool(8);

    // dt = 1e8;
    float indiv_dt;
    // Phase 1: Minmod Slope and initialize dUdt
    // this->grav_mesh->reset_mesh_tree_states();
    this->grav_mesh->update_mesh_mass();
    for (CFDMesh* mesh_ele : this->mesh) {
        mesh_ele->fx = 0;
        mesh_ele->fy = 0;
        // this->grav_mesh->applyGrav(mesh_ele);
        pool.enqueue([this,mesh_ele](){
            this->grav_mesh->applyGrav(mesh_ele);
        });
    }
    pool.stopPool();

    pool.restart(8);
    for (CFDMesh* mesh_ele : this->mesh) {
        const int ind = mesh_ele->ind;
        const int x_index = ind % nx_total;
        const int y_index = ind / nx_total;


        if (x_index < x_boundary_min || x_index >= x_boundary_max ||
            y_index < y_boundary_min || y_index >= y_boundary_max) {
            continue;
            }
        u = mesh_ele->U[1]/mesh_ele->U[0]; v = mesh_ele->U[2]/mesh_ele->U[0];
        V = (u * u + v * v);
        p = (mesh_ele->gamma - 1) * (mesh_ele->U[3] - 0.5f * mesh_ele->U[0] * V);
        indiv_dt = 0.8*mesh_ele->dx/(V+sqrt(mesh_ele->gamma*p/mesh_ele->U[0]));
        dt = min(dt,indiv_dt);

        // mesh_ele->minmodSlope();
        // mesh_ele->calcSlope();
        pool.enqueue([mesh_ele](){
            mesh_ele->minmodSlope();
        });
        for (int ind=0;ind<4;ind++){
            mesh_ele->dUdt[ind] = 0;
        }

    }
    pool.stopPool();

    // Phase 2: Compute fluxes
    pool.restart(8);
    for (CFDMesh* mesh_ele : this->mesh) {
        const int ind = mesh_ele->ind;
        const int x_index = ind % nx_total;
        const int y_index = ind / nx_total;

        if (x_index < x_boundary_min - 1 || x_index >= x_boundary_max ||
            y_index < y_boundary_min - 1 || y_index >= y_boundary_max) {
            continue;
        }

        // if (mesh_ele->ind == (25+2)*(104)+3){
        //     cout << "tmp" << endl;
        // }
        // mesh_ele->computeHLLFlux();
        pool.enqueue([mesh_ele](){
            mesh_ele->computeHLLFlux();
        });

        // mesh_ele->computeFluxCentralDifference();

    }
    pool.stopPool();

    // Phase 3: Update solution
    for (CFDMesh* mesh_ele : this->mesh) {
        const int ind = mesh_ele->ind;
        const int x_index = ind % nx_total;
        const int y_index = ind / nx_total;

        if (x_index < x_boundary_min || x_index >= x_boundary_max ||
            y_index < y_boundary_min || y_index >= y_boundary_max) {
            continue;
        }

        for (int ind=0;ind<4;ind++){
            mesh_ele->U[ind] += mesh_ele->dUdt[ind] * dt;
        }
    }

    this->initializer->updateBC(this->mesh,this->nx,this->ny,this->ghost_cell);
}

CFDMesh::CFDMesh(int ind, glm::vec2 center, float dx, vector <float> U):
center(center), dx(dx), U(U), ind(ind) {
    this->vol = dx * dx;
    this->prim = Eigen::Matrix<float, 5, 1>::Zero();

    this->neighbor.resize(4, nullptr);
}

void CFDMesh::cons2prim() {
    this->prim[0] = this->U[0];
    this->prim[1] = this->U[1]/this->U[0];
    this->prim[2] = this->U[2]/this->U[0];

    this->vel_squared = this->prim(seq(1, 2)).squaredNorm();
    this->prim(3) = this->prim(3) / this->prim(0) - this->vel_squared / 2;

    float tmp_p = (gamma - 1) * (this->U[3] - 0.5 * this->U[0] * this->prim(seq(1, 2)).squaredNorm());
    this->prim(4) = tmp_p;
}

void CFDMesh::prim2cons() {
    this->U[0] = this->prim[0];
    this->U[1] = this->prim[0]*this->prim[1];
    this->U[2] = this->prim[0]*this->prim[2];
    this->U[3] = this->prim[4] / (this->gamma - 1) + (this->prim(seq(1, 2)).squaredNorm() / 2) * this->U[0];
}


// Minmod Limiter
float CFDMesh::minmod(float a, float b, float c) {
    float v = min(min(abs(a),abs(b)),abs(c));

    if ( (a<0&&b<0&&c<0) || (a>0&&b>0&&c>0)){
        if (a>0&&b>0&&c>0){
            return v;
        }
        else {
            return -v;
        }
    }
    else {
        return 0;
    }

}

void CFDMesh::minmodSlope() {
    const float dx_inv = 1.0f / this->dx; // Precompute reciprocal of dx for reuse.

    // // Compute slopes for the x-direction
    // const Vector4f slope_R_x = (this->neighbor[0]->U - this->U) * dx_inv;
    // const Vector4f slope_L_x = (this->U - this->neighbor[2]->U) * dx_inv;
    // const Vector4f slope_C_x = (this->neighbor[0]->U - this->neighbor[2]->U) * dx_inv / 2.0f;

    // for (int i = 0; i < 4; ++i) {
    //     this->slope_x[i] = minmod(slope_L_x[i], slope_R_x[i], slope_C_x[i]);
    // }

    // // Compute slopes for the y-direction
    // const Vector4f slope_R_y = (this->neighbor[1]->U - this->U) * dx_inv;
    // const Vector4f slope_L_y = (this->U - this->neighbor[3]->U) * dx_inv;
    // const Vector4f slope_C_y = (this->neighbor[1]->U - this->neighbor[3]->U) * dx_inv / 2.0f;

    // for (int i = 0; i < 4; ++i) {
    //     this->slope_y[i] = minmod(slope_L_y[i], slope_R_y[i], slope_C_y[i]);
    // }


    float slope_R,slope_L,slope_C;

    for (int ind=0;ind<4;ind++){

        // Compute slopes for the x-direction
        slope_R = (this->neighbor[0]->U[ind] - this->U[ind]) * dx_inv;
        slope_L = (this->U[ind] - this->neighbor[2]->U[ind]) * dx_inv;
        slope_C = (this->neighbor[0]->U[ind] - this->neighbor[2]->U[ind]) * dx_inv / 2.0f;

        this->slope_x[ind] = minmod(slope_L, slope_R, slope_C);

        // Compute slopes for the y-direction
        slope_R = (this->neighbor[1]->U[ind] - this->U[ind]) * dx_inv;
        slope_L = (this->U[ind] - this->neighbor[3]->U[ind]) * dx_inv;
        slope_C = (this->neighbor[1]->U[ind] - this->neighbor[3]->U[ind]) * dx_inv / 2.0f;

        this->slope_y[ind] = minmod(slope_L, slope_R, slope_C);
    }
}

void CFDMesh::calcSlope() {
    // Compute slopes for the x-direction
    for (int ind=0;ind<4;ind++){
        this->slope_x[ind] = (this->neighbor[2]->U[ind] - this->neighbor[0]->U[ind]) / (2*this->dx);

        // Compute slopes for the y-direction
        this->slope_y[ind] = (this->neighbor[3]->U[ind] - this->neighbor[1]->U[ind]) / (2*this->dx);
    }
}


void CFDMesh::computeHLLFlux() {
    const float dx_half = this->dx * 0.5f; // Precompute for efficiency.

    for (size_t i = 0; i < this->neighbor.size(); ++i) {
        CFDMesh* neighbor = this->neighbor[i];
        // if (!neighbor || neighbor->ind <= this->ind) continue;

        vector <float> U_L(4), U_R(4);
        int normal_x = 0, normal_y = 0;

        for (int ind=0;ind<4;ind++){
            // Select face and compute U_L, U_R based on the neighbor index
            switch (i) {
            case 0: // Left face
                U_L[ind] = this->U[ind] + this->slope_x[ind] * dx_half;
                U_R[ind] = neighbor->U[ind] - neighbor->slope_x[ind] * dx_half;
                normal_x = 1;
                break;
            case 1: // Bottom face
                U_L[ind] = this->U[ind] + this->slope_y[ind] * dx_half;
                U_R[ind] = neighbor->U[ind] - neighbor->slope_y[ind] * dx_half;
                normal_y = 1;
                break;
            case 2: // Right face
                U_L[ind] = neighbor->U[ind] + neighbor->slope_x[ind] * dx_half;
                U_R[ind] = this->U[ind] - this->slope_x[ind] * dx_half;
                normal_x = -1;
                break;
            case 3: // Top face
                U_L[ind] = neighbor->U[ind] + neighbor->slope_y[ind] * dx_half;
                U_R[ind] = this->U[ind] - this->slope_y[ind] * dx_half;
                normal_y = -1;
                break;
            }
        }

        // Compute HLLE flux
        vector <float> flux(4,0);
        this->HLLE(flux, U_L, U_R, normal_x, normal_y);
        // this->Rusanov(flux, U_L, U_R, normal_x, normal_y);

        // Update dUdt for this cell and the neighbor
        for (int ind=0;ind<4;ind++){
            this->dUdt[ind] -= flux[ind] * this->dx;
        }
        // neighbor->dUdt += flux * this->dx;
    }

    for (int ind=0;ind<4;ind++){
        this->dUdt[ind] /= this->vol; // volume normalize
    }
    this->bodyForce();
}

void CFDMesh::HLLE(vector <float>& flux, const vector <float>& qL, const vector <float>& qR, int nx, int ny) {
    // Extract state variables for left and right
    const float rhoL = qL[0], rhoR = qR[0];
    const float uL = qL[1] / rhoL, uR = qR[1] / rhoR;
    const float vL = qL[2] / rhoL, vR = qR[2] / rhoR;
    const float pL = (gamma - 1) * (qL[3] - 0.5f * rhoL * (uL * uL + vL * vL));
    const float pR = (gamma - 1) * (qR[3] - 0.5f * rhoR * (uR * uR + vR * vR));
    const float vnL = uL * nx + vL * ny;
    const float vnR = uR * nx + vR * ny;

    // Compute sound speed and enthalpy
    const float aL = std::sqrt(gamma * pL / rhoL);
    const float aR = std::sqrt(gamma * pR / rhoR);
    const float HL = (qL[3] + pL) / rhoL;
    const float HR = (qR[3] + pR) / rhoR;

    // Roe averages
    const float RT = std::sqrt(rhoR / rhoL);
    const float u = (uL + RT * uR) / (1.0f + RT);
    const float v = (vL + RT * vR) / (1.0f + RT);
    const float H = (HL + RT * HR) / (1.0f + RT);
    const float a = std::sqrt((gamma - 1) * (H - 0.5f * (u * u + v * v)));
    const float vn = u * nx + v * ny;

    // Wave speed estimates
    const float SLm = std::min(vnL - aL, vn - a);
    const float SRp = std::max(vnR + aR, vn + a);

    if (SLm >= SRp) {
        // flux.setZero();
        for (int ind=0;ind<4;ind++){
            flux[ind] = 0;
        }
        return;
    }

    // // Compute fluxes for left and right states
    // const Vector4f FL = {rhoL * vnL,
    //                      rhoL * vnL * uL + pL * nx,
    //                      rhoL * vnL * vL + pL * ny,
    //                      rhoL * vnL * HL};
    // const Vector4f FR = {rhoR * vnR,
    //                      rhoR * vnR * uR + pR * nx,
    //                      rhoR * vnR * vR + pR * ny,
    //                      rhoR * vnR * HR};

    // // HLLE flux
    // int sign_change = 1;
    // if (ny<0 || nx<0){
    //     sign_change = -1;
    // }
    // flux = (SRp * FL - SLm * FR + sign_change*SLm * SRp * (qR - qL)) / (SRp - SLm);

    int sign_change = 1;
    if (ny<0 || nx<0){
        sign_change = -1;
    }
    flux[0] = SRp*(rhoL*vnL) - SLm*(rhoR*vnR) + sign_change*SLm*SRp*(qR[0]-qL[0])/(SRp-SLm);
    flux[1] = SRp*(rhoL * vnL * uL + pL * nx) - SLm*(rhoR * vnR * uR + pR * nx) + sign_change*SLm*SRp*(qR[1]-qL[1])/(SRp-SLm);
    flux[2] = SRp*(rhoL * vnL * vL + pL * ny) - SLm*(rhoR * vnR * vR + pR * ny) + sign_change*SLm*SRp*(qR[2]-qL[2])/(SRp-SLm);
    flux[3] = SRp*(rhoL * vnL * HL) - SLm*(rhoR * vnR * HR) + sign_change*SLm*SRp*(qR[3]-qL[3])/(SRp-SLm);
}


void CFDMesh::bodyForce() {
    // float fx = 0;
    // float fy = -1;
    this->dUdt[1] += this->U[0]*fx;
    this->dUdt[2] += this->U[0]*fy;
    this->dUdt[3] += this->U[0]*(this->U[1]/this->U[0]*fx+this->U[2]/this->U[0]*fy);

    this->dUdt[3] += 0.1*(this->U[0]-2);
}


// void CFDMesh::Rusanov(Vector4f& flux, const Vector4f& qL, const Vector4f& qR, int nx, int ny) {
//     // Extract state variables for left and right
//     const float rhoL = qL[0], rhoR = qR[0];
//     const float uL = qL[1] / rhoL, uR = qR[1] / rhoR;
//     const float vL = qL[2] / rhoL, vR = qR[2] / rhoR;
//     const float pL = (gamma - 1) * (qL[3] - 0.5f * rhoL * (uL * uL + vL * vL));
//     const float pR = (gamma - 1) * (qR[3] - 0.5f * rhoR * (uR * uR + vR * vR));
//     const float vnL = uL * nx + vL * ny;
//     const float vnR = uR * nx + vR * ny;

//     // Compute sound speed and enthalpy
//     const float aL = std::sqrt(gamma * pL / rhoL);
//     const float aR = std::sqrt(gamma * pR / rhoR);
//     const float HL = (qL[3] + pL) / rhoL;
//     const float HR = (qR[3] + pR) / rhoR;

//     // Roe averages
//     const float RT = std::sqrt(rhoR / rhoL);
//     const float u = (uL + RT * uR) / (1.0f + RT);
//     const float v = (vL + RT * vR) / (1.0f + RT);
//     const float H = (HL + RT * HR) / (1.0f + RT);
//     const float a = std::sqrt((gamma - 1) * (H - 0.5f * (u * u + v * v)));
//     const float vn = u * nx + v * ny;

//     // Wave speed estimates
//     const float SLm = std::min(vnL - aL, vn - a);
//     const float SRp = std::max(vnR + aR, vn + a);

//     if (SLm >= SRp) {
//         flux.setZero();
//         return;
//     }

//     // Compute fluxes for left and right states
//     const Vector4f FL = {rhoL * vnL,
//                          rhoL * vnL * uL + pL * nx,
//                          rhoL * vnL * vL + pL * ny,
//                          rhoL * vnL * HL};
//     const Vector4f FR = {rhoR * vnR,
//                          rhoR * vnR * uR + pR * nx,
//                          rhoR * vnR * vR + pR * ny,
//                          rhoR * vnR * HR};

//     // HLLE flux
//     // flux = (SRp * FL - SLm * FR + SLm * SRp * (qR - qL)) / (SRp - SLm);

//     flux = 0.5*(FL+FR) - 0.5*max(abs(SLm),abs(SRp))*(qR-qL);
// }

// void CFDMesh::computeFlux() {
//     // Initialize flux computation
//     this->dUdt = Vector4f::Zero();
//     bool this_upwind;


//     for (size_t i = 0; i < this->neighbor.size(); i++) {
//         int dir = (i == 0 || i == 1) ? 1 : -1;
//         float vel_avg;
//         int face_orientation = i % 2; // 0: horizontal, 1: vertical

//         if (this->neighbor[i] == nullptr) {
//             // Handle boundary conditions (e.g., adiabatic wall)
//             float p = (gamma - 1) * (this->U[3] - 0.5 * (this->U[1] * this->U[1] + this->U[2] * this->U[2]) / this->U[0]);
//             this->dUdt[1 + face_orientation] -= dir * (this->U[face_orientation+1]*this->U[face_orientation+1] /this->U[0] + p) * dx;
//             continue;
//         }


//         // Compute average velocity and pressure
//         vel_avg = 0.5 * (this->U[1 + face_orientation] / this->U[0] +
//                          this->neighbor[i]->U[1 + face_orientation] / this->neighbor[i]->U[0]);
//         // vel_avg = 0.5 * (this->U[1 + face_orientation]  +
//         //                  this->neighbor[i]->U[1 + face_orientation]);
//         if (vel_avg == 0){
//             float p1 = (gamma - 1) * (this->U[3] - 0.5 * (this->U[1] * this->U[1] + this->U[2] * this->U[2]) / this->U[0]);
//             float p2 = (gamma - 1) * (this->neighbor[i]->U[3] - 0.5 * (this->neighbor[i]->U[1] * this->neighbor[i]->U[1] + this->neighbor[i]->U[2] * this->neighbor[i]->U[2]) / this->neighbor[i]->U[0]);
//             this->dUdt -= dir * computeFluxForFace((this->U+this->neighbor[i]->U)/2, (p1+p2)/2, vel_avg, dx, face_orientation);
//             continue;
//         }
//         this_upwind = (vel_avg > 0 && i<=1) || (vel_avg<0 && i>=2);


//         // this_upwind = (vel_avg < 0 && i<=1) || (vel_avg>0 && i>=2);

//         Vector4f flux;
//         if (this_upwind) {
//             float p = (gamma - 1) * (this->U[3] - 0.5 * (this->U[1] * this->U[1] + this->U[2] * this->U[2]) / this->U[0]);
//             flux = computeFluxForFace(this->U, p, vel_avg, dx, face_orientation);
//             // flux = computeFluxForFace(this->U, p, this->neighbor[i]->U[1+face_orientation]/this->neighbor[i]->U[0], dx, face_orientation);
//             // flux = computeFluxForFace(this->U, p, this->U[1+face_orientation]/this->U[0], dx, face_orientation);
//         } else {
//             float p = (gamma - 1) * (this->neighbor[i]->U[3] - 0.5 * (this->neighbor[i]->U[1] * this->neighbor[i]->U[1] + this->neighbor[i]->U[2] * this->neighbor[i]->U[2]) / this->neighbor[i]->U[0]);
//             flux = computeFluxForFace(this->neighbor[i]->U, p, vel_avg, dx, face_orientation);
//             // flux = computeFluxForFace(this->neighbor[i]->U, p, this->U[1+face_orientation]/this->U[0], dx, face_orientation);
//             // flux = computeFluxForFace(this->neighbor[i]->U, p, this->neighbor[i]->U[1+face_orientation]/this->neighbor[i]->U[0], dx, face_orientation);
//         }

//         this->dUdt -= dir * flux;
//     }

//     this->dUdt /= this->vol; // Normalize by cell volume
// }


// void CFDMesh::computeFluxCentralDifference() {
//     // Initialize the flux derivative term
//     this->dUdt = Vector4f::Zero();

//     for (size_t i = 0; i < this->neighbor.size(); i++) {
//         if (this->neighbor[i] == nullptr) {
//             // Apply boundary conditions for adiabatic no-slip walls
//             float p = (gamma - 1) * (this->U[3] - 0.5 * (this->U[1] * this->U[1] + this->U[2] * this->U[2]) / this->U[0]);
//             int dir = (i >= 2) ? 1 : -1;
//             if (i == 0 || i == 2) {
//                 this->dUdt[1] += dir * p * dx;
//             } else {
//                 this->dUdt[2] += dir * p * dx;
//             }
//             continue;
//         }

//         // Central difference flux computation
//         Vector4f U_left = this->U;                    // Current cell's state
//         Vector4f U_right = this->neighbor[i]->U;      // Neighboring cell's state
//         // Vector4f U_face = 0.5 * (U_left + U_right);   // State at the face (central average)
//         Vector4f U_face = 0.5 * U_left + 0.5 * U_right;   // State at the face (central average)

//         // Compute flux components based on Euler equations
//         float rho = U_face[0];
//         float u = U_face[1] / rho;
//         float v = U_face[2] / rho;
//         float energy = U_face[3];
//         float pressure = (gamma - 1) * (energy - 0.5 * rho * (u * u + v * v));

//         Vector4f flux_face;
//         if (i == 0 || i == 2) { // Horizontal flux
//             flux_face[0] = rho * u;                    // Mass flux
//             flux_face[1] = rho * u * u + pressure;     // Momentum in x
//             flux_face[2] = rho * u * v;                // Momentum in y
//             flux_face[3] = (energy + pressure) * u;    // Energy flux
//         } else { // Vertical flux
//             flux_face[0] = rho * v;                    // Mass flux
//             flux_face[1] = rho * u * v;                // Momentum in x
//             flux_face[2] = rho * v * v + pressure;     // Momentum in y
//             flux_face[3] = (energy + pressure) * v;    // Energy flux
//         }

//         // Adjust flux derivative for this cell
//         int dir = (i == 0 || i == 1) ? 1 : -1; // Direction: positive for right/bottom, negative for left/top
//         this->dUdt -= dir * flux_face * dx;
//     }

//     // Divide by the volume to get the rate of change of conserved quantities
//     this->dUdt /= this->vol;
// }

// Vector4f CFDMesh::computeFluxForFace(const Vector4f& U, float p, float vel, float dx, int face_orientation) {
//     Vector4f flux;
//     if (face_orientation == 0) { // Horizontal flux
//         flux[0] = U[1];
//         flux[1] = U[1] * vel + p;
//         flux[2] = U[2] * vel;
//         flux[3] = U[3] * vel + p * vel;
//     } else { // Vertical flux
//         flux[0] = U[2];
//         flux[1] = U[1] * vel;
//         flux[2] = U[2] * vel + p;
//         flux[3] = U[3] * vel + p * vel;
//     }
//     return flux * dx;
// }



BarnesHutCFDMesh::BarnesHutCFDMesh(std::vector<CFDMesh*> cfd_meshes,float width, float height, glm::vec2 center, std::vector<int> active_ind, int active, BarnesHutCFDMesh* parent_mesh):
cfd_meshes(cfd_meshes),center(center),com(glm::vec2(0.,0.)),tot_mass(0.),parent_mesh(parent_mesh),width(width),height(height),active_ind(active_ind),npart_active(active)
{
    // define corner positions
    this->corners = {
        glm::vec2(center[0] - width / 2, center[1] - height / 2),  // lower left
        glm::vec2(center[0] + width / 2, center[1] - height / 2),  // lower right
        glm::vec2(center[0] + width / 2, center[1] + height / 2),  // upper right
        glm::vec2(center[0] - width / 2, center[1] + height / 2)   // upper left
    };

    this->char_len = std::max(width,height);

}

void BarnesHutCFDMesh::initializeMesh(){
    // this->deleteTree();
    // remove references in cfd_meshes to old mesh

    // create new mesh
    // ThreadPool pool(8);
    this->checkActive();
    this->splitMesh();
}

void BarnesHutCFDMesh::update_mesh_mass(){
    // update which mesh cfd_meshes are in
    // generates or modifies mesh
    // this -> modify_mesh();

    // update mass and com
    this->tot_mass = 0;
    this->com = glm::vec2(0,0);
    if (!this->sub_mesh.empty()){
        for (int i = 0; i < 4; i++) {
            // Update mass and center of mass
            this->sub_mesh[i]->update_mesh_mass();
            this->tot_mass += this->sub_mesh[i]->tot_mass;
            this->com += this->sub_mesh[i]->tot_mass * this->sub_mesh[i]->com;
        }
        this -> com /= tot_mass;
    }
    else {
        // this -> com = glm::vec2(0,0);
        // this -> tot_mass = 0;
        for (int ind: this->active_ind){
            float mass = this->cfd_meshes[ind]->U[0]*this->cfd_meshes[ind]->vol;
            this -> com += this->cfd_meshes[ind]->center * mass;
            this -> tot_mass += mass;
        }
        this -> com /= tot_mass;
    }
}


void BarnesHutCFDMesh::deleteTree(){
    if (!this->sub_mesh.empty()){
        for (BarnesHutCFDMesh* mesh: this->sub_mesh){
            mesh->deleteTree();
            delete mesh;
        }
    }
    this->sub_mesh.clear();  // Clear the vector after deletion
}


void BarnesHutCFDMesh::splitMesh(){
    this -> tot_mass = 0;
    this -> com = glm::vec2(0,0);
    // if (this->active>100){
    if (this->npart_active>16){
        if (this->sub_mesh.empty()) {
            this->sub_mesh.resize(4, nullptr);  // Ensure submesh is resized only once
        }        
        
        // split mesh
        float submesh_width = this->width/2;
        float submesh_height = this->height/2;    
        

        // for (int i=0; i<4; i++){
        //     glm::vec2 submesh_center = 0.5f*(this->center+this->corners[i]);

        //     // this->sub_mesh[i] = new BarnesHutMesh(this->cfd_meshes,submesh_width,submesh_height,submesh_center,this->active_ind);
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

            threads.emplace_back([this, i, submesh_width, submesh_height, submesh_center]() {
                this->sub_mesh[i] = new BarnesHutCFDMesh(this->cfd_meshes,submesh_width,submesh_height,submesh_center,this->active_ind,this->npart_active,this);
    
                this->sub_mesh[i]->checkActive();

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

        if (this->tot_mass!=0){
            this->com /= tot_mass;
        }

    }
    else {
        // assign values

        if (npart_active!=0){

            // if (!this->sub_mesh.empty()){
            //     for (BarnesHutMesh* mesh: this->sub_mesh){
            //         mesh->deleteTree();
            //         delete mesh;
            //     }
            //     this->sub_mesh.clear();
            // }
            
            this -> com = glm::vec2(0,0);
            this -> tot_mass = 0;
            for (int ind: this->active_ind){
                float mass = this->cfd_meshes[ind]->U[0]*this->cfd_meshes[ind]->vol;
                this -> com += this->cfd_meshes[ind]->center * mass;
                this -> tot_mass += mass;
            }
            this -> com /= tot_mass;
        }
    }
}


void BarnesHutCFDMesh::checkActive(){

    this->npart_active = 0;

    auto isInsideBoundary = [this](int ind) {
        return this->cfd_meshes[ind] != nullptr &&
               this->cfd_meshes[ind]->center.x > this->corners[0].x &&
               this->cfd_meshes[ind]->center.x <= this->corners[2].x &&
               this->cfd_meshes[ind]->center.y > this->corners[0].y &&
               this->cfd_meshes[ind]->center.y <= this->corners[2].y;
    };

    // Remove indices for cfd_meshes that are outside the boundary
    std::vector<int>::iterator it = std::remove_if(
        this->active_ind.begin(), this->active_ind.end(),
        [this, &isInsideBoundary](int ind) -> bool {
            if (!isInsideBoundary(ind)) {
                // Set cfd_meshes outside the boundary to nullptr
                this->cfd_meshes[ind] = nullptr;
                return true; // Mark for removal
            }
            return false; // Keep if inside the boundary
        });

    // Erase the removed elements
    this->active_ind.erase(it, this->active_ind.end());

    // Count remaining active cfd_meshes
    this->npart_active = int(this->active_ind.size());
}

void BarnesHutCFDMesh::applyGrav(CFDMesh* cfd_mesh){
    // Check if particle is outside the current mesh
    float dist;
    float G = 0.4;
    bool isOutside = cfd_mesh->center.x < this->corners[0].x ||
                     cfd_mesh->center.x > this->corners[2].x ||
                     cfd_mesh->center.y < this->corners[0].y ||
                     cfd_mesh->center.y > this->corners[2].y;
    if (isOutside){
        glm::vec2 rel_pos =  this->center - cfd_mesh->center;
        bool use_grav_mesh = false;

        // if (this->char_len/(abs(rel_pos.x)+abs(rel_pos.y))<1.0){
        //     dist = sqrt(dot(rel_pos,rel_pos));
        //     if (this->char_len/dist<1.0){
        //         use_grav_mesh = true;
        //     }
        // }

        // dist = sqrt(dot(rel_pos,rel_pos));
        // if (this->char_len/dist<0.8){
        //     use_grav_mesh = true;
        // }

        // if (use_grav_mesh){
        // rel_pos = this->com - cfd_mesh->center;
        dist = sqrt(dot(rel_pos,rel_pos));
        if (this->char_len/dist<1.0){
            rel_pos = this->com - cfd_mesh->center;
            dist = sqrt(dot(rel_pos,rel_pos));
            cfd_mesh->fx += G*this->tot_mass/(dist*dist*dist)*rel_pos.x;
            cfd_mesh->fy += G*this->tot_mass/(dist*dist*dist)*rel_pos.y;
        }
        else {
            // check submeshes
            if (!this->sub_mesh.empty()){
                for (BarnesHutCFDMesh* sub_mesh_curr: this->sub_mesh){
                    sub_mesh_curr->applyGrav(cfd_mesh);
                }
            }
            else {
                // loop apply force if not within tolerance
                float mass_j;
                for (CFDMesh* mesh_j:this->cfd_meshes){
                    if (mesh_j!= nullptr){
                        glm::vec2 rel_pos = mesh_j->center - cfd_mesh->center;
                        dist = sqrt(dot(rel_pos,rel_pos));
                        mass_j = mesh_j->U[0]*mesh_j->vol;
                        cfd_mesh->fx += G*mass_j/(dist*dist*dist)*rel_pos.x;
                        cfd_mesh->fy += G*mass_j/(dist*dist*dist)*rel_pos.y;
                    }
                }
            }
        }
    }
    else {
        // check submeshes
        if (!this->sub_mesh.empty()){
            for (BarnesHutCFDMesh* sub_mesh_curr: this->sub_mesh){
                sub_mesh_curr->applyGrav(cfd_mesh);
            }
        }
        else {
            // loop apply force if not within tolerance
            float mass_j;
            for (CFDMesh* mesh_j:this->cfd_meshes){
                if (mesh_j!= nullptr && mesh_j->ind != cfd_mesh->ind){
                    glm::vec2 rel_pos = mesh_j->center - cfd_mesh->center;
                    dist = sqrt(dot(rel_pos,rel_pos));
                    mass_j = mesh_j->U[0]*mesh_j->vol;
                    cfd_mesh->fx += G*mass_j/(dist*dist*dist)*rel_pos.x;
                    cfd_mesh->fy += G*mass_j/(dist*dist*dist)*rel_pos.y;
                }
            }
        }
    }
}

void BarnesHutCFDMesh::reset_mesh_tree_states(){
    // reset all particle defined states to 0
    this -> com = glm::vec2(0,0);
    this -> tot_mass = 0;
    // this -> tot_energy = 0;
    // this -> tot_ang_momentum = 0;
    // this -> tot_momentum = 0;

    if (!this->sub_mesh.empty()){
        for (size_t i=0; i<this->sub_mesh.size(); i++){
            this->sub_mesh[i]->reset_mesh_tree_states();
        }
    }
}