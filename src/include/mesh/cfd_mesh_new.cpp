#include <mesh.h>
#include <initialization.h>

// CFDMeshWrapper::CFDMeshWrapper(int n, float l) {
//     float dx = l / n;
//     Vector4f U = Vector4f::Zero();

//     // Resize mesh and initialize cells
//     this->mesh.resize(n * n, nullptr);
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < n; j++) {
//             glm::vec2 center(dx / 2 + dx * j, dx / 2 + dx * i);
//             this->mesh[i * n + j] = new CFDMesh(i * n + j, center, dx, U);
//         }
//     }

//     // Create connections
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < n; j++) {
//             for (int k = 0; k < 4; k++) {
//                 // Setup boundary conditions (adiabatic walls)
//                 if ((i == 0 && k == 3) || (i == n - 1 && k == 1) ||
//                     (j == 0 && k == 2) || (j == n - 1 && k == 0)) {
//                     this->mesh[i * n + j]->neighbor[k] = nullptr;
//                     continue;
//                 }

//                 // Internal connectivity
//                 if (k == 0) {
//                     this->mesh[i * n + j]->neighbor[k] = this->mesh[i * n + j + 1];
//                 } else if (k == 1) {
//                     this->mesh[i * n + j]->neighbor[k] = this->mesh[i * n + j + n];
//                 } else if (k == 2) {
//                     this->mesh[i * n + j]->neighbor[k] = this->mesh[i * n + j - 1];
//                 } else if (k == 3) {
//                     this->mesh[i * n + j]->neighbor[k] = this->mesh[i * n + j - n];
//                 }
//             }
//         }
//     }

//     // Initialize the mesh using a specific initializer
//     CFDInitializer* initializer = new ObliqueShock();
//     initializer->init(this->mesh, n);
// }

// void CFDMeshWrapper::update_mesh() {
//     float dt = 1e-5;
//     for (CFDMesh* mesh_ele : this->mesh) {
//         mesh_ele->computeFlux();
//         mesh_ele->U += mesh_ele->dUdt * dt;

//     }
// }

// CFDMesh::CFDMesh(int ind, glm::vec2 center, float dx, Vector4f U)
//     : center(center), dx(dx), U(U), ind(ind) {
//     this->vol = dx * dx;
//     this->prim = Eigen::Matrix<float, 5, 1>::Zero(5);

//     this->neighbor.resize(4, nullptr);
//     this->face_vals.resize(4, Vector4f::Zero());
//     this->face_vals_prim.resize(4, Eigen::Matrix<float, 5, 1>::Zero());
//     this->flux.resize(4, Vector4f::Zero());

//     if (this->mu != 0) {
//         this->face_gradient.resize(4, Matrix2f::Zero());
//     }
// }

// void CFDMesh::cons2prim() {
//     this->prim(seq(0, 3)) = this->U;
//     this->prim(seq(1, 2)) /= this->prim(0);
//     this->vel_squared = this->prim(seq(1, 2)).squaredNorm();
//     this->prim(3) = this->prim(3) / this->prim(0) - this->vel_squared / 2;

//     // float gamma = (this->R + this->cv) / this->cv;
//     float tmp_p = (gamma - 1) * (this->U[3] - 0.5 * this->U[0] * this->prim(seq(1, 2)).squaredNorm());
//     this->prim(4) = tmp_p;
// }

// void CFDMesh::prim2cons() {
//     this->U = this->prim(seq(0, 3));
//     this->U(seq(1, 2)) *= this->U[0];
//     // float temp = this->prim(4) / this->R / this->U[0];
//     this->U(3) = this->prim[4]/(gamma-1) + (this->prim(seq(1, 2)).squaredNorm() / 2) * this->U[0];
// }


// void CFDMesh::computeFlux() {
//     // compute conv flux for euler equations on dUdt
//     this->dUdt = Vector4f::Zero();
//     for (size_t i = 0; i < this->neighbor.size(); i++) {
//         if (this->neighbor[i] == nullptr) {
//             // for adiabatic no slip wall
//             // rho, p, e is the same
//             // u, v opposite
//             float p = (gamma - 1)*((this->U[3] - (this->U[1]*this->U[1] + this->U[2]*this->U[2])/2/this->U[0]));
//             int dir = -1;
//             if (i >= 2){
//                 dir = 1;
//             }
//             if (i == 0 || i == 2){
//                 this->dUdt[1] = dir*p*dx;
//             }
//             else {
//                 this->dUdt[2] = dir*p*dx;
//             }
//             continue;
//         }
        
//         // catesian grid with square grid points
//         // horizontal handling
//         if (i==0 || i == 2){
//             // direction either -1 for entering or 1 for exiting
//             int dir = i-1;
//             float u = this->U[1]/this->U[0];
//             if (u+this->neighbor[i]->U[1]/this->neighbor[i]->U[i]>0){
//                 // if averaged velocity at face is > 0, this cell is upwind
//                 float p = (gamma - 1)*((this->U[3] - (this->U[1]*this->U[1] + this->U[2]*this->U[2])/2/this->U[0]));
//                 this->dUdt[0] -= dir*this->U[1]*dx;
//                 this->dUdt[1] -= dir*(this->U[1]*u + p)*dx;
//                 this->dUdt[2] -= dir*(this->U[2]*u)*dx;
//                 this->dUdt[3] -= dir*(this->U[3]*u + p*u)*dx;
//             }
//             else {
//                 float p = (gamma - 1)*((this->neighbor[i]->U[3] - (this->neighbor[i]->U[1]*this->neighbor[i]->U[1] + this->neighbor[i]->U[2]*this->neighbor[i]->U[2])/2/this->neighbor[i]->U[0]));
//                 this->dUdt[0] -= dir*this->neighbor[i]->U[1]*dx;
//                 this->dUdt[1] -= dir*(this->neighbor[i]->U[1]*u + p)*dx;
//                 this->dUdt[2] -= dir*(this->neighbor[i]->U[2]*u)*dx;
//                 this->dUdt[3] -= dir*(this->neighbor[i]->U[3]*u + p*u)*dx;
//             }
//         }
//         // vertical handling
//         else {
//             // direction either -1 for entering or 1 for exiting
//             int dir = i-2;
//             float v = this->U[2]/this->U[0];
//             if (v+this->neighbor[i]->U[2]/this->neighbor[i]->U[i]>0){
//                 // if averaged velocity at face is > 0, this cell is upwind
//                 float p = (gamma - 1)*((this->U[3] - (this->U[1]*this->U[1] + this->U[2]*this->U[2])/2/this->U[0]));
//                 this->dUdt[0] -= dir*this->U[1]*dx;
//                 this->dUdt[1] -= dir*(this->U[1]*v)*dx;
//                 this->dUdt[2] -= dir*(this->U[2]*v + p)*dx;
//                 this->dUdt[3] -= dir*(this->U[3]*v + p*v)*dx;
//             }
//             else {
//                 float p = (gamma - 1)*((this->neighbor[i]->U[3] - (this->neighbor[i]->U[1]*this->neighbor[i]->U[1] + this->neighbor[i]->U[2]*this->neighbor[i]->U[2])/2/this->neighbor[i]->U[0]));
//                 this->dUdt[0] -= dir*this->neighbor[i]->U[1]*dx;
//                 this->dUdt[1] -= dir*(this->neighbor[i]->U[1]*v)*dx;
//                 this->dUdt[2] -= dir*(this->neighbor[i]->U[2]*v + p)*dx;
//                 this->dUdt[3] -= dir*(this->neighbor[i]->U[3]*v + p*v)*dx;
//             }
//         }
//     }

//     this->dUdt /= this->vol;
// }

#include <mesh.h>
#include <initialization.h>

CFDMeshWrapper::CFDMeshWrapper(int nx, int ny, int ghost_cell, float dx) {
    // float dx = l / n;
    Vector4f U = Vector4f::Zero();

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
    this->initializer->init(this->mesh, nx, ny, ghost_cell);
    this->ghost_cell = ghost_cell;
    this->nx = nx; this->ny = ny;
    // delete initializer; // Cleanup to prevent memory leaks
}

void CFDMeshWrapper::update_mesh() {
    float dt = 1e-5;
    dt = 0.11;
    // dt /= 10;
    for (CFDMesh* mesh_ele : this->mesh) {
        if (mesh_ele->ind%(this->nx+2*this->ghost_cell)<this->ghost_cell-1 || 
        mesh_ele->ind%(this->nx+2*this->ghost_cell) >= this->nx+this->ghost_cell+1 ||
        mesh_ele->ind<(nx+ghost_cell*2)*(ghost_cell-1) || mesh_ele->ind>=(nx+ghost_cell*2)*(ghost_cell+ny+1))
        {
            continue;
        }
        mesh_ele->minmodSlope();
        mesh_ele->dUdt = Vector4f::Zero();
    }
 
    for (CFDMesh* mesh_ele : this->mesh) {
        // need to compute some of ghost cell fluxes
        if (mesh_ele->ind%(this->nx+2*this->ghost_cell)<this->ghost_cell-1 || 
        mesh_ele->ind%(this->nx+2*this->ghost_cell) >= this->nx+this->ghost_cell ||
        mesh_ele->ind<(nx+ghost_cell*2)*(ghost_cell-1) || mesh_ele->ind>=(nx+ghost_cell*2)*(ghost_cell+ny))
        {
            continue;
        }

        mesh_ele->computeHLLFlux();

        // mesh_ele->computeFluxSecondOrder();
    

    }

    for (CFDMesh* mesh_ele : this->mesh) {
        if (mesh_ele->ind%(this->nx+2*this->ghost_cell)<this->ghost_cell || 
        mesh_ele->ind%(this->nx+2*this->ghost_cell) >= this->nx+this->ghost_cell ||
        mesh_ele->ind<(nx+ghost_cell*2)*ghost_cell || mesh_ele->ind>=(nx+ghost_cell*2)*(ghost_cell+ny))
        {
            continue;
        }
        // mesh_ele->U = mesh_ele->U + mesh_ele->dUdt * dt;
        mesh_ele->U += mesh_ele->dUdt * dt;
        if (mesh_ele->U.array().isNaN().any()) {
            cout << "nan" << endl;
        }
    }

    this->initializer->updateBC(this->mesh,this->nx,this->ny,this->ghost_cell);
}

// void CFDMeshWrapper::update_mesh() {
//     float dt = 1e-5;
//     dt = 0.0011;
//     dt /= 10;

//     // Stage 1: Compute k1
//     for (CFDMesh* mesh_ele : this->mesh) {
//         if (mesh_ele->ind % (this->nx + 2 * this->ghost_cell) < this->ghost_cell - 1 || 
//             mesh_ele->ind % (this->nx + 2 * this->ghost_cell) >= this->nx + this->ghost_cell + 1 ||
//             mesh_ele->ind < (nx + ghost_cell * 2) * (ghost_cell - 1) || 
//             mesh_ele->ind >= (nx + ghost_cell * 2) * (ghost_cell + ny + 1)) {
//             continue;
//         }
//         mesh_ele->minmodSlope();
//     }

//     for (CFDMesh* mesh_ele : this->mesh) {
//         if (mesh_ele->ind % (this->nx + 2 * this->ghost_cell) < this->ghost_cell || 
//             mesh_ele->ind % (this->nx + 2 * this->ghost_cell) >= this->nx + this->ghost_cell ||
//             mesh_ele->ind < (nx + ghost_cell * 2) * ghost_cell || 
//             mesh_ele->ind >= (nx + ghost_cell * 2) * (ghost_cell + ny)) {
//             continue;
//         }
//         mesh_ele->computeHLLFlux();
//     }

//     for (CFDMesh* mesh_ele : this->mesh) {
//         mesh_ele->k1 = mesh_ele->dUdt * dt;
//     }

//     // Apply k1 to get intermediate state
//     for (CFDMesh* mesh_ele : this->mesh) {
//         mesh_ele->U_temp = mesh_ele->U + mesh_ele->k1 * 0.5; // U_temp = U + 0.5 * k1
//     }

//     this->initializer->updateBC(this->mesh, this->nx, this->ny, this->ghost_cell);

//     // Stage 2: Compute k2 using intermediate state
//     for (CFDMesh* mesh_ele : this->mesh) {
//         if (mesh_ele->ind % (this->nx + 2 * this->ghost_cell) < this->ghost_cell - 1 || 
//             mesh_ele->ind % (this->nx + 2 * this->ghost_cell) >= this->nx + this->ghost_cell + 1 ||
//             mesh_ele->ind < (nx + ghost_cell * 2) * (ghost_cell - 1) || 
//             mesh_ele->ind >= (nx + ghost_cell * 2) * (ghost_cell + ny + 1)) {
//             continue;
//         }
//         mesh_ele->minmodSlope();
//     }

//     for (CFDMesh* mesh_ele : this->mesh) {
//         if (mesh_ele->ind % (this->nx + 2 * this->ghost_cell) < this->ghost_cell || 
//             mesh_ele->ind % (this->nx + 2 * this->ghost_cell) >= this->nx + this->ghost_cell ||
//             mesh_ele->ind < (nx + ghost_cell * 2) * ghost_cell || 
//             mesh_ele->ind >= (nx + ghost_cell * 2) * (ghost_cell + ny)) {
//             continue;
//         }
//         mesh_ele->computeHLLFlux();
//     }

//     for (CFDMesh* mesh_ele : this->mesh) {
//         mesh_ele->k2 = mesh_ele->dUdt * dt;
//     }

//     // Final update: Combine k1 and k2
//     for (CFDMesh* mesh_ele : this->mesh) {
//         mesh_ele->U += 0.5 * (mesh_ele->k1 + mesh_ele->k2); // U = U + 0.5 * (k1 + k2)
//         if (mesh_ele->U.array().isNaN().any()) {
//             cout << "nan" << endl;
//         }
//     }

//     this->initializer->updateBC(this->mesh, this->nx, this->ny, this->ghost_cell);
// }

CFDMesh::CFDMesh(int ind, glm::vec2 center, float dx, Vector4f U)
    : center(center), dx(dx), U(U), ind(ind) {
    this->vol = dx * dx;
    this->prim = Eigen::Matrix<float, 5, 1>::Zero();

    this->neighbor.resize(4, nullptr);
    this->face_vals.resize(4, Vector4f::Zero());
    this->face_vals_prim.resize(4, Eigen::Matrix<float, 5, 1>::Zero());
    this->flux.resize(4, Vector4f::Zero());

    // if (this->mu != 0) {
    //     this->face_gradient.resize(4, Matrix2f::Zero());
    // }
}

void CFDMesh::cons2prim() {
    this->prim(seq(0, 3)) = this->U;
    this->prim(seq(1, 2)) /= this->prim(0);
    this->vel_squared = this->prim(seq(1, 2)).squaredNorm();
    this->prim(3) = this->prim(3) / this->prim(0) - this->vel_squared / 2;

    float tmp_p = (gamma - 1) * (this->U[3] - 0.5 * this->U[0] * this->prim(seq(1, 2)).squaredNorm());
    this->prim(4) = tmp_p;
}

void CFDMesh::prim2cons() {
    this->U = this->prim(seq(0, 3));
    this->U(seq(1, 2)) *= this->U[0];
    this->U(3) = this->prim[4] / (gamma - 1) + (this->prim(seq(1, 2)).squaredNorm() / 2) * this->U[0];
}

// void CFDMesh::computeFlux() {
//     // Initialize flux computation
//     this->dUdt = Vector4f::Zero();

//     for (size_t i = 0; i < this->neighbor.size(); i++) {
//         if (this->neighbor[i] == nullptr) {
//             // Adiabatic no-slip wall conditions
//             float p = (gamma - 1) * (this->U[3] - 0.5 * (this->U[1] * this->U[1] + this->U[2] * this->U[2]) / this->U[0]);
//             int dir = (i >= 2) ? 1 : -1;
//             if (i == 0 || i == 2) {
//                 this->dUdt[1] = dir * p * dx;
//             } else {
//                 this->dUdt[2] = dir * p * dx;
//             }
//             continue;
//         }

//         // Handle flux contributions
//         int dir;
//         float vel_avg;
//         if (i == 0 || i == 2) { // Horizontal flux
//             dir = (i == 0) ? 1 : -1;
//             vel_avg = (this->U[1] / this->U[0] + this->neighbor[i]->U[1] / this->neighbor[i]->U[0]) / 2;
//         }
//         // Add similar handling for vertical flux (i == 1 || i == 3)
//         else {
//             dir = (i == 1) ? 1 : -1;
//             vel_avg = (this->U[2] / this->U[0] + this->neighbor[i]->U[2] / this->neighbor[i]->U[0]) / 2;
//         }

//         bool this_upwind;
//         if (vel_avg == 0) {
//             float p1 = (gamma - 1) * (this->U[3] - 0.5 * (this->U[1] * this->U[1] + this->U[2] * this->U[2]) / this->U[0]);
//             float p2 = (gamma - 1) * (this->neighbor[i]->U[3] - 0.5 * (this->neighbor[i]->U[1] * this->neighbor[i]->U[1] + this->neighbor[i]->U[2] * this->neighbor[i]->U[2]) / this->neighbor[i]->U[0]);
//             this->dUdt -= dir * computeFluxForFace((this->U+this->neighbor[i]->U)/2, (p1+p2)/2, vel_avg, dx, i%2);
//             continue;
//         }
//         else if (vel_avg > 0 && i<=1) {
//             // Upwind scheme (current cell is upwind)
//             this_upwind = true;
//         } else if (vel_avg < 0 && i<=1) {
//             // Neighbor is upwind
//             this_upwind = false;
//         }
//         else if (vel_avg < 0 && i>=2) {
//             // Upwind scheme (current cell is upwind)
//             this_upwind = true;
//         } else if (vel_avg > 0 && i>=2) {
//             // Neighbor is upwind
//             this_upwind = false;
//         }

//         if (this_upwind) {
//             float p = (gamma - 1) * (this->U[3] - 0.5 * (this->U[1] * this->U[1] + this->U[2] * this->U[2]) / this->U[0]);
//             this->dUdt -= dir * computeFluxForFace(this->U, p, vel_avg, dx, i%2);
//         }
//         else {
//             float p = (gamma - 1) * (this->neighbor[i]->U[3] - 0.5 * (this->neighbor[i]->U[1] * this->neighbor[i]->U[1] + this->neighbor[i]->U[2] * this->neighbor[i]->U[2]) / this->neighbor[i]->U[0]);
//             this->dUdt -= dir * computeFluxForFace(this->neighbor[i]->U, p, vel_avg, dx, i%2);
//         }
        
//     }

//     this->dUdt /= this->vol;
// }

// // Helper function for flux calculation
// Vector4f CFDMesh::computeFluxForFace(const Vector4f& U, float p, float vel, float dx, int face_orientation) {
//     Vector4f flux;
//     if (face_orientation == 0){
//         flux[0] = U[1] * dx;
//         flux[1] = (U[1] * vel + p) * dx;
//         flux[2] = (U[2] * vel) * dx;
//         flux[3] = (U[3] * vel + p * vel) * dx;
//     }
//     else if (face_orientation == 0){
//         flux[0] = U[1] * dx;
//         flux[1] = (U[1] * vel) * dx;
//         flux[2] = (U[2] * vel + p) * dx;
//         flux[3] = (U[3] * vel + p * vel) * dx;
//     }
//     return flux;
// }

void CFDMesh::computeFlux() {
    // Initialize flux computation
    this->dUdt = Vector4f::Zero();
    bool this_upwind;

    if (this->ind >=9 && this->ind <= 11){
        cout << "eee" << endl;
    }

    for (size_t i = 0; i < this->neighbor.size(); i++) {
        int dir = (i == 0 || i == 1) ? 1 : -1;
        float vel_avg;
        int face_orientation = i % 2; // 0: horizontal, 1: vertical

        if (this->neighbor[i] == nullptr) {
            // Handle boundary conditions (e.g., adiabatic wall)
            float p = (gamma - 1) * (this->U[3] - 0.5 * (this->U[1] * this->U[1] + this->U[2] * this->U[2]) / this->U[0]);
            this->dUdt[1 + face_orientation] -= dir * (this->U[face_orientation+1]*this->U[face_orientation+1] /this->U[0] + p) * dx;
            continue;
        }


        // Compute average velocity and pressure
        vel_avg = 0.5 * (this->U[1 + face_orientation] / this->U[0] +
                         this->neighbor[i]->U[1 + face_orientation] / this->neighbor[i]->U[0]);
        // vel_avg = 0.5 * (this->U[1 + face_orientation]  +
        //                  this->neighbor[i]->U[1 + face_orientation]);
        if (vel_avg == 0){
            float p1 = (gamma - 1) * (this->U[3] - 0.5 * (this->U[1] * this->U[1] + this->U[2] * this->U[2]) / this->U[0]);
            float p2 = (gamma - 1) * (this->neighbor[i]->U[3] - 0.5 * (this->neighbor[i]->U[1] * this->neighbor[i]->U[1] + this->neighbor[i]->U[2] * this->neighbor[i]->U[2]) / this->neighbor[i]->U[0]);
            this->dUdt -= dir * computeFluxForFace((this->U+this->neighbor[i]->U)/2, (p1+p2)/2, vel_avg, dx, face_orientation);
            continue;
        }
        this_upwind = (vel_avg > 0 && i<=1) || (vel_avg<0 && i>=2);


        // this_upwind = (vel_avg < 0 && i<=1) || (vel_avg>0 && i>=2);

        Vector4f flux;
        if (this_upwind) {
            float p = (gamma - 1) * (this->U[3] - 0.5 * (this->U[1] * this->U[1] + this->U[2] * this->U[2]) / this->U[0]);
            flux = computeFluxForFace(this->U, p, vel_avg, dx, face_orientation);
            // flux = computeFluxForFace(this->U, p, this->neighbor[i]->U[1+face_orientation]/this->neighbor[i]->U[0], dx, face_orientation);
            // flux = computeFluxForFace(this->U, p, this->U[1+face_orientation]/this->U[0], dx, face_orientation);
        } else {
            float p = (gamma - 1) * (this->neighbor[i]->U[3] - 0.5 * (this->neighbor[i]->U[1] * this->neighbor[i]->U[1] + this->neighbor[i]->U[2] * this->neighbor[i]->U[2]) / this->neighbor[i]->U[0]);
            flux = computeFluxForFace(this->neighbor[i]->U, p, vel_avg, dx, face_orientation);
            // flux = computeFluxForFace(this->neighbor[i]->U, p, this->U[1+face_orientation]/this->U[0], dx, face_orientation);
            // flux = computeFluxForFace(this->neighbor[i]->U, p, this->neighbor[i]->U[1+face_orientation]/this->neighbor[i]->U[0], dx, face_orientation);
        }

        this->dUdt -= dir * flux;
    }

    this->dUdt /= this->vol; // Normalize by cell volume
}


void CFDMesh::computeFluxCentralDifference() {
    // Initialize the flux derivative term
    this->dUdt = Vector4f::Zero();

    for (size_t i = 0; i < this->neighbor.size(); i++) {
        if (this->neighbor[i] == nullptr) {
            // Apply boundary conditions for adiabatic no-slip walls
            float p = (gamma - 1) * (this->U[3] - 0.5 * (this->U[1] * this->U[1] + this->U[2] * this->U[2]) / this->U[0]);
            int dir = (i >= 2) ? 1 : -1;
            if (i == 0 || i == 2) {
                this->dUdt[1] += dir * p * dx;
            } else {
                this->dUdt[2] += dir * p * dx;
            }
            continue;
        }

        // Central difference flux computation
        Vector4f U_left = this->U;                    // Current cell's state
        Vector4f U_right = this->neighbor[i]->U;      // Neighboring cell's state
        // Vector4f U_face = 0.5 * (U_left + U_right);   // State at the face (central average)
        Vector4f U_face = 0.5 * U_left + 0.5 * U_right;   // State at the face (central average)

        // int face_orientation = i % 2; // 0: horizontal, 1: vertical
        // // Compute average velocity and pressure
        // float vel_avg = 0.5 * (this->U[1 + face_orientation] / this->U[0] +
        //                  this->neighbor[i]->U[1 + face_orientation] / this->neighbor[i]->U[0]);

        // bool this_upwind = (vel_avg > 0 && i<=1) || (vel_avg<0 && i>=2);

        // if (this_upwind){
        //     U_face = 0.6 * U_left + 0.4 * U_right;
        // }
        // else {
        //     U_face = 0.6 * U_left + 0.4 * U_right;;
        // }


        // Compute flux components based on Euler equations
        float rho = U_face[0];
        float u = U_face[1] / rho;
        float v = U_face[2] / rho;
        float energy = U_face[3];
        float pressure = (gamma - 1) * (energy - 0.5 * rho * (u * u + v * v));

        Vector4f flux_face;
        if (i == 0 || i == 2) { // Horizontal flux
            flux_face[0] = rho * u;                    // Mass flux
            flux_face[1] = rho * u * u + pressure;     // Momentum in x
            flux_face[2] = rho * u * v;                // Momentum in y
            flux_face[3] = (energy + pressure) * u;    // Energy flux
        } else { // Vertical flux
            flux_face[0] = rho * v;                    // Mass flux
            flux_face[1] = rho * u * v;                // Momentum in x
            flux_face[2] = rho * v * v + pressure;     // Momentum in y
            flux_face[3] = (energy + pressure) * v;    // Energy flux
        }

        // float a = sqrt(gamma*pressure/rho);
        // // cout << a << endl;
        // float M = (abs(U_face[1+i%2]/rho)/abs(a))*U_face[1+i%2]/abs(U_face[1+i%2]);
        // cout << (abs(U_face[1+i%2]/rho)/abs(a)) << endl;
        // cout << U_face[1+i%2]/abs(U_face[1+i%2]) << endl;
        // cout << M << endl;
        // float phi = (std::min(1.f,abs(M)))*M/abs(M);

        // U_face = 0.5 * (U_left + U_right) - 0.5 * phi * (U_right - U_left);   // State at the face (central average)
        // Vector4f flux_face;
        // if (i == 0 || i == 2) { // Horizontal flux
        //     flux_face[0] = rho * u;                    // Mass flux
        //     flux_face[1] = rho * u * u + pressure;     // Momentum in x
        //     flux_face[2] = rho * u * v;                // Momentum in y
        //     flux_face[3] = (energy + pressure) * u;    // Energy flux
        // } else { // Vertical flux
        //     flux_face[0] = rho * v;                    // Mass flux
        //     flux_face[1] = rho * u * v;                // Momentum in x
        //     flux_face[2] = rho * v * v + pressure;     // Momentum in y
        //     flux_face[3] = (energy + pressure) * v;    // Energy flux
        // }

        // Adjust flux derivative for this cell
        int dir = (i == 0 || i == 1) ? 1 : -1; // Direction: positive for right/bottom, negative for left/top
        this->dUdt -= dir * flux_face * dx;
    }

    // Divide by the volume to get the rate of change of conserved quantities
    this->dUdt /= this->vol;
}

// void CFDMesh::computeFluxCentralDifference() {
//     // Initialize the flux derivative term
//     this->dUdt = Vector4f::Zero();

//     // Helper function for minmod limiter
//     auto minmod = [](float a, float b) -> float {
//         if (a * b > 0) {
//             return std::copysign(std::min(std::fabs(a), std::fabs(b)), a);
//         } else {
//             return 0.0f;
//         }
//     };

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
//         // Vector4f U_upwind;

//         // int face_orientation = i % 2; // 0: horizontal, 1: vertical
//         // // Compute average velocity and pressure
//         // float vel_avg = 0.5 * (this->U[1 + face_orientation] / this->U[0] +
//         //                  this->neighbor[i]->U[1 + face_orientation] / this->neighbor[i]->U[0]);

//         // bool this_upwind = (vel_avg > 0 && i<=1) || (vel_avg<0 && i>=2);

//         // if (this_upwind){
//         //     if (this->neighbor[(i+2)%4] != nullptr){

//         //     U_upwind = this->neighbor[(i+2)%4]->U;
//         //     }
//         //     else {
//         //         U_upwind = 
//         //     }
//         // }

//         // Compute flux limiter for each conserved quantity (rho, rho*u, etc.)
//         Vector4f r; // Slope ratio
//         for (int j = 0; j < 4; ++j) {
//             // float delta_left = U_left[j] - (neighbor_left ? neighbor_left->U[j] : U_left[j]);
//             // float delta_right = (neighbor_right ? neighbor_right->U[j] : U_right[j]) - U_right[j];
//             float delta_left = U_left[j] - (neighbor_left ? neighbor_left->U[j] : U_left[j]);
//             float delta_right = (neighbor_right ? neighbor_right->U[j] : U_right[j]) - U_right[j];
//             r[j] = (delta_right != 0.0f) ? delta_left / delta_right : 0.0f;
//         }

//         // Apply the minmod limiter to each component
//         Vector4f phi; // Flux limiter for each component
//         for (int j = 0; j < 4; ++j) {
//             phi[j] = minmod(1.0f, r[j]); // Compute the limited slope
//         }

//         // Modify U_face with the limiter
//         Vector4f U_face = 0.5 * (U_left + U_right) + 0.5 * phi.cwiseProduct(U_right - U_left);

//         // // // Compute flux limiter between U_left and U_right
//         // // float flux_limiter = minmod(U_right[0] - U_left[0], U_right[1] - U_left[1]);

//         // // State at the face with limiter applied (central average with limiter)
//         // Vector4f U_face = 0.5 * U_left + 0.5 * U_right + flux_limiter;

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

void CFDMesh::computeFlux2() {
    // Initialize flux computation
    this->dUdt = Vector4f::Zero();

    for (size_t i = 0; i < this->neighbor.size(); i++) {
        int dir = (i == 0 || i == 1) ? 1 : -1;
        float vel_avg;
        int face_orientation = i % 2; // 0: horizontal, 1: vertical

        if (this->neighbor[i] == nullptr) {
            // Handle boundary conditions (e.g., adiabatic wall)
            float p = (gamma - 1) * (this->U[3] - 0.5 * (this->U[1] * this->U[1] + this->U[2] * this->U[2]) / this->U[0]);
            this->dUdt[1 + face_orientation] += dir * p * dx;
            continue;
        }

        // Compute average velocity and pressure
        vel_avg = 0.5 * (this->U[1 + face_orientation] / this->U[0] +
                         this->neighbor[i]->U[1 + face_orientation] / this->neighbor[i]->U[0]);

        // Upwind decision
        bool this_upwind = (vel_avg > 0 && i <= 1) || (vel_avg < 0 && i > 1);

        // Compute gradients for second-order linear upwind scheme        
        


        // Compute the state at the face using second-order linear upwind differencing
        Vector4f U_face,grad_U;
        if (this_upwind) {
            if (this->neighbor[(i+2)%4] != nullptr){
                grad_U = (this->U - this->neighbor[(i+2)%4]->U) / dx;
                U_face = this->U + 0.5 * grad_U * dx;
            }
            else {
                U_face = this->U;
            }
        }
        else {
            if (this->neighbor[i]->neighbor[i] != nullptr){

                grad_U = (this->neighbor[i]->U - this->neighbor[i]->neighbor[i]->U) / dx;
                U_face = this->neighbor[i]->U + 0.5 * grad_U * dx;
            }
            else {
                U_face = this->neighbor[i]->U;
            }
        }
                             
        // Compute pressure at the face
        float p_face = (gamma - 1) * (U_face[3] - 0.5 * (U_face[1] * U_face[1] + U_face[2] * U_face[2]) / U_face[0]);

        // Compute the flux at the face
        Vector4f flux = computeFluxForFace(U_face, p_face, vel_avg, dx, face_orientation);

        // Accumulate flux contribution
        this->dUdt -= dir * flux;
    }

    this->dUdt /= this->vol; // Normalize by cell volume
}


Vector4f CFDMesh::computeFluxForFace(const Vector4f& U, float p, float vel, float dx, int face_orientation) {
    Vector4f flux;
    if (face_orientation == 0) { // Horizontal flux
        flux[0] = U[1];
        flux[1] = U[1] * vel + p;
        flux[2] = U[2] * vel;
        flux[3] = U[3] * vel + p * vel;
    } else { // Vertical flux
        flux[0] = U[2];
        flux[1] = U[1] * vel;
        flux[2] = U[2] * vel + p;
        flux[3] = U[3] * vel + p * vel;
    }
    return flux * dx;
}


// Vector4f CFDMesh::fluxFunction(const Vector4f& U) {
//     float rho = U[0];
//     float u = U[1] / rho;
//     float v = U[2] / rho;
//     float p = (gamma - 1) * (U[3] - 0.5 * rho * (u * u + v * v));
//     return Vector4f{
//         rho * u,                      // Mass flux
//         rho * u * u + p,              // Momentum in x
//         rho * u * v,                  // Momentum in y
//         (U[3] + p) * u                // Energy flux
//     };
// }

Vector4f CFDMesh::fluxFunction(const Vector4f& U, int face_orientation) {
    float rho = U[0];
    float u = U[1] / rho;
    float v = U[2] / rho;
    float p = (gamma - 1) * (U[3] - 0.5 * rho * (u * u + v * v));

    // Flux calculation depends on face orientation
    if (face_orientation == 0) {
        // x-face flux
        return Vector4f{
            rho * u,               // Mass flux
            rho * u * u + p,       // Momentum flux in x
            rho * u * v,           // Momentum flux in y
            (U[3] + p) * u         // Energy flux
        };
    } else if (face_orientation == 1) {
        // y-face flux
        return Vector4f{
            rho * v,               // Mass flux
            rho * u * v,           // Momentum flux in x
            rho * v * v + p,       // Momentum flux in y
            (U[3] + p) * v         // Energy flux
        };
    } else {
        throw std::invalid_argument("Invalid face orientation. Must be 0 (x-face) or 1 (y-face).");
    }
}




// void CFDMesh::computeFluxMUSCL() {
//     this->dUdt = Vector4f::Zero();

//     for (size_t i = 0; i < this->neighbor.size(); i++) {
//         if (this->neighbor[i] == nullptr) {
//             // Boundary conditions for adiabatic no-slip walls
//             float p = (gamma - 1) * (this->U[3] - 0.5 * (this->U[1] * this->U[1] + this->U[2] * this->U[2]) / this->U[0]);
//             int dir = (i >= 2) ? 1 : -1;
//             if (i == 0 || i == 2) {
//                 this->dUdt[1] += dir * p * dx;
//             } else {
//                 this->dUdt[2] += dir * p * dx;
//             }
//             continue;
//         }

//         // MUSCL reconstruction for face values
//         Vector4f U_left, U_right;
//         if (i<=1){
//             U_left = this->U;
//             U_right = this->neighbor[i]->U;
//         }
//         else {
//             U_left = this->neighbor[i]->U;
//             U_right = this->U;
//         }

//         // Vector4f delta_left = U_right-U_left;
//         // Vector4f delta_right = U_right;

//         // Slope calculation (van Leer limiter)
//         // Vector4f delta_left = U_left - this->neighbor[mod(i - 1, 4)]->U;
//         Vector4f delta_left = U_left - this->neighbor[(i+2)%4]->U;
//         // Vector4f delta_right = this->neighbor[i]->U - U_left;
//         Vector4f delta_center = (U_right - this->neighbor[(i+2)%4]->U)/2;
//         Vector4f delta_right = U_right - U_left;
//         // Vector4f slope = 0.5 * (delta_left + delta_right);
//         Vector4f slope;

//         for (int j = 0; j < 4; ++j) {
//             if (delta_left[j] * delta_right[j] <= 0) {
//                 slope[j] = 0;
//             } else {
//                 slope[j] = std::min(std::abs(delta_left[j]), std::abs(delta_right[j])) * (delta_left[j] > 0 ? 1 : -1);
//                 // slope[j] = std::min(std::abs(delta_left[j]), std::abs(delta_right[j]));
//             }
//             // if (delta_left[j]<0 && delta_right[j]<0 && delta_center[j]<0) {
//             //     slope[j] = 0;
//             // }
//             // else if (delta_left[j]>0 && delta_right[j]>0 && delta_center[j]>0) {
//             //     slope[j] = 0;
//             // }
//             // else {
//             //     slope[j] = std::min(std::min(std::abs(delta_left[j]), std::abs(delta_right[j])), std::abs(delta_left[j])) * (delta_left[j] > 0 ? 1 : -1);
//             // }
//         }

//         // Left and right face states using MUSCL
//         Vector4f U_face_left = U_left + 0.5 * slope;
//         Vector4f U_face_right = U_right - 0.5 * slope;

//         // Determine face orientation (0 for x-face, 1 for y-face)
//         int face_orientation = (i == 0 || i == 2) ? 0 : 1;

//         // Compute flux at the face using a Riemann solver (HLL)
//         Vector4f flux_face = computeHLLFlux(U_face_left, U_face_right, face_orientation);

//         // Adjust flux derivative
//         int dir = (i == 0 || i == 1) ? 1 : -1;
//         this->dUdt += dir * flux_face * dx;
//         // this->dUdt -= -flux_face * dx;
//     }

//     // Normalize by the volume
//     this->dUdt /= this->vol;
// }

// void CFDMesh::computeFluxMUSCL() {
//     this->dUdt = Vector4f::Zero();

//     for (size_t i = 0; i < this->neighbor.size(); i++) {
//         if (this->neighbor[i] == nullptr) {
//             // Boundary conditions for adiabatic no-slip walls
//             float p = (gamma - 1) * (this->U[3] - 0.5 * (this->U[1] * this->U[1] + this->U[2] * this->U[2]) / this->U[0]);
//             int dir = (i == 2 || i == 3) ? 1 : -1;
//             if (i == 0 || i == 2) {
//                 this->dUdt[1] += dir * p * dx;
//             } else {
//                 this->dUdt[2] += dir * p * dx;
//             }
//             continue;
//         }

//         // MUSCL reconstruction for face values
//         Vector4f U_left, U_right;
//         if (i == 0) { 
//             // Right face: current cell (U) on the left, neighbor (right) on the right
//             U_left = this->U;
//             U_right = this->neighbor[i]->U;
//         } else if (i == 1) { 
//             // Top face: current cell (U) below, neighbor (top) above
//             U_left = this->U;
//             U_right = this->neighbor[i]->U;
//         } else if (i == 2) { 
//             // Left face: neighbor (left) on the left, current cell (U) on the right
//             U_left = this->neighbor[i]->U;
//             U_right = this->U;
//         } else if (i == 3) { 
//             // Bottom face: neighbor (bottom) below, current cell (U) above
//             U_left = this->neighbor[i]->U;
//             U_right = this->U;
//         }

//         // Slope calculation (van Leer limiter)
//         Vector4f delta_left = U_left - this->neighbor[(i+2)%4]->U;
//         Vector4f delta_right = this->neighbor[i]->U - U_left;
//         Vector4f slope = 0.5 * (delta_left + delta_right);

//         for (int j = 0; j < 4; ++j) {
//             if (delta_left[j] * delta_right[j] <= 0) {
//                 slope[j] = 0;
//             } else {
//                 slope[j] = std::min(std::abs(delta_left[j]), std::abs(delta_right[j])) * (delta_left[j] > 0 ? 1 : -1);
//             }
//         }

//         // Left and right face states using MUSCL
//         Vector4f U_face_left = U_left + 0.5 * slope;
//         Vector4f U_face_right = U_right - 0.5 * slope;

//         // Determine face orientation (0 for x-face, 1 for y-face)
//         int face_orientation = (i == 0 || i == 2) ? 0 : 1;

//         // Compute flux at the face using a Riemann solver (HLL)
//         Vector4f flux_face = computeHLLFlux(U_face_left, U_face_right, face_orientation);

//         // Adjust flux derivative
//         int dir = (i == 0 || i == 1) ? 1 : -1;
//         this->dUdt += dir * flux_face * dx;
//     }

//     // Normalize by the volume
//     this->dUdt /= this->vol;
// }

// Find minimum of three numbers
float mins(float a, float b, float c){
    float k;
    k = a;
    if (b<k){
        k=b;
    }
    if(c<k){
        k=c;
    }
    return k;
}

// Sign(x) function
int sgn(float v) {
    if (v < 0) return -1;
    if (v > 0) return 1;
    return 0;
}

// Minmod Limiter
float CFDMesh::minmod(float a, float b, float c) {
    float vs[] = {abs(a), abs(b), abs(c)};
    float v = mins(vs[0],vs[1],vs[2]);
    float mm;
    float s;
    // Using Harten's generalized definition
    // minmod: zero if any opposing sign, otherwise the one of smaller magnitude.
    s = (sgn(a) + sgn(b) + sgn(c))/3;

    if (abs(s) == 1.0){
        mm = s * v;
    }

    if (abs(s) != 1){
        mm = 0;
    }
    return mm;
}

void CFDMesh::minmodSlope(){
    Vector4f slope_L, slope_R, slope_C;

    slope_L = 2*(this->U - this->neighbor[0]->U)/this->dx;
    slope_R = 2*(this->neighbor[2]->U - this->U)/this->dx;
    slope_C = (this->neighbor[2]->U - this->neighbor[0]->U)/(2*this->dx);

    for (int i=0;i<4;i++){
        this->slope_x[i] = minmod(slope_L[i],slope_R[i],slope_C[i]);
    }

    slope_L = 2*(this->U - this->neighbor[1]->U)/this->dx;
    slope_R = 2*(this->neighbor[3]->U - this->U)/this->dx;
    slope_C = (this->neighbor[3]->U - this->neighbor[1]->U)/(2*this->dx);

    for (int i=0;i<4;i++){
        this->slope_y[i] = minmod(slope_L[i],slope_R[i],slope_C[i]);
    }
}

void CFDMesh::computeHLLFlux() {
    Vector4f U_L,U_R;
    int normal_x, normal_y;
    for (size_t i=0; i<this->neighbor.size(); i++){
        if (this->neighbor[i]->ind > this->ind){
            if (i==0){
                U_L = this->U + this->slope_x * this->dx/2;
                U_R = this->neighbor[i]->U - this->neighbor[i]->slope_x * this->neighbor[i]->dx/2;
                normal_x = 1;
                normal_y = 0;
            }
            else if (i==2){
                U_L = this->neighbor[i]->U + this->neighbor[i]->slope_x * this->neighbor[i]->dx/2;
                U_R = this->U - this->slope_x * this->dx/2;
                normal_x = -1;
                normal_y = 0;
            }
            else if (i==1){
                U_L = this->U + this->slope_y * this->dx/2;
                U_R = this->neighbor[i]->U - this->neighbor[i]->slope_y * this->neighbor[i]->dx/2;
                normal_x = 0;
                normal_y = 1;
            }
            else if (i==3){
                U_L = this->neighbor[i]->U + this->neighbor[i]->slope_y * this->neighbor[i]->dx/2;
                U_R = this->U - this->slope_y * this->dx/2;
                normal_x = 0;
                normal_y = -1;
            }

            Vector4f flux;

            // Vector4f& flux, const Vector4f& qL, const Vector4f& qR, int normal_x, int normal_y
            this->HLLE(flux,U_L,U_R,normal_x,normal_y);

            this->dUdt -= flux*this->dx;
            this->neighbor[i]->dUdt += flux*this->dx;
        }
    }
    this->dUdt /= this->vol;

    // // Compute primitive variables
    // float rho_L = U_left[0], rho_R = U_right[0];
    // float u_L = U_left[1] / rho_L, u_R = U_right[1] / rho_R;
    // float v_L = U_left[2] / rho_L, v_R = U_right[2] / rho_R;
    // float p_L = (gamma - 1) * (U_left[3] - 0.5 * rho_L * (u_L * u_L + v_L * v_L));
    // float p_R = (gamma - 1) * (U_right[3] - 0.5 * rho_R * (u_R * u_R + v_R * v_R));

    // float u_face = (face_orientation == 0) ? u_L : v_L; // Velocity normal to face
    // float c_L = sqrt(gamma * p_L / rho_L), c_R = sqrt(gamma * p_R / rho_R);

    // // Wave speeds
    // float S_L = std::min(u_face - c_L, u_face - c_R);
    // float S_R = std::max(u_face + c_L, u_face + c_R);

    // // HLL flux calculation
    // if (S_L >= 0) {
    //     return fluxFunction(U_left, face_orientation);
    // } else if (S_R <= 0) {
    //     return fluxFunction(U_right, face_orientation);
    // } else {
    //     Vector4f F_L = fluxFunction(U_left, face_orientation);
    //     Vector4f F_R = fluxFunction(U_right, face_orientation);
    //     return (S_R * F_L - S_L * F_R + S_L * S_R * (U_right - U_left)) / (S_R - S_L);
    // }
}


// void CFDMesh::computeFluxMUSCL() {
//     // Initialize the flux derivative term
//     this->dUdt = Vector4f::Zero();

//     for (size_t i = 0; i < this->neighbor.size(); i++) {
//         if (this->neighbor[i] == nullptr) {
//             // Boundary conditions for adiabatic no-slip walls
//             float p = (gamma - 1) * (this->U[3] - 0.5 * (this->U[1] * this->U[1] + this->U[2] * this->U[2]) / this->U[0]);
//             int dir = (i >= 2) ? 1 : -1;
//             if (i == 0 || i == 2) {
//                 this->dUdt[1] += dir * p * dx;
//             } else {
//                 this->dUdt[2] += dir * p * dx;
//             }
//             continue;
//         }

//         // MUSCL reconstruction for face values
//         Vector4f U_left = this->U; // Current cell state
//         Vector4f U_right = this->neighbor[i]->U; // Neighboring cell state

//         // Slope calculation (van Leer limiter)
//         // Vector4f delta_left = U_left - this->neighbor[mod(i - 1, 4)]->U;  // Left-side slope
//         // Vector4f delta_left = U_left - this->neighbor[(i-1)%4]->U;  // Left-side slope
//         // Vector4f delta_left = U_left - this->neighbor[(i-2)%4]->U;  // Left-side slope
//         // Vector4f delta_right = this->neighbor[i]->U - U_left;  // Right-side slope
//         // Vector4f slope = 0.5 * (delta_left + delta_right);
//         Vector4f delta_left = U_left - this->neighbor[(i+2)%4]->U;  // Left-side slope
//         Vector4f delta_right = U_right - U_left;  // Right-side slope
//         Vector4f slope = 0.5 * (delta_left + delta_right);

//         // Apply the slope limiter
//         for (int j = 0; j < 4; ++j) {
//             if (delta_left[j] * delta_right[j] <= 0) {
//                 slope[j] = 0;  // Prevent oscillations
//             } else {
//                 slope[j] = std::min(std::abs(delta_left[j]), std::abs(delta_right[j])) * (delta_left[j] > 0 ? 1 : -1);
//             }
//         }

//         // Left and right face states using MUSCL
//         Vector4f U_face_left = U_left + 0.5 * slope;
//         Vector4f U_face_right = U_right - 0.5 * slope;

//         // Compute flux at the face using a Riemann solver (HLL)
//         Vector4f flux_face = computeHLLFlux(U_face_left, U_face_right);

//         // Adjust flux derivative
//         int dir = (i == 0 || i == 1) ? 1 : -1; // Direction for face contributions
//         this->dUdt += dir * flux_face * dx;
//     }

//     // Normalize by the volume
//     this->dUdt /= this->vol;
// }

// Vector4f CFDMesh::computeHLLFlux(const Vector4f& U_left, const Vector4f& U_right) {
//     // Compute primitive variables
//     float rho_L = U_left[0], rho_R = U_right[0];
//     float u_L = U_left[1] / rho_L, u_R = U_right[1] / rho_R;
//     float v_L = U_left[2] / rho_L, v_R = U_right[2] / rho_R;
//     float p_L = (gamma - 1) * (U_left[3] - 0.5 * rho_L * (u_L * u_L + v_L * v_L));
//     float p_R = (gamma - 1) * (U_right[3] - 0.5 * rho_R * (u_R * u_R + v_R * v_R));

//     // Speed of sound
//     float c_L = sqrt(gamma * p_L / rho_L);
//     float c_R = sqrt(gamma * p_R / rho_R);

//     // Wave speeds
//     float S_L = std::min(u_L - c_L, u_R - c_R);
//     float S_R = std::max(u_L + c_L, u_R + c_R);

//     // HLL flux calculation
//     if (S_L >= 0) {
//         return fluxFunction(U_left); // Use flux from left state
//     } else if (S_R <= 0) {
//         return fluxFunction(U_right); // Use flux from right state
//     } else {
//         // Intermediate HLL flux
//         Vector4f F_L = fluxFunction(U_left);
//         Vector4f F_R = fluxFunction(U_right);
//         return (S_R * F_L - S_L * F_R + S_L * S_R * (U_right - U_left)) / (S_R - S_L);
//     }
// }

void CFDMesh::computeFluxHLLE() {
    // this->dUdt = Vector4f::Zero();

    // Temporary storage for the fluxes
    Vector4f flux;

    for (size_t i = 0; i < this->neighbor.size(); i++) {
        // if (this->neighbor[i] == nullptr) {
        //     // Boundary condition: Adiabatic no-slip walls
        //     Vector4f qL = this->U;
        //     Vector4f qR = qL;  // Reflective boundary

        //     int normal_x = (i == 0 || i == 2) ? ((i == 0) ? 1 : -1) : 0;
        //     int normal_y = (i == 1 || i == 3) ? ((i == 1) ? 1 : -1) : 0;

        //     HLLE(flux, qL, qR, normal_x, normal_y);
        //     this->dUdt -= flux / this->vol;
        //     continue;
        // }

        // Calculate slopes using Minmod limiter
        Vector4f dqL, dqR;
        for (int k = 0; k < 4; k++) {
            float dqw, dqe, dqc;

            // x-direction slopes
            if (i == 0 || i == 2) {
                dqw = 2 * (this->U[k] - this->neighbor[(i + 2) % 4]->U[k]) / dx;
                dqe = 2 * (this->neighbor[i]->U[k] - this->U[k]) / dx;
                dqc = (this->neighbor[i]->U[k] - this->neighbor[(i + 2) % 4]->U[k]) / dx;
                dqL[k] = minmod(dqw, dqe, dqc);
            }

            // y-direction slopes
            if (i == 1 || i == 3) {
                dqw = 2 * (this->U[k] - this->neighbor[(i + 2) % 4]->U[k]) / dx;
                dqe = 2 * (this->neighbor[i]->U[k] - this->U[k]) / dx;
                dqc = (this->neighbor[i]->U[k] - this->neighbor[(i + 2) % 4]->U[k]) / dx;
                dqL[k] = minmod(dqw, dqe, dqc);
            }
        }

        // Construct left and right states
        Vector4f qL = this->U + 0.5 * dqL;
        // Vector4f qR = this->neighbor[i]->U - 0.5 * dqR;
        Vector4f qR = this->neighbor[i]->U - 0.5 * dqL;

        // Determine the normal direction
        float normal_x = (i == 0 || i == 2) ? ((i == 0) ? dx : -dx) : 0;
        float normal_y = (i == 1 || i == 3) ? ((i == 1) ? dx : -dx) : 0;

        // Compute flux using HLLE
        HLLE(flux, qL, qR, normal_x, normal_y);
        this->dUdt -= flux / this->vol;
        if (flux.array().isNaN().any()) {
            cout << "nan" << endl;
        }
    }
}

// float CFDMesh::minmod(float a, float b, float c) {
//     float s = (std::copysign(1.0, a) + std::copysign(1.0, b) + std::copysign(1.0, c)) / 3.0;
//     if (std::fabs(s) == 1.0) {
//         return s * std::min({std::fabs(a), std::fabs(b), std::fabs(c)});
//     }
//     return 0.0;
// }


void CFDMesh::HLLE(Vector4f& flux, const Vector4f& qL, const Vector4f& qR, int nx, int ny) {
    // Normal vectors

    // Extract left state variables
    float rL = qL[0];
    float uL = qL[1] / rL;
    float vL = qL[2] / rL;
    float vnL = uL * nx + vL * ny;
    float pL = (gamma - 1) * (qL[3] - rL * (uL * uL + vL * vL) / 2);
    float aL = std::sqrt(gamma * pL / rL);
    float HL = (qL[3] + pL) / rL;

    // Extract right state variables
    float rR = qR[0];
    float uR = qR[1] / rR;
    float vR = qR[2] / rR;
    float vnR = uR * nx + vR * ny;
    float pR = (gamma - 1) * (qR[3] - rR * (uR * uR + vR * vR) / 2);
    float aR = std::sqrt(gamma * pR / rR);
    float HR = (qR[3] + pR) / rR;

    // Compute Roe averages
    float RT = std::sqrt(rR / rL);
    float u = (uL + RT * uR) / (1 + RT);
    float v = (vL + RT * vR) / (1 + RT);
    float H = (HL + RT * HR) / (1 + RT);
    float a = std::sqrt((gamma - 1) * (H - (u * u + v * v) / 2));
    float vn = u * nx + v * ny;

    // Wave speed estimates
    float SLm = std::min({vnL - aL, vn - a, 0.0f});
    float SRp = std::max({vnR + aR, vn + a, 0.0f});

    if (SLm >= SRp) {
        std::cerr << "Warning: Degenerate wave speeds detected (SLm >= SRp)" << std::endl;
        flux.setZero();
        return;
    }

    // Left and right fluxes
    Vector4f FL, FR;
    FL << rL * vnL,
          rL * vnL * uL + pL * nx,
          rL * vnL * vL + pL * ny,
          rL * vnL * HL;
    FR << rR * vnR,
          rR * vnR * uR + pR * nx,
          rR * vnR * vR + pR * ny,
          rR * vnR * HR;

    // int sign_change = 1;
    // if (nx==-1 || ny==-1){
    //     sign_change = -1;
    // }

    // Compute HLLE flux
    // flux = (SRp * FL - SLm * FR + sign_change*SLm * SRp * (qR - qL)) / (SRp - SLm);
    flux = (SRp * FL - SLm * FR + SLm * SRp * (qR - qL)) / (SRp - SLm);
    if (flux.array().isNaN().any()) {
        cout << "nan" << endl;
    }
    // this->dUdt -= flux;
}


// Reconstruct left and right states at the interface
void CFDMesh::reconstructStates(Vector4f& U_L, Vector4f& U_R, const Vector4f& grad) const {
    U_L = this->U + 0.5f * grad * dx;
    U_R = this->U - 0.5f * grad * dx;
}

// Placeholder for the Riemann solver (e.g., Roe's flux)
// Vector4f CFDMesh::computeRiemannFlux(const Vector4f& U_L, const Vector4f& U_R, int face_orientation) const {
//     // Implement Roe's flux or another appropriate Riemann solver here
//     // This is a placeholder implementation
//     Vector4f flux = Vector4f::Zero();
//     // ... Riemann solver logic ...
//     return flux;
// }

// Vector4f CFDMesh::computeRiemannFlux(const Vector4f& U_L, const Vector4f& U_R, int face_orientation) const {
//     // Constants
//     // const float gamma = 1.4f; // Specific heat ratio for ideal gas

//     // Extract primitive variables from conserved variables
//     auto conservedToPrimitive = [](const Vector4f& U) -> Vector4f {
//         float rho = U[0];
//         float u = U[1] / rho;
//         float v = U[2] / rho;
//         float p = (1.4 - 1) * (U[3] - 0.5f * (U[1] * U[1] + U[2] * U[2]) / rho);
//         return Vector4f(rho, u, v, p);
//     };

//     Vector4f prim_L = conservedToPrimitive(U_L);
//     Vector4f prim_R = conservedToPrimitive(U_R);

//     float rho_L = prim_L[0], u_L = prim_L[1], v_L = prim_L[2], p_L = prim_L[3];
//     float rho_R = prim_R[0], u_R = prim_R[1], v_R = prim_R[2], p_R = prim_R[3];

//     // Compute normal velocity and direction
//     float u_L_norm = (face_orientation == 0) ? u_L : v_L;
//     float u_R_norm = (face_orientation == 0) ? u_R : v_R;

//     // Compute Roe averages
//     float sqrt_rho_L = std::sqrt(rho_L);
//     float sqrt_rho_R = std::sqrt(rho_R);
//     float rho_avg = sqrt_rho_L * sqrt_rho_R;
//     float u_avg = (sqrt_rho_L * u_L_norm + sqrt_rho_R * u_R_norm) / (sqrt_rho_L + sqrt_rho_R);
//     float H_L = (U_L[3] + p_L) / rho_L; // Enthalpy
//     float H_R = (U_R[3] + p_R) / rho_R; // Enthalpy
//     float H_avg = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) / (sqrt_rho_L + sqrt_rho_R);
//     float c_avg = std::sqrt((gamma - 1) * (H_avg - 0.5f * u_avg * u_avg)); // Speed of sound

//     // Compute wave speeds
//     float lambda1 = u_avg - c_avg;
//     float lambda2 = u_avg;
//     float lambda3 = u_avg + c_avg;

//     // Roe entropy fix
//     auto entropyFix = [](float lambda, float epsilon = 1e-2f) -> float {
//         return (std::abs(lambda) > epsilon) ? lambda : 0.5f * (lambda * lambda / epsilon + epsilon);
//     };

//     lambda1 = entropyFix(lambda1);
//     lambda2 = entropyFix(lambda2);
//     lambda3 = entropyFix(lambda3);

//     // Compute fluxes
//     Vector4f F_L, F_R;
//     F_L[0] = rho_L * u_L_norm;
//     F_L[1] = rho_L * u_L_norm * u_L + p_L * (face_orientation == 0);
//     F_L[2] = rho_L * u_L_norm * v_L + p_L * (face_orientation == 1);
//     F_L[3] = u_L_norm * (U_L[3] + p_L);

//     F_R[0] = rho_R * u_R_norm;
//     F_R[1] = rho_R * u_R_norm * u_R + p_R * (face_orientation == 0);
//     F_R[2] = rho_R * u_R_norm * v_R + p_R * (face_orientation == 1);
//     F_R[3] = u_R_norm * (U_R[3] + p_R);

//     Vector4f delta_U = U_R - U_L;

//     // Compute Roe flux
//     Vector4f flux = 0.5f * (F_L + F_R) -
//                     0.5f * (std::abs(lambda1) * delta_U +
//                             std::abs(lambda2) * delta_U +
//                             std::abs(lambda3) * delta_U);

//     return flux;
// }


// // Second-order upwind flux computation
// void CFDMesh::computeFluxSecondOrder() {
//     // Initialize flux computation
//     this->dUdt = Vector4f::Zero();

//     for (size_t i = 0; i < this->neighbor.size(); i++) {
//         int dir = (i == 0 || i == 1) ? 1 : -1;
//         int face_orientation = i % 2; // 0: horizontal, 1: vertical

//         // Handle boundary conditions
//         if (this->neighbor[i] == nullptr) {
//             // For second-order, use ghost cells or appropriate boundary reconstruction
//             // Here, we assume reflective (adiabatic wall) boundary conditions
//             Vector4f U_reflected = this->U;
//             U_reflected[1 + face_orientation] = -U_reflected[1 + face_orientation]; // Reverse velocity component

//             // Compute pressure
//             float p = (gamma - 1.0f) * (U_reflected[3] - 0.5f * 
//                         (U_reflected[1] * U_reflected[1] + U_reflected[2] * U_reflected[2]) / U_reflected[0]);

//             // Compute flux for boundary face
//             Vector4f flux = computeFluxForFace(U_reflected, p, 
//                                 U_reflected[1 + face_orientation] / U_reflected[0], dx, face_orientation);

//             this->dUdt[1 + face_orientation] += dir * flux[1 + face_orientation];
//             continue;
//         }

//         // Compute gradients using neighboring cells
//         Vector4f grad = computeGradient(i); // Implemented to compute gradient with limiter

//         // Reconstruct left and right states at the interface
//         Vector4f U_L, U_R;
//         reconstructStates(U_L, U_R, grad);

//         // Similarly, compute gradients for the neighbor cell for its reconstruction
//         Vector4f neighbor_grad = this->neighbor[i]->computeGradient(i);
//         Vector4f neighbor_U_L, neighbor_U_R;
//         this->neighbor[i]->reconstructStates(neighbor_U_L, neighbor_U_R, neighbor_grad);

//         // Choose appropriate states based on velocity direction
//         float vel_avg = 0.5f * (this->U[1 + face_orientation] / this->U[0] +
//                                 this->neighbor[i]->U[1 + face_orientation] / this->neighbor[i]->U[0]);

//         Vector4f flux;
//         if (vel_avg > 0) {
//             flux = computeRiemannFlux(U_L, this->neighbor[i]->U_R, face_orientation);
//         } else {
//             flux = computeRiemannFlux(this->neighbor[i]->U_L, U_R, face_orientation);
//         }

//         // Update dUdt
//         this->dUdt -= dir * flux;
//     }

//     // Normalize by cell volume
//     this->dUdt /= this->vol;
// }
