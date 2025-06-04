#include <mesh.h>
#include <initialization.h>

CFDMeshWrapper::CFDMeshWrapper(int n,float l)
{
    // this->constructMatrices();
    float dx = l/n;
    
    Vector4f U = Vector4f::Zero();

    this->mesh.resize(n*n,nullptr);
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            glm::vec2 center(dx/2+dx*j,dx/2+dx*i);
            this->mesh[i*n+j] = new CFDMesh(i*n+j,center,dx,U);
        }
    }

    // create connections
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            for (int k=0;k<4;k++){

                // setup bcs adiabatic wall
                if (i==0 && k==3){
                    this->mesh[i*n+j]->neighbor[k] = nullptr;
                    continue;
                }
                else if (i==n-1 && k==1){
                    this->mesh[i*n+j]->neighbor[k] = nullptr;
                    continue;
                }
                else if (j==0 && k==2){
                    this->mesh[i*n+j]->neighbor[k] = nullptr;
                    continue;
                }
                else if (j==n-1 && k==0){
                    this->mesh[i*n+j]->neighbor[k] = nullptr;
                    continue;
                }

                // internal connectivity
                if (k==0){
                    this->mesh[i*n+j]->neighbor[k] = this->mesh[i*n+j+1];
                }
                else if (k==1){
                    this->mesh[i*n+j]->neighbor[k] = this->mesh[i*n+j+n];
                }
                else if (k==2){
                    this->mesh[i*n+j]->neighbor[k] = this->mesh[i*n+j-1];
                }
                else if (k==3){
                    this->mesh[i*n+j]->neighbor[k] = this->mesh[i*n+j-n];
                }

            }
        }
    }

    // CFDInitializer* initializer = new KelvinHelholtz();
    // initializer->init(this->mesh,n);

    // CFDInitializer* initializer = new ObliqueShock();
    // initializer->init(this->mesh,n);

    
    CFDInitializer* initializer = new GravPressureFull();
    initializer->init(this->mesh,n);


}



void CFDMeshWrapper::update_mesh(){
    // run one iteration
    // temporarily use euler
    // float dt = 6e-4;
    float max_temp = -1e3;
    float max_vel = -1e-3;
    for (CFDMesh* mesh_ele: this->mesh){
        mesh_ele->cons2prim();
        max_temp = std::max(mesh_ele->temp,max_temp);
        max_vel = std::max(mesh_ele->vel_squared,max_vel);
    }
    max_vel = sqrt(max_vel);
    float cp = this->mesh[0]->R+this->mesh[0]->cv;
    float gamma = cp/this->mesh[0]->cv;
    // float dt = 1e-4;
    float dt = 0.45*this->mesh[0]->dx/(max_vel+sqrt(gamma*this->mesh[0]->R*max_temp));
    float cfl = (max_vel+sqrt(gamma*this->mesh[0]->R*max_temp))*dt/this->mesh[0]->dx;
    // std::cout << cfl << std::endl;
    // std::cout << dt << std::endl;
    dt = 1e-3;
    // if (cfl<0.5){
    //     std::cout << cfl << std::endl;
    // }
    for (CFDMesh* mesh_ele: this->mesh){
        if (mesh_ele->mu!=0){
            mesh_ele->computeGradient();
        }
        mesh_ele->computeFaceVals();
    }
    for (CFDMesh* mesh_ele: this->mesh){
        mesh_ele->computeFlux();
        mesh_ele->computeSource();
        mesh_ele->computeDE();
        if (mesh_ele->U[0] + mesh_ele->dUdt[0]*dt < 0){
            cout<<"neg dens" << endl;
        }

        mesh_ele->U += mesh_ele->dUdt*dt;
        // if (mesh_ele->U[0] < 0){
        //     mesh_ele->U[0] = 1e-12;
        // }
    }
}

// void CFDMeshWrapper::constructMatrices(){
//     // vector = [rho;rho_u;rho_v;t;p]

//     // construct vector of values using dereferenced pointers

//     // construct vector of derivitives

//     // construct matrix of cfd operations
// }

CFDMesh::CFDMesh(int ind, glm::vec2 center, float dx, Vector4f U):
center(center),dx(dx),U(U),ind(ind)
{
    this->vol = dx*dx;
    this->prim = Eigen::Matrix<float, 5, 1>::Zero(5);

    this->neighbor.resize(4,nullptr);
    this->face_vals.resize(4,Vector4f::Zero());

    // this->face_vals_prim.resize(4,VectorXf::Zero(5));
    this->face_vals_prim.resize(4,Eigen::Matrix<float, 5, 1>::Zero());
    this->flux.resize(4,Vector4f::Zero());
    if (this->mu!=0){
        this->face_gradient.resize(4,Matrix2f::Zero());
    }
}

void CFDMesh::cons2prim(){
    // this->prim = this->U;
    // this->prim = Eigen::Matrix<float, 5, 1>::Zero(5);
    this->prim(seq(0,3)) = this->U;
    // recover velocities
    this->prim(seq(1,2)) /= this->prim(0);
    // calc vel mag
    this->vel_squared = this->prim(seq(1,2)).squaredNorm();
    // recovre energy
    this->prim(3) = this->prim(3)/this->prim(0)-this->vel_squared/2;
    // recover pressure
    // this->p = this->prim(0)*this->prim(3);

    this->temp =this->prim(3)/this->cv;
    float gamma = (this->R+this->cv)/this->cv;
    float tmp_p = (gamma-1) * (this->U[3]-0.5*this->U[0]*(this->prim(seq(1,2)).squaredNorm()));
    // this->prim(4) = this->R*this->prim(0)*this->temp;
    this->prim(4) = tmp_p;
    if (tmp_p != this->prim(4)){
        cout << "p mismatch" << endl;
    }
}

void CFDMesh::prim2cons(){
    // this->prim = this->U;
    this->U = this->prim(seq(0,3));
    this->U(seq(1,2)) *= this->U[0];
    float temp = this->prim(4)/this->R/this->U[0];
    this->U(3) = (this->cv*temp+this->prim(seq(1,2)).squaredNorm()/2)*this->U[0];
}


void CFDMesh::computeFaceVals(){
    // center difference method to compute face value
    // 2nd order but may induce oscillations

    // need to implement for weighted average if differing size meshes
    Vector4f face_val;
    Eigen::Matrix<float, 5, 1> face_vals_prim;
    for (size_t i=0; i<this->neighbor.size();i++){
        if (this->neighbor[i] == nullptr) {
            continue;
        }

        // iterated through a list, assigned value to both
        if (this->neighbor[i]->ind > this->ind){
            // face_val = (this->U + this->neighbor[i]->U)/2;
            // // face_val = this->U;
            // face_vals_prim = (this->prim + this->neighbor[i]->prim)/2;

            // determine upwind 
            if (this->U[i%2+1]+this->neighbor[i]->U[i%2+1]>0){
                // this cell upwind
                face_val = this->U;
                face_vals_prim = this->prim;
            }
            else {
                // neighbor cell downwind
                face_val = this->neighbor[i]->U;
                face_vals_prim = this->neighbor[i]->prim;
            }

            // // determine upwind 
            // if (this->U[i%2+1]>0){
            //     // this cell upwind
            //     face_val = 0.8*this->U+0.2*this->neighbor[i]->U;
            //     face_vals_prim = 0.8*this->prim+0.2*this->neighbor[i]->prim;
            // }
            // else {
            //     // neighbor cell downwind
            //     face_val = 0.2*this->U+0.8*this->neighbor[i]->U;
            //     face_vals_prim = 0.2*this->prim+0.8*this->neighbor[i]->prim;
            // }

            // conservative face vals
            this->face_vals[i] = face_val;
            this->neighbor[i]->face_vals[(i+2)%4] = face_val;

            // face_vals_prim = this->prim;

            this->face_vals_prim[i] = face_vals_prim;
            this->neighbor[i]->face_vals_prim[(i+2)%4] = face_vals_prim;

            if (face_val.array().isNaN().any()){
                cout << "nan" << endl;
            }
            if (face_vals_prim.array().isNaN().any()){
                cout << "nan" << endl;
            }

            // if (i%2==0){
            //     // horizontal

            //     // determine upwind 
            //     if (this->U[i%2+1]>0){
            //         // this cell upwind
            //         face_val = this->U;
            //         face_vals_prim = this->prim;
            //     }
            //     else {
            //         // neighbor cell downwind
            //         face_val = this->neighbor[i]->U;
            //         face_vals_prim = this->neighbor[i]->prim;
            //     }
            // }

            if (this->mu!=0){
                Matrix2f face_gradient = (this->gradient + this->neighbor[i]->gradient)/2;

                this->face_gradient[i] = face_gradient;
                this->neighbor[i]->face_gradient[(i+2)%4] = face_gradient;
            }
        }
    }
}

void CFDMesh::computeFlux(){
    this->computeConvFlux();

    if (this->mu!=0){
        this->computeViscFlux();
    }
}

void CFDMesh::computeConvFlux(){
    Vector4f face_flux;
    Vector4f face_vals;
    Eigen::Matrix<float, 5, 1> face_vals_prim;

    for (size_t i=0; i<this->neighbor.size();i++){
        if (this->neighbor[i] == nullptr) {
            this->flux[i] = Vector4f::Zero(); // No flux at boundaries
            // this->flux[i](seq(1,2)) = this->prim[4];
            // this->flux[i][i%2+1] = (this->prim[4])*dx;
            this->flux[i][i%2+1] = (this->prim[4]+this->U[i%2+1]*this->prim[i%2+1])*this->dx;
            if (i>1){
                this->flux[i][i%2+1] *= -1;
            }
            continue;
        }
        // right, up, left, down
        // iterated through a list, assigned value to both
        if (this->neighbor[i]->ind > this->ind){
            face_flux = Vector4f::Zero();

            face_vals = this->face_vals[i];
            face_vals_prim = this->face_vals_prim[i];
            // i%2+1 gives either rhou or rhov velocity if i==0, gives rhou, if i==1, gives rhov
            face_flux[0] = face_vals[i%2+1]*dx;
            // gives convective fluxes for normal direction (cartesian grid) otherwise 0
            face_flux[i%2+1] = (face_vals[i%2+1]*face_vals_prim[i%2+1]+face_vals_prim[4])*dx;
            // need to include rhouv terms?
            face_flux[(i+1)%2+1] = (face_vals[(i+1)%2+1]*face_vals_prim[i%2+1])*dx;

            float temp_grad = (this->neighbor[i]->temp-this->temp)/(this->neighbor[i]->center[i%2]-this->center[i%2]);
            temp_grad = 0;

            // energy transfer from mass transfer, pressure work, conduction
            face_flux[3] = (face_vals[3]*face_vals_prim[i%2+1] + this->prim(4)*face_vals_prim[i%2+1] - this->k*temp_grad)*dx;

            if (i>1){
                // flip direction due to normal face vector negative
                face_flux *= -1;
            }

            this->flux[i] = face_flux;
            this->neighbor[i]->flux[(i+2)%4] = -face_flux;

            if (face_flux.array().isNaN().any()){
                cout << "nan" << endl;
            }
        }
    }
}

void CFDMesh::computeGradient(){
    // compute gradients at cell centers
    // du/dx
    this->gradient(0,0) = (this->neighbor[0]->prim[1] - this->neighbor[2]->prim[1])/this->dx/2;
    // dv/dy
    this->gradient(1,1) = (this->neighbor[1]->prim[2] - this->neighbor[3]->prim[2])/this->dx/2;
    // du/dy
    this->gradient(1,0) = (this->neighbor[1]->prim[1] - this->neighbor[3]->prim[1])/this->dx/2;
    // dv/dx
    this->gradient(0,1) = (this->neighbor[0]->prim[2] - this->neighbor[2]->prim[2])/this->dx/2;

}


void CFDMesh::computeViscFlux(){
    // compute visc function
    Vector4f face_flux;
    Vector4f face_vals;
    Eigen::Matrix<float, 5, 1> face_vals_prim;

    for (size_t i=0; i<this->neighbor.size();i++){
        // iterated through a list, assigned value to both
        if (this->neighbor[i]->ind > this->ind){
            face_flux = Vector4f::Zero();

            // recompute face gradients of orthogonal velocity gradients
            if (i==0){
                this->face_gradient[i](0,0) = (this->neighbor[i]->prim(1) - this->prim(1))/dx;
            }
            else if (i==2){
                this->face_gradient[i](0,0) = (this->prim(1) - this->neighbor[i]->prim(1))/dx;
            }
            else if (i==1){
                this->face_gradient[i](1,1) = (this->neighbor[i]->prim(2) - this->prim(2))/dx;
            }
            else if (i==3){
                this->face_gradient[i](1,1) = (this->prim(2) - this->neighbor[i]->prim(2))/dx;
            }

            float txy = this->mu*(this->face_gradient[i](0,1)+this->face_gradient[i](1,0));

            if (i==0 || i==2){
                float txx = this->mu*(2*(this->face_gradient[i](0,0))-2/3*(this->face_gradient[i](0,0)+this->face_gradient[i](1,1)));
                face_flux[1] = txx*dx;
                face_flux[2] = txy*dx;
                face_flux[3] = this->prim(1)*txx+this->prim(2)*txy;
            }
            else if (i==1 || i==3){
                float tyy = this->mu*(2*(this->face_gradient[i](1,1))-2/3*(this->face_gradient[i](0,0)+this->face_gradient[i](1,1)));
                face_flux[1] = txy*dx;
                face_flux[2] = tyy*dx;
                face_flux[3] = this->prim(1)*txy+this->prim(2)*tyy;
            }

            if (i>2){
                // flip direction due to normal face vector negative
                face_flux *= -1;
            }

            this->flux[i] -= face_flux;
            this->neighbor[i]->flux[(i+2)%4] -= -face_flux;
        }
    }
}

void CFDMesh::computeSource(){
    // compute source term
    
    // no source term
    this->J = Vector4f::Zero();
}

void CFDMesh::computeDE(){
    // compute dUdt using fluxes
    this -> dUdt = Vector4f::Zero();

    for (Vector4f flux_indiv: this->flux){
        this->dUdt -= flux_indiv; 
    }
    // add source contribution
    this->dUdt += this->J * this->vol;
    this->dUdt/=this->vol;
}
