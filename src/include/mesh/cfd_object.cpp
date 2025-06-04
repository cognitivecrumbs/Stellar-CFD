#include <cfd_object.h>
#include <initialization.h>
#include <iomanip> // Required for std::setprecision

// std::mutex mtx; // Global mutex

// Find maximum of three numbers
float CFDTools::maxs(float a, float b, float c) {
    float k = a;
    if (b > k) k = b;
    if (c > k) k = c;
    return k;
}

// Find minimum of three numbers
float CFDTools::mins(float a, float b, float c) {
    float k = a;
    if (b < k) k = b;
    if (c < k) k = c;
    return k;
}

// Sign(x) function
int CFDTools::sgn(float v) {
    if (v < 0) return -1;
    if (v > 0) return 1;
    return 0;
}

// Minmod Limiter
float CFDTools::minmod(float a, float b, float c) {
    float vs[] = {std::abs(a), std::abs(b), std::abs(c)};
    float v = mins(vs[0], vs[1], vs[2]);
    float s = (sgn(a) + sgn(b) + sgn(c)) / 3;
    return (std::abs(s) == 1.0) ? s * v : 0.0f;
}

// HLLE Solver
template <size_t size_t>
void CFDTools::HLLE(float (&flux)[size_t], float qL[4], float qR[4], float gamma, int normal_x, int normal_y) {
    int nx;
    int ny;
    float rL,uL,vL,pL;
    float rR,uR,vR,pR;
    float vnL;
    float vnR;
    float vn;
    float aL,HL,HR;
    float aR,RT;
    float u,v,H,a;

    // normal vectors
    nx = normal_x;
    ny = normal_y;

    // Left State
    rL = qL[0];
    uL = qL[1]/rL;
    vL = qL[2]/rL;
    vnL = uL*nx + vL*ny;
    pL = (gamma-1)*(qL[3] - rL*(uL*uL + vL*vL)/2);
    // pL = std::max(pL,float(1e-16));
    aL= sqrt((gamma*pL)/rL);
    HL = (qL[3] + pL)/rL;

    // Right State
    rR = qR[0];
    uR = qR[1]/rR;
    vR = qR[2]/rR;
    vnR= uR*nx + vR*ny;
    pR = (gamma-1)*(qR[3] - rR*(uR*uR + vR*vR)/2);
    // pR = std::max(pR,float(1e-16));
    aR = sqrt(gamma*pR/rR);
    HR = (qR[3] + pR)/rR;


    // if (gamma*pR/rR<=0){
    //     cout << "stop" << endl;
    // }
    // if ((gamma*pL)/rL<=0){
    //     cout << "stop" << endl;
    // }

    // First compute the Roe Averages
    RT = sqrt(rR/rL);
    u = (uL+RT*uR)/(1+RT);
    v = (vL+RT*vR)/(1+RT);
    H = (HL+RT*HR)/(1+RT);
    a = sqrt((gamma-1)*(H - (u*u + v*v)/2));
    vn = u*nx + v*ny;

    if (std::isnan(a) || std::isnan(RT) || std::isnan(aL) || std::isnan(aR)) {
        for(int i=0; i<=3; i++) {
            flux[i] = 0;
        }
        return;
    }
    // if (rR/rL<=0){
    //     cout << "stop" << endl;
    // }
    // if ((gamma-1)*(H - (u*u + v*v)/2)<=0){
    //     cout << "stop" << endl;
    // }


    float SLm;
    float SRp;

    // Wave speed estimates
    SLm = mins((vnL-aL),(vn - a),0.0);
    SRp = maxs((vnR+aR), (vn + a), 0.0);

    // Left and Right fluxes
    float FL[4] = {rL*vnL, rL*vnL*uL + pL*nx, rL*vnL*vL + pL*ny, rL*vnL*HL};
    float FR[4] = {rR*vnR, rR*vnR*uR + pR*nx, rR*vnR*vR + pR*ny, rR*vnR*HR};

    if ((SRp-SLm) == 0) {
        printf("DIVIDINGS BY 0: SLm = %lf, SRp = %lf  \n", SLm, SRp);
        for(int j=0; j<=3; j++) {
            // flux[j] = 0;
            flux[j] = (FL[j]+FR[j])/2;
        }
        return;
    }

    // Compute HLL Flux
    for(int i=0; i<=3; i++) {
        flux[i] = ((SRp)*FL[i] - (SLm) * FR[i] + (SLm)* (SRp)*(qR[i] - qL[i])) / ((SRp) - (SLm));
        if (std::isinf(flux[i]) || std::isnan(flux[i])){
            // cout << "HJL" << endl;
            for(int j=0; j<=3; j++) {
                flux[j] = 0;
            }
            return;
        }
    }

    // for(int i=0; i<=3; i++) {
    //     flux[i] = 0.5*(FL[i]+FR[i]) - 0.5*max(abs(SLm),abs(SRp))*(qR[i]-qL[i]);
    // }

}

template <size_t size_t>
void CFDTools::VISC(float (&flux)[size_t], float tii_face[3], float dtdi_face[2], float qL[4], float qR[4], int nx, int ny)
{
    float rL,uL,vL;
    float rR,uR,vR;
    float vnL;
    float vnR;
    // float vn;
    float RT;
    float u,v;

    // Left State
    rL = qL[0];
    uL = qL[1]/rL;
    vL = qL[2]/rL;

    // Right State
    rR = qR[0];
    uR = qR[1]/rR;
    vR = qR[2]/rR;

    // First compute the Roe Averages
    RT = sqrt(rR/rL);
    u = (uL+RT*uR)/(1+RT);
    v = (vL+RT*vR)/(1+RT);

    // u = (uL+uR)/2;
    // v = (vL+vR)/2;

    // vn = u*nx + v*ny;

    // if (tii_face[0]*nx + tii_face[2]*ny != 0){
    //     cout << "here" << endl;
    // }
    // if (tii_face[1]*ny + tii_face[2]*nx != 0){
    //     cout << "here" << endl;
    // }
    // if ((tii_face[0]*u + tii_face[2]*v)*nx + (tii_face[2]*u + tii_face[1]*v)*ny + dtdi_face[0]*nx + dtdi_face[0]*ny != 0){
    //     cout << "here" << endl;
    // }

    flux[1] -= tii_face[0]*nx + tii_face[2]*ny;
    flux[2] -= tii_face[1]*ny + tii_face[2]*nx;
    flux[3] -= (tii_face[0]*u + tii_face[2]*v)*nx + (tii_face[2]*u + tii_face[1]*v)*ny + dtdi_face[0]*nx + dtdi_face[0]*ny;
}

// Flux solver
// float*** CFDTools::FLUXSOLVER(float*** U, int nx, int ny, int ghost_cell, float dx, float dy, float gamma) {
// void CFDTools::FLUXSOLVER(float*** dUdt, float*** U, float***dUdx, float***dUdy, int nx, int ny, int ghost_cell, float dx, float dy, float gamma) {
// void CFDTools::FLUXSOLVER(float*** U, float*** dUdt, float***dUdx, float***dUdy, int nx, int ny, int ghost_cell, float dx, float dy, float gamma) {
void CFDTools::FLUXSOLVER(float*** U, 
                          float** e, 
                          float*** dUdt, 
                          float***dUdx, 
                          float***dUdy, 
                          float***tii, 
                          float***dtdi, 
                          int nx, 
                          int ny, 
                          int ghost_cell, 
                          float dx, 
                          float dy, 
                          float gamma) {
    float dqw;
    float dqe;
    float dqc;
    float dqn;
    float dqs;
    // float*** dUdt;
    // float dqdx[ny+2*ghost_cell][nx+2*ghost_cell][4];
    // float dqdy[ny+2*ghost_cell][nx+2*ghost_cell][4];
    // float*** dqdx;
    // float*** dqdy;
    float qxL[4];
    float qxR[4];
    float qyL[4];
    float qyR[4];
    float qL[4];
    float qR[4];
    float flux[4];

    // // dUdt = new float**[ny+2*ghost_cell];
    // dqdx = new float**[ny+2*ghost_cell];
    // dqdy = new float**[ny+2*ghost_cell];
    // for (int i=0; i<ny+2*ghost_cell; i++) {
    //     // dUdt[i] = new float*[nx+2*ghost_cell];
    //     dqdx[i] = new float*[nx+2*ghost_cell];
    //     dqdy[i] = new float*[nx+2*ghost_cell];
    // }
    // for (int i=0; i<ny+2*ghost_cell; i++) {
    //     for (int j = 0; j < nx+2*ghost_cell; j++) {
    //         // dUdt[i][j] = new float[4];
    //         dqdx[i][j] = new float[4];
    //         dqdy[i][j] = new float[4];
    //     }
    // }
    float gen_offset = 4;

    for (int i = 0; i < ny+2*ghost_cell; i++) {
        for (int j = 0; j < nx+2*ghost_cell; j++) {
            for (int k = 0; k < 4; k++) {
                dUdx[i][j][k] = 0;
                dUdy[i][j][k] = 0;
                dUdt[i][j][k] = 0;
            }

            e[i][j] = (U[i][j][3] - (U[i][j][1]*U[i][j][1] + U[i][j][2]*U[i][j][2])/(2*U[i][j][0]))/U[i][j][0];

            // dUdt[i][j][3] -= (e[i][j]*e[i][j]*e[i][j]*e[i][j]/(1e4)-(max(U[i][j][0],0.4f)-0.4f)*e[i][j]*e[i][j]*sqrt(e[i][j]))/200000000.f;
        }
    }

    for (int k = 0; k < 4; k++) {
        qxL[k] = 0;
        qxR[k] = 0;
        qyL[k] = 0;
        qyR[k] = 0;
        qL[k]  = 0;
        qR[k]  = 0;
    }

    // to find viscous effects compute stresses at node centers (viscous shears, normals, and heat transfer)
    // interpolate to get face values for stresses(central method)
    // may not be most accurate but viscous effects are diffusive so should not introduce oscillation
    // similar to overdamped systems

    float k_cv = 1e-2;
    float mu = 1e-10;
    k_cv = 1e-4;
    mu = 1e-4;
    // k_cv = 0;

    float dudx, dudy, dvdx, dvdy;
    float dudxL, dudyL, dvdxL, dvdyL;
    float dudxR, dudyR, dvdxR, dvdyR;
    float tii_face[3];
    float dtdi_face[2];


    // // float e;
    // for (int i = ghost_cell; i < ny+2*ghost_cell-1; i++) {
    //     for (int j = ghost_cell; j < nx+2*ghost_cell-1; j++) {
    //         e[i][j] = (U[i][j][3]-(U[i][j][1]*U[i][j][1]+U[i][j][2]*U[i][j][2])/(2*U[i][j][0]))/U[i][j][0];
    //     }
    // }

    for (int i = 1; i < ny+2*ghost_cell-1; i++) {
        for (int j = 1; j < nx+2*ghost_cell-1; j++ ) {
            dtdi[i][j][0] = (e[i][j+1] - e[i][j-1])/(2*dy);
            dtdi[i][j][1] = (e[i+1][j] - e[i-1][j])/(2*dx);

            dudx = (U[i][j+1][1]/U[i][j+1][0] - U[i][j-1][1]/U[i][j-1][0])/(2*dx);
            dudy = (U[i+1][j][1]/U[i+1][j][0] - U[i-1][j][1]/U[i-1][j][0])/(2*dy);
            dvdx = (U[i][j+1][2]/U[i][j+1][0] - U[i][j-1][2]/U[i][j-1][0])/(2*dx);
            dvdy = (U[i+1][j][2]/U[i+1][j][0] - U[i-1][j][2]/U[i-1][j][0])/(2*dy);

            dudxL = (U[i][j][1]/U[i][j][0] - U[i][j-1][1]/U[i][j-1][0])/dx;
            dudyL = (U[i][j][1]/U[i][j][0] - U[i-1][j][1]/U[i-1][j][0])/dy;
            dvdxL = (U[i][j][2]/U[i][j][0] - U[i][j-1][2]/U[i][j-1][0])/dx;
            dvdyL = (U[i][j][2]/U[i][j][0] - U[i-1][j][2]/U[i-1][j][0])/dy;

            dudxR = (U[i][j+1][1]/U[i][j+1][0] - U[i][j][1]/U[i][j][0])/dx;
            dudyR = (U[i+1][j][1]/U[i+1][j][0] - U[i][j][1]/U[i][j][0])/dy;
            dvdxR = (U[i][j+1][2]/U[i][j+1][0] - U[i][j][2]/U[i][j][0])/dx;
            dvdyR = (U[i+1][j][2]/U[i+1][j][0] - U[i][j][2]/U[i][j][0])/dy;

            dudx = minmod(dudx,dudxL,dudxR);
            dvdx = minmod(dvdx,dvdxL,dvdxR);
            dudy = minmod(dudy,dudyL,dudyR);
            dvdy = minmod(dvdy,dvdyL,dvdyR);

            tii[i][j][0] = -2/3*mu*(dudx+dvdy) + 2*mu*dudx;
            tii[i][j][1] = -2/3*mu*(dudx+dvdy) + 2*mu*dvdy;
            tii[i][j][2] = mu*(dvdx+dudy);

            // if (tii[i][j][0]!=0 ||tii[i][j][1]!=0 ||tii[i][j][2]!=0 ){
            //     cout << "here" << endl;
            // }
        }
    }

    // Compute and limit slopes at cells (i,j)
    // for (int i = ghost_cell; i < ny+2*ghost_cell-1; i++) {
    //     for (int j = ghost_cell; j < nx+2*ghost_cell-1; j++) {
    for (int i = ghost_cell; i < ny+ghost_cell; i++) {
        for (int j = ghost_cell; j < nx+ghost_cell; j++) {
            for (int k = 0; k <= 3; k++) {

                dqw = (U[i][j][k] - U[i][j - 1][k])/dx;
                dqe = (U[i][j + 1][k] - U[i][j][k])/dx;
                // dqc = (U[i][j + 1][k] - U[i][j - 1][k]) /(2 * dx);
                dqc = (U[i][j + 1][k] - U[i][j - 1][k]) /(2 * dx);
                // dqdx[i][j][k] = minmod(dqw, dqe, dqc)/2;
                dUdx[i][j][k] = minmod(dqw, dqe, dqc);

                dqs = (U[i][j][k] - U[i - 1][j][k])/ dy;
                dqn = (U[i + 1][j][k] - U[i][j][k])/ dy;
                // dqc = (U[i + 1][j][k] - U[i - 1][j][k])/(2 * dy);
                dqc = (U[i + 1][j][k] - U[i - 1][j][k])/(2 * dy);
                // dqdy[i][j][k] = minmod(dqs, dqn, dqc)/2;
                dUdy[i][j][k] = minmod(dqs, dqn, dqc);



                // // if (i==60 && j != 1){
                // //     cout << "here" << endl;
                // // }
                // dqw = 2*(U[i][j][k] - U[i][j - 1][k])/dx;
                // dqe = 2*(U[i][j + 1][k] - U[i][j][k])/dx;
                // // dqc = (U[i][j + 1][k] - U[i][j - 1][k]) /(dx);
                // // dqc = 2*(U[i][j + 1][k] - U[i][j - 1][k]) /(4 * dx);
                // dqc = 2*(U[i][j + 1][k] - U[i][j - 1][k]) /(2 * dx);
                // // dqdx[i][j][k] = minmod(dqw, dqe, dqc)/2;
                // dqdx[i][j][k] = minmod(dqw, dqe, dqc);

                // dqs = 2*(U[i][j][k] - U[i - 1][j][k])/ dy;
                // dqn = 2*(U[i + 1][j][k] - U[i][j][k])/ dy;
                // // dqc = (U[i + 1][j][k] - U[i - 1][j][k])/(dy);
                // // dqc = 2*(U[i + 1][j][k] - U[i - 1][j][k])/(4 * dy);
                // dqc = 2*(U[i + 1][j][k] - U[i - 1][j][k])/(2 * dy);
                // // dqdy[i][j][k] = minmod(dqs, dqn, dqc)/2;
                // dqdy[i][j][k] = minmod(dqs, dqn, dqc);
            }
        }
    }


    // ThreadPool pool = ThreadPool(8);
    // Compute fluxes x-direction
    // face values
    for (int i = ghost_cell; i < ny + ghost_cell; i++) {
        for (int j = ghost_cell-1; j < nx + ghost_cell; j++) {
            
            for (int k = 0; k <= 3; k++) {
                qxL[k] = U[i][j][k] + (dUdx[i][j][k])*dx/2;
                qxR[k] = U[i][j + 1][k] - (dUdx[i][j + 1][k])*dx/2;
            }

            HLLE(flux, qxL, qxR, gamma, 1, 0);

            // compute stress tensor and energy gradient at face
            for (int k=0; k<3; k++){
                tii_face[k] = (tii[i][j][k]+tii[i][j+1][k])/2;
            }
            for (int k=0; k<2; k++){
                dtdi_face[k] = (dtdi[i][j][k]+dtdi[i][j+1][k])/2*k_cv;
            }
               
            VISC(flux,tii_face,dtdi_face,qxL,qxR,1,0);

            for (int k = 0; k <= 3; k++) {
                dUdt[i][j][k] = dUdt[i][j][k] + (flux[k])/ dx;
                dUdt[i][j + 1][k] = dUdt[i][j + 1][k] - (flux[k])/dx;
                if (std::isnan(dUdt[i][j][k])){
                    cout << "HJL" << endl;
                }
            }

            // pool.enqueue([i,j,U,dUdt,dUdx,dUdy,dx,dy,gamma,this,tii,dtdi,k_cv](){
            //     float tii_face[3];
            //     float dtdi_face[2];
            //     float qxL[4];
            //     float qxR[4];
            //     float flux[4];

            //     for (int k = 0; k <= 3; k++) {
            //         qxL[k] = U[i][j][k] + (dUdx[i][j][k])*dx/2;
            //         qxR[k] = U[i][j + 1][k] - (dUdx[i][j + 1][k])*dx/2;
            //     }

            //     this->HLLE(flux, qxL, qxR, gamma, 1, 0);

            //     // compute stress tensor and energy gradient at face
            //     for (int k=0; k<3; k++){
            //         tii_face[k] = (tii[i][j][k]+tii[i][j+1][k])/2;
            //     }
            //     for (int k=0; k<2; k++){
            //         dtdi_face[k] = (dtdi[i][j][k]+dtdi[i][j+1][k])/2*k_cv;
            //     }
                
            //     this->VISC(flux,tii_face,dtdi_face,qxL,qxR,1,0);

            //     // std::lock_guard<std::mutex> lock(mtx);
            //     for (int k = 0; k <= 3; k++) {
            //         dUdt[i][j][k] = dUdt[i][j][k] + (flux[k])/ dx;
            //         dUdt[i][j + 1][k] = dUdt[i][j + 1][k] - (flux[k])/dx;
            //         if (std::isnan(dUdt[i][j][k])){
            //             cout << "HJL" << endl;
            //         }
            //     }
            //     // std::lock_guard<std::mutex> unlock(mtx);
            // });
        }
    }
    // pool.stopPool();

    // pool.restart(8);
    // Compute fluxes y-direction
    for (int i = ghost_cell-1; i < ny + ghost_cell; i++) {
        for (int j = ghost_cell; j < nx + ghost_cell; j++) {
            for (int k = 0; k <= 3; k++) {
                qyL[k] = U[i][j][k] + dUdy[i][j][k]*dy/2;
                qyR[k] = U[i + 1][j][k] - dUdy[i + 1][j][k]*dy/2;
            }
            HLLE(flux, qyL, qyR, gamma, 0, 1);
            // if (U[i][j][0] < 1e-6 && flux[0]>0 || U[i+1][j][0] < 1e-6 && flux[0]<0){
            //     continue;
            // }

            // compute stress tensor and energy gradient at face
            for (int k=0; k<3; k++){
                tii_face[k] = (tii[i][j][k]+tii[i+1][j][k])/2;
            }
            for (int k=0; k<2; k++){
                dtdi_face[k] = (dtdi[i][j][k]+dtdi[i+1][j][k])/2*k_cv;
            }
               
            VISC(flux,tii_face,dtdi_face,qyL,qyR,0,1);

            for (int k = 0; k <= 3; k++) {
                dUdt[i][j][k] = dUdt[i][j][k] + (flux[k])/dy;
                dUdt[i + 1][j][k] = dUdt[i + 1][j][k] - (flux[k])/dy;
                if (std::isnan(dUdt[i][j][k])){
                    cout << "HJL" << endl;
                }
            }

            // pool.enqueue([i,j,U,dUdt,dUdx,dUdy,dx,dy,gamma,this,tii,dtdi,k_cv](){
            //     float tii_face[3];
            //     float dtdi_face[2];
            //     float qyL[4];
            //     float qyR[4];
            //     float flux[4];

            //     for (int k = 0; k <= 3; k++) {
            //         qyL[k] = U[i][j][k] + dUdy[i][j][k]*dy/2;
            //         qyR[k] = U[i + 1][j][k] - dUdy[i + 1][j][k]*dy/2;
            //     }
            //     this->HLLE(flux, qyL, qyR, gamma, 0, 1);
            //     // if (U[i][j][0] < 1e-6 && flux[0]>0 || U[i+1][j][0] < 1e-6 && flux[0]<0){
            //     //     continue;
            //     // }

            //     // compute stress tensor and energy gradient at face
            //     for (int k=0; k<3; k++){
            //         tii_face[k] = (tii[i][j][k]+tii[i+1][j][k])/2;
            //     }
            //     for (int k=0; k<2; k++){
            //         dtdi_face[k] = (dtdi[i][j][k]+dtdi[i+1][j][k])/2*k_cv;
            //     }
                
            //     this->VISC(flux,tii_face,dtdi_face,qyL,qyR,0,1);

            //     // std::lock_guard<std::mutex> lock(mtx);
            //     for (int k = 0; k <= 3; k++) {
            //         dUdt[i][j][k] = dUdt[i][j][k] + (flux[k])/dy;
            //         dUdt[i + 1][j][k] = dUdt[i + 1][j][k] - (flux[k])/dy;
            //         if (std::isnan(dUdt[i][j][k])){
            //             cout << "HJL" << endl;
            //         }
            //     }
            //     // std::lock_guard<std::mutex> unlock(mtx);
            // });
        }
    }
    // pool.stopPool();
}

void CFDSimulation::solvePoissonEquation(float***& U, float**& phi, int nx, int ny, int ghost_cell, float dx, float dy, float tol, int max_iters) {
    const float G = 2.0f; // Gravitational constant
    float** rho = new float*[ny + 2 * ghost_cell];

    // Allocate memory for rho
    for (int i = 0; i < ny + 2 * ghost_cell; i++) {
        rho[i] = new float[nx + 2 * ghost_cell]();
    }

    // Compute mass density from U
    for (int i = ghost_cell; i < ny + ghost_cell; i++) {
        for (int j = ghost_cell; j < nx + ghost_cell; j++) {
            rho[i][j] = U[i][j][0]; // Extract mass density
        }
    }

    // Precompute scaling factors
    const float dx2 = dx * dx;
    const float dy2 = dy * dy;
    const float denom = 2.0f * (1.0f / dx2 + 1.0f / dy2);

    // Iteratively solve the Poisson equation using Gauss-Seidel iteration
    for (int iter = 0; iter < max_iters; iter++) {
        float max_diff = 0.0f;

        // Update interior points
        for (int i = ghost_cell; i < ny + ghost_cell; i++) {
            for (int j = ghost_cell; j < nx + ghost_cell; j++) {
                float old_phi = phi[i][j];

                phi[i][j] = ((phi[i + 1][j] + phi[i - 1][j]) / dx2 +
                             (phi[i][j + 1] + phi[i][j - 1]) / dy2 -
                             4.0f * M_PI * G * rho[i][j]) /
                            denom;

                max_diff = std::max(max_diff, std::fabs(phi[i][j] - old_phi));
            }
        }

        // Enforce periodic boundary conditions after updating interior points
        for (int i = ghost_cell; i < ny + ghost_cell; i++) {
            for (int k = 0; k < ghost_cell; k++) {
                phi[i][ghost_cell - k - 1] = phi[i][nx + ghost_cell - k - 1]; // Left boundary
                phi[i][nx + ghost_cell + k] = phi[i][ghost_cell + k];        // Right boundary
            }
        }
        for (int j = ghost_cell; j < nx + ghost_cell; j++) {
            for (int k = 0; k < ghost_cell; k++) {
                phi[ghost_cell - k - 1][j] = phi[ny + ghost_cell - k - 1][j]; // Bottom boundary
                phi[ny + ghost_cell + k][j] = phi[ghost_cell + k][j];         // Top boundary
            }
        }

        // Check for convergence
        if (max_diff < tol) {
            // std::cout << "Poisson solver converged in " << iter + 1 << " iterations.\n";
            break;
        }

        // cout << max_diff << endl;
    }

    // Clean up allocated memory
    for (int i = 0; i < ny + 2 * ghost_cell; i++) {
        delete[] rho[i];
    }
    delete[] rho;
}


#include <iostream>
#include <cmath>
#include <vector>
#include <kiss_fftnd.h>

void CFDSimulation::solvePoissonEquationFFT(float***& U, float**& phi, int nx, int ny, int ghost_cell, float dx, float dy) {
    // const float G = 4.0e0; // Gravitational constant
    const float G = 15; // Gravitational constant
    // const float G = 2; // Gravitational constant
    // const float G = 5; // Gravitational constant

    // ghost_cell = 1;
    // nx = 242;
    // ny = 242;

    // Compute rho from U
    std::vector<std::vector<float>> rho(ny + 2 * ghost_cell, std::vector<float>(nx + 2 * ghost_cell, 0.0f));
    for (int i = ghost_cell; i < ny + ghost_cell; i++) {
        for (int j = ghost_cell; j < nx + ghost_cell; j++) {
            rho[i][j] = U[i][j][0]; // Extract mass density
            // if (rho[i][j]<0){
            //     cout << "stop" << endl;
            // }
        }
    }

    // Define effective domain sizes (no ghost cells for FFT)
    const int effective_nx = nx;
    const int effective_ny = ny;

    // Prepare input and output arrays for FFT
    std::vector<kiss_fft_cpx> input(effective_nx * effective_ny), output(effective_nx * effective_ny);
    std::vector<kiss_fft_cpx> inverse_input(effective_nx * effective_ny), inverse_output(effective_nx * effective_ny);

    // Copy rho into the input array for FFT
    for (int i = 0; i < effective_ny; i++) {
        for (int j = 0; j < effective_nx; j++) {
            input[i * effective_nx + j].r = rho[i + ghost_cell][j + ghost_cell];
            input[i * effective_nx + j].i = 0.0f;
        }
    }

    // Create FFT plans
    int dims[2] = {effective_ny, effective_nx};
    kiss_fftnd_cfg forward_plan = kiss_fftnd_alloc(dims, 2, 0, nullptr, nullptr);
    kiss_fftnd_cfg inverse_plan = kiss_fftnd_alloc(dims, 2, 1, nullptr, nullptr);

    // Perform forward FFT
    kiss_fftnd(forward_plan, input.data(), output.data());

    // Solve in Fourier space
    for (int i = 0; i < effective_ny; i++) {
        for (int j = 0; j < effective_nx; j++) {
            int kx = (j <= effective_nx / 2) ? j : j - effective_nx;
            int ky = (i <= effective_ny / 2) ? i : i - effective_ny;

            float kx_scaled = (2.0f * M_PI * kx) / (effective_nx * dx);
            float ky_scaled = (2.0f * M_PI * ky) / (effective_ny * dy);

            float k_squared = kx_scaled * kx_scaled + ky_scaled * ky_scaled;

            if (k_squared != 0.0f) {
                float scale = -4.0f * M_PI * G / k_squared;
                output[i * effective_nx + j].r *= scale;
                output[i * effective_nx + j].i *= scale;
            } else {
                // Handle the zero-frequency mode (average value)
                output[i * effective_nx + j].r = 0.0f;
                output[i * effective_nx + j].i = 0.0f;
            }
        }
    }

    // Perform inverse FFT
    kiss_fftnd(inverse_plan, output.data(), inverse_output.data());

    // Copy the result back to phi
    for (int i = 0; i < effective_ny; i++) {
        for (int j = 0; j < effective_nx; j++) {
            phi[i + ghost_cell][j + ghost_cell] = inverse_output[i * effective_nx + j].r / (effective_nx * effective_ny);
        }
    }

    // Cleanup
    kiss_fft_free(forward_plan);
    kiss_fft_free(inverse_plan);

    // for (int i = 0; i < effective_ny; i++) {
    //     phi[i + ghost_cell][ghost_cell-1] = phi[i + ghost_cell][effective_nx+ghost_cell-1];
    //     phi[i + ghost_cell][effective_nx+ghost_cell] = phi[i + ghost_cell][ghost_cell];
    // }
    // for (int j = 0; j < effective_nx; j++) {
    //     phi[ghost_cell-1][j+ghost_cell] = phi[effective_ny+ghost_cell-1][j+ghost_cell];
    //     phi[effective_ny+ghost_cell][j+ghost_cell] = phi[ghost_cell][j+ghost_cell];
    // }  
}





void CFDSimulation::computeGravitationalAcceleration(float***& U, float**& phi, int nx, int ny, int ghost_cell, float dx, float dy) {
    for (int i = ghost_cell; i < ny + ghost_cell; i++) {
        for (int j = ghost_cell; j < nx + ghost_cell; j++) {
            // if (j == 47 && i == 52 || j == 57 && i == 52){
            // if (j == 52 && i == 52){
            //     cout << "stop" << endl;
            // }

            if (U[i][j][0]<=1e-6){
                continue;;
            }


            // float gx, gy;
            // if (i == ghost_cell){
            //     gy = -(phi[i+1][j] - phi[i][j]) / (dy);
            // }
            // else if (i == ny+ghost_cell-1) {
            //     gy = -(phi[i][j] - phi[i-1][j]) / (dy);
            // }
            // else {
            //     gy = -(phi[i+1][j] - phi[i-1][j]) / (2 * dy);
            // }

            // if (j == ghost_cell){
            //     gx = -(phi[i][j+1] - phi[i][j]) / (dx);
            // }
            // else if (j == nx+ghost_cell-1) {
            //     gx = -(phi[i][j] - phi[i][j-1]) / (dx);
            // }
            // else {
            //     gx = -(phi[i][j+1] - phi[i][j-1]) / (2 * dx);
            // }


            // Compute gravitational acceleration
            float gy = -(phi[i+1][j] - phi[i-1][j]) / (2 * dy);
            float gx = -(phi[i][j+1] - phi[i][j-1]) / (2 * dx);

            // gy = -80;
            // gx = 0;

            // float tmp1 = U[i][j][1];
            // float tmp2 = U[i][j][2];
            // float tmp3 = U[i][j][3];



            // Update momentum and energy
            // U[i][j][1] += gx * dt/2; // x-momentum
            // U[i][j][2] += gy * dt/2; // y-momentum
            // U[i][j][3] += (gx * (U[i][j][1] / U[i][j][0]) + gy * (U[i][j][2] / U[i][j][0])) * dt/2; // energy

            U[i][j][1] += gx*U[i][j][0] * dt; // x-momentum
            U[i][j][2] += gy*U[i][j][0] * dt; // y-momentum
            U[i][j][3] += (gx * (U[i][j][1] / U[i][j][0]) + gy * (U[i][j][2] / U[i][j][0]))*U[i][j][0] * dt; // energy

            // for (int c = 0; c < 4; c++) {
            //     if (std::isnan(U[i][j][c]) || std::isinf(U[i][j][c])){
            //         cout << "stop" << endl;
            //     }
            // }
        }
    }
}

CFDSimulation::CFDSimulation(CFDInitializerFull* initializer, int nx, int ny, int ghost_cell, float dx, float dy, float gamma):
initializer(initializer), nx(nx), ny(ny), ghost_cell(ghost_cell), dx(dx), dy(dy), gamma(gamma), CFL(0.8)
{
    this-> dt = 1e-12;
    // this-> dt = 0.0041;

    this-> xc = new float[nx+2*ghost_cell];                      // x array
    this-> yc = new float[ny+2*ghost_cell];                      // x array


    // Initialization
    xc[0]= dx/2;
    yc[0]= dy/2;
    // Fill x,y array
    for (int i=0; i<nx+2*ghost_cell; i++){
        xc[i] = dx/2 + dx*(i-ghost_cell);
    }
    for (int i=0; i<ny+2*ghost_cell; i++){
        yc[i] = dx/2 + dy*(i-ghost_cell);
    }

    this-> e = new float*[ny+2*ghost_cell];

    this->U = new float**[ny+2*ghost_cell];
    this->U_buffer = new float**[ny+2*ghost_cell];
    this->tii = new float**[ny+2*ghost_cell];
    this->dtdi = new float**[ny+2*ghost_cell];
    for (int i=0; i<ny+2*ghost_cell; i++) {
        this->e[i] = new float[nx+2*ghost_cell];
        
        this->U[i] = new float*[nx+2*ghost_cell];
        this->U_buffer[i] = new float*[nx+2*ghost_cell];
        this->tii[i] = new float*[nx+2*ghost_cell];
        this->dtdi[i] = new float*[nx+2*ghost_cell];
    }
    for (int i=0; i<ny+2*ghost_cell; i++) {
        for (int j = 0; j < nx+2*ghost_cell; j++) {
            this->U[i][j] = new float[4];
            this->U_buffer[i][j] = new float[4];
            this->tii[i][j] = new float[3];
            this->dtdi[i][j] = new float[2];
        }
    }
    for (int i=0; i<ny+2*ghost_cell; i++) {
        for (int j = 0; j < nx+2*ghost_cell; j++) {
            for (int k = 0; k < 4; k++){
                // set default to 1. to avoid divide by zero nan in corners of U_buffer
                this->U[i][j][k] = 1.;
                this->U_buffer[i][j][k] = 1.;
            }
        }
    }


    initializer->init(this->U,xc,yc,nx,ny,ghost_cell,gamma);
    initializer->updateBC(U,nx,ny,ghost_cell);

    // float*** dUdt;
    dUdt = new float**[ny+2*ghost_cell];
    dUdx = new float**[ny+2*ghost_cell];
    dUdy = new float**[ny+2*ghost_cell];
    for (int i=0; i<ny+2*ghost_cell; i++) {
        dUdt[i] = new float*[nx+2*ghost_cell];
        dUdx[i] = new float*[nx+2*ghost_cell];
        dUdy[i] = new float*[nx+2*ghost_cell];
    }
    for (int i=0; i<ny+2*ghost_cell; i++) {
        for (int j = 0; j < nx+2*ghost_cell; j++) {
            dUdt[i][j] = new float[4];
            dUdx[i][j] = new float[4];
            dUdy[i][j] = new float[4];
        }
    }
    // for (int i=0; i<ny+2*ghost_cell; i++) {
    //     for (int j = 0; j < nx+2*ghost_cell; j++) {
    //         for (int c = 0; c < 4; c++)
    //         {
    //             dUdt[i][j][c] = 0;
    //         }
    //     }
    // }

    this->phi = new float*[ny + 2 * ghost_cell];
    for (int i = 0; i < ny + 2 * ghost_cell; i++) {
        this->phi[i] = new float[nx + 2 * ghost_cell]();
    }
}

// void CFDSimulation::update(float***& U, int nx, int ny, float dx, float dy, float gamma, float dt) {
void CFDSimulation::update() {
    // float dt = 0.0011;
    // euler method for now
    this->dt_buffer = 1e8;
    float vel,p;

    // dt *= 2;
    // dt = 1e-3;
    // tools.FLUXSOLVER(U, dUdt, dUdx, dUdy, nx, ny, ghost_cell, dx, dy, gamma);
    initializer->updateBC(U,nx,ny,ghost_cell);
    tools.FLUXSOLVER(this->U, 
                     this->e,
                     this->dUdt, 
                     this->dUdx, 
                     this->dUdy, 
                     this->tii, 
                     this->dtdi, 
                     this->nx, 
                     this->ny, 
                     this->ghost_cell, 
                     this->dx, 
                     this->dy, 
                     this->gamma);

    // // Define gravity constant

    for (int y = 0; y < ny+2*ghost_cell; y++) {
        for (int x = 0; x < nx+2*ghost_cell; x++) {
            this->phi[y][x] = 0;
        }
    }
    // this->solvePoissonEquation(U,this->phi,nx,ny,ghost_cell,dx,dy,5e-6,100000);

    this->solvePoissonEquationFFT(U,this->phi,nx,ny,ghost_cell,dx,dy);


    this->computeGravitationalAcceleration(U,this->phi,nx,ny,ghost_cell,dx,dy);

    // for (int y = 0; y < ny+2*ghost_cell; y++) {
    //     for (int x = 0; x < nx+2*ghost_cell; x++) {
    //         for (int c = 0; c < 4; c++) {
    //             // if (std::isnan(U[y][x][c]) || std::isinf(U[y][x][c])){
    //             //     cout << "stop" << endl;
    //             // }

    //             // if (std::isnan(dUdt[y][x][c]) || std::isinf(dUdt[y][x][c])){
    //             //     cout << "stop" << endl;
    //             // }
    //         }
    //     }
    // }

    // for (int y = 0; y < ny+ghost_cell; y++) {
    //     for (int x = 0; x < nx+ghost_cell; x++) {
    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            // if (U[y+ghost_cell][x+ghost_cell][0]<1e-16 && dUdt[y+ghost_cell][x+ghost_cell][0]<0){
            //     continue;
            // }
            for (int c = 0; c < 4; c++) {
                // U[y+ghost_cell][x+ghost_cell][c] = ((U[y+ghost_cell][x+ghost_cell][c] + U[y+ghost_cell][x+ghost_cell][c] - dt * dUdt[y+ghost_cell][x+ghost_cell][c]))/2;
                
                // float tmp = U[y+ghost_cell][x+ghost_cell][c];

                // U[y+ghost_cell][x+ghost_cell][c] = (U[y+ghost_cell][x+ghost_cell][c] - dt * dUdt[y+ghost_cell][x+ghost_cell][c]/2);
                // U[y+ghost_cell][x+ghost_cell][c] -= dt * dUdt[y+ghost_cell][x+ghost_cell][c] / 2;

                U_buffer[y+ghost_cell][x+ghost_cell][c] = (U[y+ghost_cell][x+ghost_cell][c] - dt * dUdt[y+ghost_cell][x+ghost_cell][c]);

                if (std::isnan(U[y+ghost_cell][x+ghost_cell][c]) || std::isinf(U[y+ghost_cell][x+ghost_cell][c])){
                    cout << "stop" << endl;

                }
            }

        }
    }
    
    initializer->updateBC(U_buffer,nx,ny,ghost_cell);
    // float*** dUdt1 = tools.FLUXSOLVER(U_buffer, nx, ny, ghost_cell, dx, dy, gamma);
    // tools.FLUXSOLVER(dUdt, U_buffer, nx, ny, ghost_cell, dx, dy, gamma);
    // tools.FLUXSOLVER(U_buffer, dUdt, dUdx, dUdy, nx, ny, ghost_cell, dx, dy, gamma);
    tools.FLUXSOLVER(this->U_buffer, 
                     this->e,
                     this->dUdt, 
                     this->dUdx, 
                     this->dUdy, 
                     this->tii,
                     this->dtdi,
                     this->nx, 
                     this->ny, 
                     this->ghost_cell, 
                     this->dx, 
                     this->dy, 
                     this->gamma);

    for (int y = 0; y < ny+ghost_cell; y++) {
        for (int x = 0; x < nx+ghost_cell; x++) {
            for (int c = 0; c < 4; c++) {
                // U[y+ghost_cell][x+ghost_cell][c] = (U_buffer[y+ghost_cell][x+ghost_cell][c] + U[y+ghost_cell][x+ghost_cell][c] - dt * dUdt1[y+ghost_cell][x+ghost_cell][c])/2;
                U[y+ghost_cell][x+ghost_cell][c] = (U_buffer[y+ghost_cell][x+ghost_cell][c] + U[y+ghost_cell][x+ghost_cell][c] - dt * dUdt[y+ghost_cell][x+ghost_cell][c])/2;
                // U[y+ghost_cell][x+ghost_cell][c] = (U_buffer[y+ghost_cell][x+ghost_cell][c]);
            }
            if (U[y+ghost_cell][x+ghost_cell][0]<= 1e-16){
                U[y+ghost_cell][x+ghost_cell][0] = 1e-16;
                U[y+ghost_cell][x+ghost_cell][1] = 0;
                U[y+ghost_cell][x+ghost_cell][2] = 0;
                U[y+ghost_cell][x+ghost_cell][3] = 1e-22;
            }
            // if (U[y+ghost_cell][x+ghost_cell][3]<= 1e-16){
            //     U[y+ghost_cell][x+ghost_cell][0] = 1e-16;
            //     U[y+ghost_cell][x+ghost_cell][1] = 0;
            //     U[y+ghost_cell][x+ghost_cell][2] = 0;
            //     U[y+ghost_cell][x+ghost_cell][3] = 0;
            // }

            vel = (U[y][x][1]/U[y][x][0])*(U[y][x][1]/U[y][x][0])+(U[y][x][2]/U[y][x][0])*(U[y][x][2]/U[y][x][0]);
            p = (gamma-1)*(U[y][x][3] - (U[y][x][1]*U[y][x][1] + U[y][x][2]*U[y][x][2])/2/U[y][x][0]);
            // p = std::max(p,float(1e-16));
            // if (vel<0){
            //     cout << "stop" << endl;
            // }
            // if (gamma*p/U[y][x][0]<=0){
            //     cout << "stop" << endl;
            // }
            vel = sqrt(vel);
            this->dt_buffer = std::min(this->dt_buffer, 0.2f*this->dx/(sqrt(gamma*p/U[y][x][0])+vel));
            // if (std::isnan(dt_buffer) || std::isinf(dt_buffer)){
            //     cout << "stop" << endl;
            // }
            
        }
    }

    cout << dt_buffer << endl;

    // initializer->updateBC(U,nx,ny,ghost_cell);

    for (int y = 0; y < ny+2*ghost_cell; y++) {
        for (int x = 0; x < nx+2*ghost_cell; x++) {
            for (int c = 0; c < 4; c++) {
                if (std::isnan(U[y][x][c]) || std::isinf(U[y][x][c])){
                    cout << "stop" << endl;
                }
            }

            if (U[y][x][0]>100){
                cout << "stop" << endl;

            }
        }
    }

    this->dt = this->dt_buffer;
}


float*** CFDSimulation::UPtr(){
    return this->U;
}