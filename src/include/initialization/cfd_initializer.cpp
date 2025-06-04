#include <initialization.h>

void KelvinHelholtzFull::init(float*** &U, float xc[], float yc[], int nx, int ny, int ghost_cell, float gamma){
// Initial Condition
    float r0[ny+2*ghost_cell][nx+2*ghost_cell];
    float u0[ny+2*ghost_cell][nx+2*ghost_cell];
    float v0[ny+2*ghost_cell][nx+2*ghost_cell];
    float p0[ny+2*ghost_cell][nx+2*ghost_cell];


    float sigma = 0.05 / sqrt(2.0);
    for (int i=0; i<ny+2*ghost_cell; i++){
        for (int j=0; j<nx+2*ghost_cell; j++){
            // First Quadrant
            if(yc[i]>=0.25 && yc[i]<=0.75){
            // if(yc[i]<=1.0/240 && xc[j]<200.0f/240*1){
                p0[i][j] = 2.5;
                r0[i][j] = 2.0;
                u0[i][j] = 0.5;
                v0[i][j] = 0.1 * sin(4 * M_PI * xc[j]) * 
                            (exp(-(yc[i] - 0.25) * (yc[i] - 0.25) / (2 * sigma * sigma)) + 
                            exp(-(yc[i] - 0.75) * (yc[i] - 0.75) / (2 * sigma * sigma)));
                // v0[i][j] = 0.0;
            }
            // Third Quadrant
            else {
                p0[i][j] = 2.5;
                r0[i][j] = 1.0;
                u0[i][j] = -0.5;
                v0[i][j] = 0.1 * sin(4 * M_PI * xc[j]) * 
                            (exp(-(yc[i] - 0.25) * (yc[i] - 0.25) / (2 * sigma * sigma)) + 
                            exp(-(yc[i] - 0.75) * (yc[i] - 0.75) / (2 * sigma * sigma)));
                // v0[i][j] = 0.0;
                // p0[i][j] = 1.0;
            }

            // p0[i][j] = 1.0;
            // r0[i][j] = 1.0;
            // u0[i][j] = 0.0;
            // v0[i][j] = 0.0; 
        }
    }
    // Initial Total and Internal Energy
    float E0[ny+2*ghost_cell][nx+2*ghost_cell];
    float c0[ny+2*ghost_cell][nx+2*ghost_cell];
    for (int i=0; i<ny+2*ghost_cell; i++){
        for (int j=0; j<nx+2*ghost_cell; j++) {
            E0[i][j] = p0[i][j]/((gamma-1)*r0[i][j]) + 0.5*(u0[i][j]*u0[i][j]+ v0[i][j]*v0[i][j]);
            c0[i][j] = sqrt(gamma*p0[i][j]/r0[i][j]);
        }
    }
    // Initial Solution vector Q_0 = [rho_0, rho_0*u_0, rho_0*v_0, rho_0*E_0 ]
    for (int i=0; i<ny+2*ghost_cell; i++){

        for (int j=0; j<nx+2*ghost_cell; j++) {
            // cout << r0[i][j] << endl;
            U[i][j][0] = r0[i][j];
            U[i][j][1] = r0[i][j]*u0[i][j];
            U[i][j][2] = r0[i][j]*v0[i][j];
            U[i][j][3] = r0[i][j]*E0[i][j];
        }
    }
}


void GravPressureFull::init(float*** &U, float xc[], float yc[], int nx, int ny, int ghost_cell, float gamma){
// Initial Condition
    float r0[ny+2*ghost_cell][nx+2*ghost_cell];
    float u0[ny+2*ghost_cell][nx+2*ghost_cell];
    float v0[ny+2*ghost_cell][nx+2*ghost_cell];
    float p0[ny+2*ghost_cell][nx+2*ghost_cell];

    float width = 0.02;
    float rho_center = 20;

    float sigma = 0.05 / sqrt(2.0);
    for (int i=0; i<ny+2*ghost_cell; i++){
        for (int j=0; j<nx+2*ghost_cell; j++){
            // if(yc[i]>0.4 && yc[i]<0.6 && xc[j]<0.6 && xc[j]>0.4){
            // if(yc[i]>0.49 && yc[i]<0.51 && xc[j]<0.51 && xc[j]>0.49){
            if(abs(yc[i]-0.35)<width && abs(xc[j]-0.5) < width){
            // if(yc[i]>0.25 && yc[i]<0.45 && xc[j]<0.4 && xc[j]>0.2){
            // if(yc[i]>0.50 && yc[i]<0.80 && xc[j]<0.45 && xc[j]>0.15){
            // if(yc[i]>0.35 && yc[i]<0.65 && xc[j]<0.65 && xc[j]>0.35){

                p0[i][j] = 1.0;
                r0[i][j] = 2.0;
                r0[i][j] = rho_center;
                u0[i][j] = 1.5;
                v0[i][j] = 0.0;
            }
            else if(abs(yc[i]-0.65)<width && abs(xc[j]-0.5) < width){
                p0[i][j] = 1.0;
                r0[i][j] = 2.0;
                r0[i][j] = rho_center;
                u0[i][j] = -1.5;
                v0[i][j] = 0.0;
            }
            else if (abs(xc[j]-0.65)<width && abs(yc[i]-0.5)<width){
                p0[i][j] = 1.0;
                r0[i][j] = 2.0;
                r0[i][j] = rho_center;
                u0[i][j] = 0.0;
                v0[i][j] = 1.5;
            }
            else if (abs(xc[j]-0.35)<width && abs(yc[i]-0.5)<width){
                p0[i][j] = 1.0;
                r0[i][j] = 2.0;
                r0[i][j] = rho_center;
                u0[i][j] = 0.0;
                v0[i][j] = -1.5;
            }
            // Third Quadrant
            else {
                p0[i][j] = 1.0;
                r0[i][j] = 1.0;
                p0[i][j] = 0.2;
                r0[i][j] = 0.4;
                u0[i][j] = 0.0;
                v0[i][j] = 0.0;
            }


            // p0[i][j] = 0.0;
            // r0[i][j] = 1.0;
            // u0[i][j] = 0.0;
            // v0[i][j] = 0.0; 
        }
    }

    // for (int i = 0; i<ny; i++){
    //     p0[i][ghost_cell] = 1e-22;
    //     // p0[i][ghost_cell] = 1;
    //     r0[i][ghost_cell] = 1e-22;
    //     p0[i][ny+ghost_cell-1] = 1e-22;
    //     // p0[i][ny+ghost_cell-1] = 1;
    //     r0[i][ny+ghost_cell-1] = 1e-22;
    // }
    // for (int j = 0; j<nx; j++){
    //     p0[ghost_cell][j] = 1e-22;
    //     // p0[ghost_cell][j] = 1;
    //     r0[ghost_cell][j] = 1e-22;
    //     p0[nx+ghost_cell-1][j] = 1e-22;
    //     // p0[nx+ghost_cell-1][j] = 1;
    //     r0[nx+ghost_cell-1][j] = 1e-22;
    // }


    // Initial Total and Internal Energy
    float E0[ny+2*ghost_cell][nx+2*ghost_cell];
    float c0[ny+2*ghost_cell][nx+2*ghost_cell];
    for (int i=0; i<ny+2*ghost_cell; i++){
        for (int j=0; j<nx+2*ghost_cell; j++) {
            E0[i][j] = p0[i][j]/((gamma-1)*r0[i][j]) + 0.5*(u0[i][j]*u0[i][j]+ v0[i][j]*v0[i][j]);
            c0[i][j] = sqrt(gamma*p0[i][j]/r0[i][j]);
        }
    }
    // Initial Solution vector Q_0 = [rho_0, rho_0*u_0, rho_0*v_0, rho_0*E_0 ]
    for (int i=0; i<ny+2*ghost_cell; i++){

        for (int j=0; j<nx+2*ghost_cell; j++) {
            // cout << r0[i][j] << endl;
            U[i][j][0] = r0[i][j];
            U[i][j][1] = r0[i][j]*u0[i][j];
            U[i][j][2] = r0[i][j]*v0[i][j];
            U[i][j][3] = r0[i][j]*E0[i][j];
        }
    }
}


void ShockTube::init(float*** &U, float xc[], float yc[], int nx, int ny, int ghost_cell, float gamma){
// Initial Condition
    float r0[ny+2*ghost_cell][nx+2*ghost_cell];
    float u0[ny+2*ghost_cell][nx+2*ghost_cell];
    float v0[ny+2*ghost_cell][nx+2*ghost_cell];
    float p0[ny+2*ghost_cell][nx+2*ghost_cell];


    float sigma = 0.05 / sqrt(2.0);
    for (int i=0; i<ny+2*ghost_cell; i++){
        for (int j=0; j<nx+2*ghost_cell; j++){
            // First Quadrant
            if(yc[i]<0.2){
                p0[i][j] = 2.5;
                r0[i][j] = 2.0;
                u0[i][j] = 0.0;
                v0[i][j] = 0.0;
                }
            // Third Quadrant
            else {
                p0[i][j] = 1.0;
                r0[i][j] = 1.0;
                u0[i][j] = 0.0;
                v0[i][j] = 0.0;
            }

            // if(yc[i]>0.3 && yc[i]<0.7 && xc[j]<0.7 && xc[j]>0.3){
            // if(yc[i]>0.2 && yc[i]<0.8 && xc[j]<0.8 && xc[j]>0.2){
            // if (yc[i]>0.2 && yc[i]<0.8){
            //     p0[i][j] = 2.5*1.2;
            //     r0[i][j] = 2.0*2;
            //     p0[i][j] = 0.1;
            //     r0[i][j] = 0.1;
            //     u0[i][j] = 0.0;
            //     v0[i][j] = 0.0;     
            if(yc[i]>0.3 && yc[i]<0.7){
                p0[i][j] = 2.5*1.2;
                r0[i][j] = 2.0*2;
                p0[i][j] = 2.0;
                r0[i][j] = 2.0;
                u0[i][j] = 0.0;
                v0[i][j] = 0.0;
                }
            // }
            // Third Quadrant
            else {
                p0[i][j] = 1e-16;
                r0[i][j] = 1e-16;
                // r0[i][j] = 0.0001;
                // p0[i][j] = 0.0001;
                r0[i][j] = 0.0001;
                p0[i][j] = 0.0001;
                u0[i][j] = 0.0;
                v0[i][j] = 0.0;
            }

        }
    }
    // Initial Total and Internal Energy
    float E0[ny+2*ghost_cell][nx+2*ghost_cell];
    float c0[ny+2*ghost_cell][nx+2*ghost_cell];
    for (int i=0; i<ny+2*ghost_cell; i++){
        for (int j=0; j<nx+2*ghost_cell; j++) {
            E0[i][j] = p0[i][j]/((gamma-1)*r0[i][j]) + 0.5*(u0[i][j]*u0[i][j]+ v0[i][j]*v0[i][j]);
            c0[i][j] = sqrt(gamma*p0[i][j]/r0[i][j]);
        }
    }
    // Initial Solution vector Q_0 = [rho_0, rho_0*u_0, rho_0*v_0, rho_0*E_0 ]
    for (int i=0; i<ny+2*ghost_cell; i++){

        for (int j=0; j<nx+2*ghost_cell; j++) {
            // cout << r0[i][j] << endl;
            U[i][j][0] = r0[i][j];
            U[i][j][1] = r0[i][j]*u0[i][j];
            U[i][j][2] = r0[i][j]*v0[i][j];
            U[i][j][3] = r0[i][j]*E0[i][j];
        }
    }
}

void ConvectionCurrent::init(float*** &U, float xc[], float yc[], int nx, int ny, int ghost_cell, float gamma){
// Initial Condition
    float r0[ny+2*ghost_cell][nx+2*ghost_cell];
    float u0[ny+2*ghost_cell][nx+2*ghost_cell];
    float v0[ny+2*ghost_cell][nx+2*ghost_cell];
    float p0[ny+2*ghost_cell][nx+2*ghost_cell];
    
    float sigma = 0.05 / sqrt(2.0);
    for (int i=0; i<ny+2*ghost_cell; i++){
        for (int j=0; j<nx+2*ghost_cell; j++){

            r0[i][j] = 1;
            p0[i][j] = 100;
            u0[i][j] = 0.0;
            v0[i][j] = 0.0;

        }
    }
    // Initial Total and Internal Energy
    float E0[ny+2*ghost_cell][nx+2*ghost_cell];
    // float c0[ny+2*ghost_cell][nx+2*ghost_cell];
    for (int i=0; i<ny+2*ghost_cell; i++){
        for (int j=0; j<nx+2*ghost_cell; j++) {
            E0[i][j] = p0[i][j]/((gamma-1)*r0[i][j]) + 0.5*(u0[i][j]*u0[i][j]+ v0[i][j]*v0[i][j]);
            // c0[i][j] = sqrt(gamma*p0[i][j]/r0[i][j]);
        }
    }
    // Initial Solution vector Q_0 = [rho_0, rho_0*u_0, rho_0*v_0, rho_0*E_0 ]
    for (int i=0; i<ny+2*ghost_cell; i++){

        for (int j=0; j<nx+2*ghost_cell; j++) {
            // cout << r0[i][j] << endl;
            U[i][j][0] = r0[i][j];
            U[i][j][1] = r0[i][j]*u0[i][j];
            U[i][j][2] = r0[i][j]*v0[i][j];
            U[i][j][3] = r0[i][j]*E0[i][j];
        }
    }
}


void ConvectionCurrent::updateBC(float*** &U, int nx, int ny, int ghost_cell){
    // for (int j = ghost_cell+nx*3/4; j < nx+ghost_cell; j++){
    //     // U[ghost_cell][j][3] += U[ghost_cell][j][0]*1e0;
    //     // U[ghost_cell][j][3] += 4e-1;
    //     U[ghost_cell][j][3] += 1e0;
    //     // U[ghost_cell][j][3] -= U[ghost_cell][j][0]*0;
    // }
    // for (int j = ghost_cell; j < nx+ghost_cell-nx*3/4; j++){
    //     // U[ny+ghost_cell-1][j][3] -= U[ny+ghost_cell][j][0]*1e0;
    //     // U[ny+ghost_cell-1][j][3] -= 4e-1;
    //     U[ny+ghost_cell-1][j][3] -= 1e0;
    //     // U[ny+ghost_cell-1][j][3] += U[ghost_cell][j][0]*0;
    // }

    float k = 1e-2;
    float e;
    float u,v;
    float dx = 1./64;
    float dy = 1./64;
    float dt = 1e-3;
    // for (int j = ghost_cell; j < nx+ghost_cell-nx*3/4; j++){
    for (int j = ghost_cell; j < nx+ghost_cell-nx; j++){
        u = U[ny+ghost_cell-1][j][1]/U[ny+ghost_cell-1][j][0];
        v = U[ny+ghost_cell-1][j][2]/U[ny+ghost_cell-1][j][0];
        e = U[ny+ghost_cell-1][j][3]/U[ny+ghost_cell-1][j][0] - (u*u+v*v)/2;
        U[ny+ghost_cell-1][j][3] += k*(100-e)/dy*dt;

    }

    // for (int j = ghost_cell+nx*3/4; j < nx+ghost_cell; j++){
    for (int j = ghost_cell+nx*1.5/4; j < nx*2.5/4+ghost_cell; j++){
        u = U[ghost_cell][j][1]/U[ghost_cell][j][0];
        v = U[ghost_cell][j][2]/U[ghost_cell][j][0];
        e = U[ghost_cell][j][3]/U[ghost_cell][j][0] - (u*u+v*v)/2;
        U[ghost_cell][j][3] += k*(10000-e)/dy*dt;
    }

    // adiabatic no slip walls

    for (int k=0;k<ghost_cell;k++){
        for (int i=ghost_cell;i<ny+ghost_cell;i++){
            for (int c=0;c<4;c++){
                // set u,v opposite, otherwise scalars equal
                if (c == 1){
                    // left
                    U[i][ghost_cell-k-1][c] = -U[i][ghost_cell+k][c];
                    // right wall
                    U[i][nx+ghost_cell+k][c] = -U[i][nx+ghost_cell-k-1][c];            
                }
                else {
                    // left
                    U[i][ghost_cell-k-1][c] = U[i][ghost_cell+k][c];
                    // right wall
                    U[i][nx+ghost_cell+k][c] = U[i][nx+ghost_cell-k-1][c];    
                }
            }
        }
    }

    // adiabatic no slip walls
    for (int k=0;k<ghost_cell;k++){
        for (int j=ghost_cell;j<nx+ghost_cell;j++){
            for (int c=0;c<4;c++){
                if (c==2){
                    // bottom wall
                    U[ghost_cell-k-1][j][c] = -U[ghost_cell+k][j][c];
                    // top wall
                    U[ghost_cell+ny+k][j][c] = -U[ghost_cell+ny-k-1][j][c];
                }
                else {
                    // bottom wall
                    U[ghost_cell-k-1][j][c] = U[ghost_cell+k][j][c];
                    // top wall
                    U[ghost_cell+ny+k][j][c] = U[ghost_cell+ny-k-1][j][c];
                }
            }
        }
    }

    // periodic handle corners
    for (int i=0;i<ghost_cell;i++){
        for (int j=0;j<ghost_cell;j++){
            for (int k=0;k<4;k++){
                // upper right corner
                U[ghost_cell+ny+i][ghost_cell+nx+j][k] = U[ghost_cell+i][ghost_cell+j][k];
                // lower left corner
                U[i][j][k] = U[ny+i][nx+j][k];
                // upper left corner
                U[ghost_cell+ny+i][j][k] = U[ghost_cell+i][nx+j][k];
                // lower right corner
                U[i][ghost_cell+nx+j][k] = U[ny+i][ghost_cell+j][k];
            }
        }
    }
}


void ShearStress::init(float*** &U, float xc[], float yc[], int nx, int ny, int ghost_cell, float gamma){
// Initial Condition
    float r0[ny+2*ghost_cell][nx+2*ghost_cell];
    float u0[ny+2*ghost_cell][nx+2*ghost_cell];
    float v0[ny+2*ghost_cell][nx+2*ghost_cell];
    float p0[ny+2*ghost_cell][nx+2*ghost_cell];


    float sigma = 0.05 / sqrt(2.0);
    for (int i=0; i<ny+2*ghost_cell; i++){
        for (int j=0; j<nx+2*ghost_cell; j++){
            // First Quadrant
            if(yc[i]>=0.25 && yc[i]<=0.75){
            // if(yc[i]<=1.0/240 && xc[j]<200.0f/240*1){
                p0[i][j] = 1;
                r0[i][j] = 1.0;
                u0[i][j] = 10.0;
                v0[i][j] = 0;
            }
            // Third Quadrant
            else {
                p0[i][j] = 1;
                r0[i][j] = 1.0;
                u0[i][j] = -10;
                v0[i][j] = 0;
            }

            // First Quadrant
            if(xc[j]>=0.25 && xc[j]<=0.75){
            // if(yc[i]<=1.0/240 && xc[j]<200.0f/240*1){
                p0[i][j] = 1;
                r0[i][j] = 1.0;
                u0[i][j] = 0.0;
                v0[i][j] = 10.0;
            }
            // Third Quadrant
            else {
                p0[i][j] = 1;
                r0[i][j] = 1.0;
                u0[i][j] = 0;
                v0[i][j] = -10.0;
            }


        }
    }
    // Initial Total and Internal Energy
    float E0[ny+2*ghost_cell][nx+2*ghost_cell];
    float c0[ny+2*ghost_cell][nx+2*ghost_cell];
    for (int i=0; i<ny+2*ghost_cell; i++){
        for (int j=0; j<nx+2*ghost_cell; j++) {
            E0[i][j] = p0[i][j]/((gamma-1)*r0[i][j]) + 0.5*(u0[i][j]*u0[i][j]+ v0[i][j]*v0[i][j]);
            c0[i][j] = sqrt(gamma*p0[i][j]/r0[i][j]);
        }
    }
    // Initial Solution vector Q_0 = [rho_0, rho_0*u_0, rho_0*v_0, rho_0*E_0 ]
    for (int i=0; i<ny+2*ghost_cell; i++){

        for (int j=0; j<nx+2*ghost_cell; j++) {
            // cout << r0[i][j] << endl;
            U[i][j][0] = r0[i][j];
            U[i][j][1] = r0[i][j]*u0[i][j];
            U[i][j][2] = r0[i][j]*v0[i][j];
            U[i][j][3] = r0[i][j]*E0[i][j];
        }
    }
}

void ShearStress::updateBC(float*** &U, int nx, int ny, int ghost_cell){
    // periodic
    for (int k=0;k<ghost_cell;k++){
        for (int i=ghost_cell;i<ny+ghost_cell;i++){
            for (int c=0;c<4;c++){
                // left
                U[i][ghost_cell-k-1][c] = U[i][nx+ghost_cell-1-k][c];
                // right wall
                U[i][nx+ghost_cell+k][c] = U[i][ghost_cell+k][c];
            }
        }
    }

    for (int k=0;k<ghost_cell;k++){
        for (int j=ghost_cell;j<nx+ghost_cell;j++){
            for (int c=0;c<4;c++){
                // bottom wall
                U[ghost_cell-k-1][j][c] = U[ghost_cell+ny-k-1][j][c];
                // top wall
                U[ghost_cell+ny+k][j][c] = U[ghost_cell+k][j][c];

            }
        }
    }

    // periodic handle corners
    for (int i=0;i<ghost_cell;i++){
        for (int j=0;j<ghost_cell;j++){
            for (int k=0;k<4;k++){
                // upper right corner
                U[ghost_cell+ny+i][ghost_cell+nx+j][k] = U[ghost_cell+i][ghost_cell+j][k];
                // lower left corner
                U[i][j][k] = U[ny+i][nx+j][k];
                // upper left corner
                U[ghost_cell+ny+i][j][k] = U[ghost_cell+i][nx+j][k];
                // lower right corner
                U[i][ghost_cell+nx+j][k] = U[ny+i][ghost_cell+j][k];
            }
        }
    }
}

void KelvinHelholtz::init(std::vector <CFDMesh*> &mesh_list, int nx, int ny, int ghost_cell){
    float sigma = 0.05/sqrt(2);
    for (CFDMesh* mesh: mesh_list){
        if (mesh->center.y>=0.25 && mesh->center.y <= 0.75){
            // mesh->prim[0] = 3.;
            mesh->prim[0] = 1.0 + ((fabs(mesh->center.y - 0.5) < 0.25) ? 1.0 : 0.0);
            mesh->prim[1] = 0.5;
            mesh->prim[2] = sin(mesh->center.x/6.28/10);
            mesh->prim[2] = 0.01 * sin(4 * M_PI * mesh->center.x) * 
                (exp(-(mesh->center.y - 0.25) * (mesh->center.y - 0.25) / (2 * sigma * sigma)) + 
                 exp(-(mesh->center.y - 0.75) * (mesh->center.y - 0.75) / (2 * sigma * sigma)));
            mesh->prim[4] = 2.5;        
        }
        else{
            // mesh->prim[0] = 2.;
            mesh->prim[0] = 1.0 + ((fabs(mesh->center.y - 0.5) < 0.25) ? 1.0 : 0.0);
            mesh->prim[1] = -0.5;
            // mesh->prim[2] = 0;
            mesh->prim[2] = 0.01 * sin(4 * M_PI * mesh->center.x) * 
                (exp(-(mesh->center.y - 0.25) * (mesh->center.y - 0.25) / (2 * sigma * sigma)) + 
                 exp(-(mesh->center.y - 0.75) * (mesh->center.y - 0.75) / (2 * sigma * sigma)));
            mesh->prim[4] = 2.5;        
        }

        // Define constants
        double w0 = 0.1;
        double sigma = 0.05 / sqrt(2.0);

        // Apply conditions to initialize mesh fields
        if (fabs(mesh->center.y - 0.5) < 0.25) {
            // Density
            mesh->prim[0] = 1.0 + 1.0; // Equivalent to (np.abs(Y - 0.5) < 0.25)
            
            // Velocity in x-direction
            mesh->prim[1] = -0.5 + 1.0; // Equivalent to (-0.5 + (np.abs(Y - 0.5) < 0.25))
            
            // Velocity in y-direction
            mesh->prim[2] = w0 * sin(4 * M_PI * mesh->center.x) * 
                            (exp(-(mesh->center.y - 0.25) * (mesh->center.y - 0.25) / (2 * sigma * sigma)) + 
                            exp(-(mesh->center.y - 0.75) * (mesh->center.y - 0.75) / (2 * sigma * sigma)));
        } else {
            // Density
            mesh->prim[0] = 1.0; // Outside the condition, only `1.0` remains
            
            // Velocity in x-direction
            mesh->prim[1] = -0.5; // No additional value added outside the condition
            
            // Velocity in y-direction
            mesh->prim[2] = w0 * sin(4 * M_PI * mesh->center.x) * 
                            (exp(-(mesh->center.y - 0.25) * (mesh->center.y - 0.25) / (2 * sigma * sigma)) + 
                            exp(-(mesh->center.y - 0.75) * (mesh->center.y - 0.75) / (2 * sigma * sigma)));
        }

        // Pressure (uniform for all conditions)
        mesh->prim[4] = 2.5; // Equivalent to (P = 2.5 * np.ones(X.shape))

        // if (mesh->center.y<=0.7 && mesh->center.x <= 0.7 &&mesh->center.y>=0.3 && mesh->center.x >= 0.3 ){
        //     mesh->prim[0] = 3.;
        //     mesh->prim[1] = 0;
        //     mesh->prim[2] = 0;
        //     mesh->prim[4] = 2.;        
        // }
        // else{
        //     mesh->prim[0] = 2.;
        //     mesh->prim[1] = 0;
        //     mesh->prim[2] = 0;
        //     mesh->prim[4] = 2;        
        // }

        // if (mesh->center.y<=0.4 && mesh->center.x <= 0.4 ){
        //     mesh->prim[0] = 3.;
        //     mesh->prim[1] = 0;
        //     mesh->prim[2] = 0;
        //     mesh->prim[4] = 2.;        
        // }
        // else{
        //     mesh->prim[0] = 2.;
        //     mesh->prim[1] = 0;
        //     mesh->prim[2] = 0;
        //     mesh->prim[4] = 2;        
        // }
  

        mesh->prim2cons();

    }

    // // setup periodic bcs
    // for (int i=0;i<ny+2*ghost_cell;i++){
    //     for (int k=0;k<ghost_cell;k++){
    //         mesh_list[i*(nx+2*ghost_cell)+ghost_cell-k-1] = mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell-1-k];
    //         mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell+k] = mesh_list[i*(nx+2*ghost_cell)+ghost_cell+k];
    //         // mesh_list[i*(nx+2*ghost_cell)+ghost_cell-k-1]->U = &mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell-1-k]->U;
    //         // mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell-1-k]->U = &mesh_list[i*(nx+2*ghost_cell)+ghost_cell-k-1]->U;
    //     }
    // }

    // for (int j=0;j<nx+2*ghost_cell;j++){
    //     for (int k=0;k<ghost_cell;k++){
    //         mesh_list[(ghost_cell-1-k)*(nx+2*ghost_cell)+j] = mesh_list[(ghost_cell+ny-1-k)*(nx+2*ghost_cell)+j];
    //         mesh_list[(ghost_cell+ny-1-k)*(nx+2*ghost_cell)+j] = mesh_list[(ghost_cell-1-k)*(nx+2*ghost_cell)+j];
    //     }
    // }
}



void ObliqueShock::init(std::vector <CFDMesh*> &mesh_list, int nx, int ny, int ghost_cell){
    float u=-1,v=-1;
    for (CFDMesh* mesh: mesh_list){
        float diag = abs(mesh->center.x - mesh->center.y);
        if (diag > 0.9055 && diag < 1.1055){
            mesh->prim[0] = 3.;
            mesh->prim[1] = u;
            mesh->prim[2] = v;
            mesh->prim[4] = 3;
        }
        else if (diag < 0.1055){
            mesh->prim[0] = 3.;
            mesh->prim[1] = u;
            mesh->prim[2] = v;
            mesh->prim[4] = 3;
        }
        else if (diag > 1.9055){
            mesh->prim[0] = 3.;
            mesh->prim[1] = u;
            mesh->prim[2] = v;
            mesh->prim[4] = 3;
        }
        else {
            mesh->prim[0] = 2.;
            mesh->prim[1] = u;
            mesh->prim[2] = v;
            mesh->prim[4] = 2;
        }

        mesh->prim2cons();

    }

    // // setup periodic bcs
    // for (int i=0;i<ny;i++){
    //     for (int j=0;j<nx;j++){
    //         for (int k=0;k<4;k++){
    //             if (i==0 && k==3){
    //                 mesh_list[i*nx+j]->neighbor[k] = mesh_list[(nx-1)*nx+j];
    //             }
    //             else if (i==nx-1 && k==1){
    //                 mesh_list[i*nx+j]->neighbor[k] = mesh_list[j];
    //             }
    //             else if (j==0 && k==2){
    //                 mesh_list[i*nx+j]->neighbor[k] = mesh_list[i*nx+nx-1];
    //             }
    //             else if (j==nx-1 && k==0){
    //                 mesh_list[i*nx+j]->neighbor[k] = mesh_list[i*nx];
    //             }
    //         }
    //     }
    // }
}

void GravPressure::init(std::vector <CFDMesh*> &mesh_list, int nx, int ny, int ghost_cell){
    for (CFDMesh* mesh: mesh_list){      
        mesh->prim[0] = 2.;
        mesh->prim[1] = 0;
        mesh->prim[2] = 0;
        mesh->prim[4] = 2;

        mesh->prim2cons();
    }
}


void ObliqueShock::updateBC(std::vector <CFDMesh*> mesh_list, int nx, int ny, int ghost_cell){
    // adiabatic no slip walls
    for (int k=0;k<ghost_cell;k++){
        for (int i=ghost_cell;i<ny+ghost_cell;i++){
            // set u,v opposite, otherwise scalars equal
            // left wall
            mesh_list[i*(nx+2*ghost_cell)+ghost_cell-k-1]->U = mesh_list[i*(nx+2*ghost_cell)+ghost_cell+k]->U;
            mesh_list[i*(nx+2*ghost_cell)+ghost_cell-k-1]->U[1] *= -1;
            // mesh_list[i*(nx+2*ghost_cell)+ghost_cell-k-1]->U[2] *= -1;

            // right wall
            mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell+k]->U = mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell-1-k]->U;
            mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell+k]->U[1] *= -1;
            // mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell+k]->U[2] *= -1;
            
            // mesh_list[i*(nx+2*ghost_cell)+ghost_cell-k-1]->U = mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell-1-k]->U;
            // // right wall
            // mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell+k]->U = mesh_list[i*(nx+2*ghost_cell)+ghost_cell+k]->U;


        }
    }

    // adiabatic no slip walls
    for (int k=0;k<ghost_cell;k++){
        for (int j=ghost_cell;j<nx+ghost_cell;j++){
            // set u,v opposite, otherwise scalars equal
            // bottom wall
            mesh_list[(ghost_cell-k-1)*(nx+2*ghost_cell)+j]->U = mesh_list[(ghost_cell+k)*(nx+2*ghost_cell)+j]->U;
            // mesh_list[(ghost_cell-k-1)*(nx+2*ghost_cell)+j]->U[1] *= -1;
            mesh_list[(ghost_cell-k-1)*(nx+2*ghost_cell)+j]->U[2] *= -1;

            // top wall
            mesh_list[(ghost_cell+ny+k)*(nx+2*ghost_cell)+j]->U = mesh_list[(ghost_cell+ny-k-1)*(nx+2*ghost_cell)+j]->U;
            // mesh_list[(ghost_cell+ny+k)*(nx+2*ghost_cell)+j]->U[1] *= -1;
            mesh_list[(ghost_cell+ny+k)*(nx+2*ghost_cell)+j]->U[2] *= -1; 
        }
    }


    // // periodic
    // for (int k=0;k<ghost_cell;k++){
    //     for (int i=ghost_cell;i<ny+ghost_cell;i++){
    //         // left
    //         mesh_list[i*(nx+2*ghost_cell)+ghost_cell-k-1]->U = mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell-1-k]->U;
    //         // right wall
    //         mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell+k]->U = mesh_list[i*(nx+2*ghost_cell)+ghost_cell+k]->U;
    //     }
    // }

    // for (int k=0;k<ghost_cell;k++){
    //     for (int j=ghost_cell;j<nx+ghost_cell;j++){
    //         // bottom wall
    //         mesh_list[(ghost_cell-k-1)*(nx+2*ghost_cell)+j]->U = mesh_list[(ghost_cell+ny-k-1)*(nx+2*ghost_cell)+j]->U;
    //         // top wall
    //         mesh_list[(ghost_cell+ny+k)*(nx+2*ghost_cell)+j]->U = mesh_list[(ghost_cell+k)*(nx+2*ghost_cell)+j]->U;

    //     }
    // }
}


void KelvinHelholtz::updateBC(std::vector <CFDMesh*> mesh_list, int nx, int ny, int ghost_cell){
    // periodic
    for (int k=0;k<ghost_cell;k++){
        for (int i=ghost_cell;i<ny+ghost_cell;i++){
            // left
            mesh_list[i*(nx+2*ghost_cell)+ghost_cell-k-1]->U = mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell-1-k]->U;
            // right wall
            mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell+k]->U = mesh_list[i*(nx+2*ghost_cell)+ghost_cell+k]->U;
        }
    }

    for (int k=0;k<ghost_cell;k++){
        for (int j=ghost_cell;j<nx+ghost_cell;j++){
            // bottom wall
            mesh_list[(ghost_cell-k-1)*(nx+2*ghost_cell)+j]->U = mesh_list[(ghost_cell+ny-k-1)*(nx+2*ghost_cell)+j]->U;
            // top wall
            mesh_list[(ghost_cell+ny+k)*(nx+2*ghost_cell)+j]->U = mesh_list[(ghost_cell+k)*(nx+2*ghost_cell)+j]->U;

        }
    }

}

void KelvinHelholtzFull::updateBC(float*** &U, int nx, int ny, int ghost_cell){
    // periodic
    for (int k=0;k<ghost_cell;k++){
        for (int i=ghost_cell;i<ny+ghost_cell;i++){
            for (int c=0;c<4;c++){
                // left
                U[i][ghost_cell-k-1][c] = U[i][nx+ghost_cell-1-k][c];
                // right wall
                U[i][nx+ghost_cell+k][c] = U[i][ghost_cell+k][c];
            }
        }
    }

    for (int k=0;k<ghost_cell;k++){
        for (int j=ghost_cell;j<nx+ghost_cell;j++){
            for (int c=0;c<4;c++){
                // bottom wall
                U[ghost_cell-k-1][j][c] = U[ghost_cell+ny-k-1][j][c];
                // top wall
                U[ghost_cell+ny+k][j][c] = U[ghost_cell+k][j][c];

            }
        }
    }
}


void ShockTube::updateBC(float*** &U, int nx, int ny, int ghost_cell){
    // periodic
    for (int k=0;k<ghost_cell;k++){
        for (int i=ghost_cell;i<ny+ghost_cell;i++){
            for (int c=0;c<4;c++){
                // set u,v opposite, otherwise scalars equal
                if (c == 1){
                    // left
                    U[i][ghost_cell-k-1][c] = -U[i][ghost_cell+k][c];
                    // right wall
                    U[i][nx+ghost_cell+k][c] = -U[i][nx+ghost_cell-k-1][c];            
                }
                else {
                    // left
                    U[i][ghost_cell-k-1][c] = U[i][ghost_cell+k][c];
                    // right wall
                    U[i][nx+ghost_cell+k][c] = U[i][nx+ghost_cell-k-1][c];    
                }
            }
        }
    }

    // adiabatic no slip walls
    for (int k=0;k<ghost_cell;k++){
        for (int j=ghost_cell;j<nx+ghost_cell;j++){
            for (int c=0;c<4;c++){
                if (c==2){
                    // bottom wall
                    U[ghost_cell-k-1][j][c] = -U[ghost_cell+k][j][c];
                    // top wall
                    U[ghost_cell+ny+k][j][c] = -U[ghost_cell+ny-k-1][j][c];
                }
                else {
                    // bottom wall
                    U[ghost_cell-k-1][j][c] = U[ghost_cell+k][j][c];
                    // top wall
                    U[ghost_cell+ny+k][j][c] = U[ghost_cell+ny-k-1][j][c];
                }
            }
        }
    }
}

void GravPressure::updateBC(std::vector <CFDMesh*> mesh_list, int nx, int ny, int ghost_cell){

    // adiabatic no slip walls
    for (int k=0;k<ghost_cell;k++){
        for (int i=ghost_cell;i<ny+ghost_cell;i++){
            // set u,v opposite, otherwise scalars equal
            // left wall
            mesh_list[i*(nx+2*ghost_cell)+ghost_cell-k-1]->U = mesh_list[i*(nx+2*ghost_cell)+ghost_cell+k]->U;
            mesh_list[i*(nx+2*ghost_cell)+ghost_cell-k-1]->U[1] *= -1;
            // mesh_list[i*(nx+2*ghost_cell)+ghost_cell-k-1]->U[2] *= -1;

            // right wall
            mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell+k]->U = mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell-1-k]->U;
            mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell+k]->U[1] *= -1;
            // mesh_list[i*(nx+2*ghost_cell)+nx+ghost_cell+k]->U[2] *= -1;
        }
    }
    // adiabatic no slip walls
    for (int k=0;k<ghost_cell;k++){
        for (int j=ghost_cell;j<nx+ghost_cell;j++){
            // set u,v opposite, otherwise scalars equal
            // bottom wall
            mesh_list[(ghost_cell-k-1)*(nx+2*ghost_cell)+j]->U = mesh_list[(ghost_cell+k)*(nx+2*ghost_cell)+j]->U;
            // mesh_list[(ghost_cell-k-1)*(nx+2*ghost_cell)+j]->U[1] *= -1;
            mesh_list[(ghost_cell-k-1)*(nx+2*ghost_cell)+j]->U[2] *= -1;

            // top wall
            mesh_list[(ghost_cell+ny+k)*(nx+2*ghost_cell)+j]->U = mesh_list[(ghost_cell+ny-k-1)*(nx+2*ghost_cell)+j]->U;
            // mesh_list[(ghost_cell+ny+k)*(nx+2*ghost_cell)+j]->U[1] *= -1;
            mesh_list[(ghost_cell+ny+k)*(nx+2*ghost_cell)+j]->U[2] *= -1; 
        }
    }
}


void GravPressureFull::updateBC(float*** &U, int nx, int ny, int ghost_cell){

    // for (int k=0;k<ghost_cell;k++){
    //     for (int i=ghost_cell;i<ny+ghost_cell;i++){
    //         for (int c=0;c<4;c++){
    //             // set u,v opposite, otherwise scalars equal
    //             if (c == 1){
    //                 // left
    //                 U[i][ghost_cell-k-1][c] = -U[i][ghost_cell+k][c];
    //                 // right wall
    //                 U[i][nx+ghost_cell+k][c] = -U[i][nx+ghost_cell-k-1][c];            
    //             }
    //             else {
    //                 // left
    //                 U[i][ghost_cell-k-1][c] = U[i][ghost_cell+k][c];
    //                 // right wall
    //                 U[i][nx+ghost_cell+k][c] = U[i][nx+ghost_cell-k-1][c];    
    //             }
    //         }
    //     }
    // }

    // // adiabatic no slip walls
    // for (int k=0;k<ghost_cell;k++){
    //     for (int j=ghost_cell;j<nx+ghost_cell;j++){
    //         for (int c=0;c<4;c++){
    //             if (c==2){
    //                 // bottom wall
    //                 U[ghost_cell-k-1][j][c] = -U[ghost_cell+k][j][c];
    //                 // top wall
    //                 U[ghost_cell+ny+k][j][c] = -U[ghost_cell+ny-k-1][j][c];
    //             }
    //             else {
    //                 // bottom wall
    //                 U[ghost_cell-k-1][j][c] = U[ghost_cell+k][j][c];
    //                 // top wall
    //                 U[ghost_cell+ny+k][j][c] = U[ghost_cell+ny-k-1][j][c];
    //             }
    //         }
    //     }
    // }

    // periodic
    for (int k=0;k<ghost_cell;k++){
        for (int i=ghost_cell;i<ny+ghost_cell;i++){
            for (int c=0;c<4;c++){
                // left
                U[i][ghost_cell-k-1][c] = U[i][nx+ghost_cell-1-k][c];
                // right wall
                U[i][nx+ghost_cell+k][c] = U[i][ghost_cell+k][c];
            }
        }
    }

    for (int k=0;k<ghost_cell;k++){
        for (int j=ghost_cell;j<nx+ghost_cell;j++){
            for (int c=0;c<4;c++){
                // bottom wall
                U[ghost_cell-k-1][j][c] = U[ghost_cell+ny-k-1][j][c];
                // top wall
                U[ghost_cell+ny+k][j][c] = U[ghost_cell+k][j][c];
            }
        }
    }

    // periodic handle corners
    for (int i=0;i<ghost_cell;i++){
        for (int j=0;j<ghost_cell;j++){
            for (int k=0;k<4;k++){
                // upper right corner
                U[ghost_cell+ny+i][ghost_cell+nx+j][k] = U[ghost_cell+i][ghost_cell+j][k];
                // lower left corner
                U[i][j][k] = U[ny+i][nx+j][k];
                // upper left corner
                U[ghost_cell+ny+i][j][k] = U[ghost_cell+i][nx+j][k];
                // lower right corner
                U[i][ghost_cell+nx+j][k] = U[ny+i][ghost_cell+j][k];
            }
        }
    }
}