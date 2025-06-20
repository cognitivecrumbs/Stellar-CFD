#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <glfw_handler.h>
#include <shader_class.h>
#include <algorithm>    // std::max
#include <cfd_object.h>
#include <initialization.h>

//  HLLE MUSCL Scheme to solve 2D compressible Euler Equations (shock-tube problem)
//
//                                      U_t + F_x + G_y = 0, for (x,y,t) in (0,1)^2x(0,0.3]
//
// The MUSCL scheme is a finite volume method that can provide highly
// accurate numerical solutions for a given system,
// even in cases where the solutions exhibit shocks, discontinuities, or large gradients.
// MUSCL stands for Monotonic Upstream-centered Scheme for Conservation Laws (van Leer, 1979),
// and the term was introduced in a seminal paper by Bram van Leer (van Leer, 1979).
//
// In this code, we integrate in time using RK-2, use the HLLE flux solver, monotonized central (MC) flux limiter (van Leer, 1977)
// Computational Region (x,y)
//
//                                                  *---------*----------*
//                                                  |         |          |
//                                                  |  reg 2  |   reg 1  |
//                                                  |         |          |
//                                                  *---------*----------*
//                                                  |         |          |
//                                                  |  reg 3  |   reg 4  |
//                                                  |         |          |
//                                                  *---------*----------*
//

using namespace std;

// Find maximum of three numbers
double maxs(double a, double b, double c){
    double k;
    k = a;
    if  (b>k){
        k=b;
    }
    if  (c>k){
        k=c;
    }
    return k;
}
// Find minimum of three numbers
double mins(double a, double b, double c){
    double k;
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
int sgn(double v) {
    if (v < 0) return -1;
    if (v > 0) return 1;
    return 0;
}

// Minmod Limiter
double minmod(double a, double b, double c) {
    double vs[] = {abs(a), abs(b), abs(c)};
    double v = mins(vs[0],vs[1],vs[2]);
    double mm;
    double s;
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

template <size_t size_t>
void HLLE(double (&flux)[size_t], double qL[4],double qR[4], double gamma, int normal_x, int normal_y){
    int nx;
    int ny;
    double rL,uL,vL,pL;
    double rR,uR,vR,pR;
    double vnL;
    double vnR;
    double vn;
    double aL,HL,HR;
    double aR,RT;
    double u,v,H,a;

    // normal vectors
    nx = normal_x;
    ny = normal_y;

    // Left State
    rL = qL[0];
    uL = qL[1]/rL;
    vL = qL[2]/rL;
    vnL = uL*nx + vL*ny;
    pL = (gamma-1)*(qL[3] - rL*(uL*uL + vL*vL)/2);
    aL= sqrt((gamma*pL)/rL);
    HL = (qL[3] + pL)/rL;

    // Right State
    rR = qR[0];
    uR = qR[1]/rR;
    vR = qR[2]/rR;
    vnR= uR*nx + vR*ny;
    pR = (gamma-1)*(qR[3] - rR*(uR*uR + vR*vR)/2);
    aR = sqrt(gamma*pR/rR);
    HR = (qR[3] + pR)/rR;

    // First compute the Roe Averages
    RT = sqrt(rR/rL);
    u = (uL+RT*uR)/(1+RT);
    v = (vL+RT*vR)/(1+RT);
    H = (HL+RT*HR)/(1+RT);
    a = sqrt((gamma-1)*(H - (u*u + v*v)/2));
    vn = u*nx + v*ny;

    double SLm;
    double SRp;

    // Wave speed estimates
    SLm = mins((vnL-aL),(vn - a),0.0);
    SRp = maxs((vnR+aR), (vn + a), 0.0);

    if ((SRp-SLm) == 0) {
        printf("DIVIDINGS BY 0: SLm = %lf, SRp = %lf  \n", SLm, SRp);
    }

    // Left and Right fluxes
    double FL[4] = {rL*vnL, rL*vnL*uL + pL*nx, rL*vnL*vL + pL*ny, rL*vnL*HL};
    double FR[4] = {rR*vnR, rR*vnR*uR + pR*nx, rR*vnR*vR + pR*ny, rR*vnR*HR};

    // Compute HLL Flux
    for(int i=0; i<=3; i++) {
        flux[i] = ((SRp)*FL[i] - (SLm) * FR[i] + (SLm)* (SRp)*(qR[i] - qL[i])) / ((SRp) - (SLm));
    }

}

double*** FLUXSOLVER(double*** q, int nx, int ny, double dx, double dy,double gamma) {
    double dqw;
    double dqe;
    double dqc;
    double dqn;
    double dqs;
    double*** res;
    double dqdx[ny][nx][4];
    double dqdy[ny][nx][4];
    double qxL[4];
    double qxR[4];
    double qyL[4];
    double qyR[4];
    double qL[4];
    double qR[4];
    double flux[4];

    res = new double**[nx];
    for (int i=0; i<ny; i++) {
        res[i] = new double*[nx];
    }
    for (int i=0; i<ny; i++) {
        for (int j = 0; j < nx; j++) {
            res[i][j] = new double[4];
        }
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < 4; k++) {
                dqdx[i][j][k] = 0;
                dqdy[i][j][k] = 0;
            }
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
    // Compute and limit slopes at cells (i,j)
    for (int i = 1; i <= ny - 2; i++) {
        for (int j = 1; j <= nx - 2; j++) {
            for (int k = 0; k <= 3; k++) {
                dqw = 2*(q[i][j][k] - q[i][j - 1][k])/dx;
                dqe = 2*(q[i][j + 1][k] - q[i][j][k])/dx;
                dqc = (q[i][j + 1][k] - q[i][j - 1][k]) /(2 * dx);
                dqdx[i][j][k] = minmod(dqw, dqe, dqc);

                dqs = 2 * (q[i][j][k] - q[i - 1][j][k])/ dy;
                dqn = 2 * (q[i + 1][j][k] - q[i][j][k])/ dy;
                dqc = (q[i + 1][j][k] - q[i - 1][j][k])/(2 * dy);
                dqdy[i][j][k] = minmod(dqs, dqn, dqc);
            }
        }
    }

    // Compute residuals x-direction
    // face values
    for (int i = 1; i <= ny - 2; i++) {
        for (int j = 1; j <= nx - 3; j++) {
            for (int k = 0; k <= 3; k++) {
                qxL[k] = q[i][j][k] + (dqdx[i][j][k])*dx/2;
                qxR[k] = q[i][j + 1][k] - (dqdx[i][j + 1][k])*dx/2;
            }
            HLLE(flux, qxL, qxR, gamma, 1, 0);
            for (int k = 0; k <= 3; k++) {
                res[i][j][k] = res[i][j][k] + (flux[k])/ dx;
                res[i][j + 1][k] = res[i][j + 1][k] - (flux[k])/dx;

            }
        }
    }

    // Compute residuals y-direction
    for (int i = 1; i <= ny - 3; i++) {
        for (int j = 1; j <= nx - 2; j++) {
            for (int k = 0; k <= 3; k++) {
                qyL[k] = q[i][j][k] + dqdy[i][j][k]*dy/2;
                qyR[k] = q[i + 1][j][k] - dqdy[i + 1][j][k]*dy/2;
            }
            HLLE(flux, qyL, qyR, gamma, 0, 1);
            for (int k = 0; k <= 3; k++) {
                res[i][j][k] = res[i][j][k] + (flux[k])/dy;
                res[i + 1][j][k] = res[i + 1][j][k] - (flux[k])/dy;
            }
        }
    }
    // Set Boundary Conditions
    for (int j = 1; j <= nx - 2; j++) {
        for (int k = 0; k <= 3; k++) {
            qR[k] = q[ny - 2][j][k] + dqdy[ny - 2][j][k] * dy / 2;
            qL[k] = qR[k];
        }
        HLLE(flux, qL, qR, gamma, 0, 1);
        for (int k = 0; k <= 3; k++) {
            res[ny - 2][j][k] = res[ny - 2][j][k] + (flux[k]) / dy;
        }
    }

    // Flux contribution of the MOST EAST FACE: east face of cell j=ny-2
    for (int i = 1; i <= ny - 2; i++) {
        for (int k = 0; k <= 3; k++) {
            qR[k] = q[i][nx - 2][k] + dqdx[i][nx - 2][k] * dx / 2;
            qL[k] = qR[k];
        }
        HLLE(flux, qL, qR, gamma, 1, 0);
        for (int k = 0; k <= 3; k++) {
            res[i][nx - 2][k] = res[i][nx - 2][k] + (flux[k]) / dx;
        }
    }
    // Flux contribution of the MOST SOUTH FACE: south face of cells j=1.
    for (int j = 1; j <= nx - 2; j++) {
        for (int k = 0; k <= 3; k++) {
            qR[k] = q[1][j][k] - dqdy[1][j][k]*dy/2;
            qL[k] = qR[k];
        }
        HLLE(flux, qL, qR, gamma, 0, -1);
        for (int k = 0; k <= 3; k++) {
            res[1][j][k] = res[1][j][k] + (flux[k])/ dy;
        }
    }
    // Flux contribution of the MOST WEST FACE: west face of cells j=1.
    for (int i = 1; i <= ny - 2; i++) {
        for (int k = 0; k <= 3; k++) {
            qR[k] = q[i][1][k] - dqdx[i][1][k] * dx / 2;
            qL[k] = qR[k];
        }
        HLLE(flux, qL, qR, gamma, -1, 0);
        for (int k = 0; k <= 3; k++) {
            res[i][1][k] = res[i][1][k] + (flux[k]) / dx;
        }
    }
    for (int i = 0; i <=nx-1; i++) {
        for (int j = 0; j <= nx - 1; j++) {
            for (int k = 0; k <= 3; k++) {
                res[i][j][k] = res[i][j][k];
            }
        }
    }
    return res;
}
// int main(int arcv, char** argv){
//     double CFL = 0.40;                  // CFL
//     double tEnd = 3.0;                  // Final time
//     int nx = 240;                       // # of points in x
//     int ny = 240;                       // # number of points in y
//     nx = 10;
//     ny = 10;
//     double gamma = 1.4;                 // Heat Capacity Ratio
//     double Lx = 1;                      // x boundary (left)
//     double Ly = 1;                      // y boundary (top)
//     double dx = Lx/nx;                  // x spatial step
//     double dy = Ly/ny;                  // y spatial step
//     double dt = 0.0011;                 // t time step
//     // double dt = 0.00011;                 // t time step
//     double xc[nx];                      // x array
//     double yc[ny];                      // y array


// // Initialization
//     xc[0]= dx/2;
//     yc[0]= dy/2;
// // Fill x,y array
//     for (int i=1; i<ny; i++){
//         xc[i] = dx/2 + dx*i;
//         yc[i] = dx/2 + dy*i;
//     }
// // Initial Condition
//     double r0[ny][nx];
//     double u0[ny][nx];
//     double v0[ny][nx];
//     double p0[ny][nx];

//     for (int i=0; i<ny; i++){
//         for (int j=0; j<nx; j++){
//             // First Quadrant
//             if(xc[j]>=0.5 && yc[i]>=0.5){
//                 p0[i][j] = 1.5;
//                 r0[i][j] = 1.5;
//                 u0[i][j] = 0.0;
//                 v0[i][j] = 0.0;
//             }
//             // Second Quadrant
//             if(xc[j]<0.5 && yc[i]>=0.5){
//                 p0[i][j] = 0.3;
//                 r0[i][j] = 0.5323;
//                 u0[i][j] = 1.206;
//                 v0[i][j] = 0.0;
//             }
//             // Third Quadrant
//             if(xc[j]<0.5 && yc[i]<0.5){
//                 p0[i][j] = 0.029;
//                 r0[i][j] = 0.138;
//                 u0[i][j] = 1.206;
//                 v0[i][j] = 1.206;
//             }
//             // Fourth Quadrant
//             if(xc[j]>=0.5 && yc[i]<0.5){
//                 p0[i][j] = 0.3;
//                 r0[i][j] = 0.5323;
//                 u0[i][j] = 0.0;
//                 v0[i][j] = 1.206;
//             }
//         }
//     }

//     for (int i=0; i<ny; i++){
//         for (int j=0; j<nx; j++){
//             // First Quadrant
//             if(xc[i]>=0.1){
//                 p0[i][j] = 1.5;
//                 r0[i][j] = 1.5;
//                 p0[i][j] = 2.0;
//                 r0[i][j] = 2.0;
//                 u0[i][j] = 0.0;
//                 v0[i][j] = 0.0;
//             }
//             // Third Quadrant
//             else {
//                 p0[i][j] = 0.029;
//                 r0[i][j] = 0.138;
//                 p0[i][j] = 3.0;
//                 r0[i][j] = 3.;
//                 u0[i][j] = 0.0;
//                 v0[i][j] = 0.0;
//             }
//         }
//     }
// // Initial Total and Internal Energy
//     double E0[ny][nx];
//     double c0[ny][nx];
//     for (int i=0; i<ny; i++){
//         for (int j=0; j<nx; j++) {
//             E0[i][j] = p0[i][j]/((gamma-1)*r0[i][j]) + 0.5*(u0[i][j]*u0[i][j]+ v0[i][j]*v0[i][j]);
//             c0[i][j] = sqrt(gamma*p0[i][j]/r0[i][j]);
//         }
//     }
// // Initial Solution vector Q_0 = [rho_0, rho_0*u_0, rho_0*v_0, rho_0*E_0 ]
//     double Q0[ny][nx][4];
//     for (int i=0; i<ny; i++){

//         for (int j=0; j<nx; j++) {
//             Q0[i][j][0] = r0[i][j];
//             Q0[i][j][1] = r0[i][j]*u0[i][j];
//             Q0[i][j][2] = r0[i][j]*v0[i][j];
//             Q0[i][j][3] = r0[i][j]*E0[i][j];

//         }
//     }
// // Introduce Ghost points
//     nx +=  2;
//     ny += 2;
//     double*** q;
//     q = new double**[nx];
//     for (int i=0; i<ny; i++) {
//         q[i] = new double*[nx];
//     }
//     for (int i=0; i<ny; i++) {
//         for (int j = 0; j < nx; j++) {
//             q[i][j] = new double[4];
//         }
//     }

//     for (int i=0; i<ny-2; i++){
//         for (int j=0; j<nx-2; j++) {
//             q[i+1][j+1][0] = Q0[i][j][0];
//             q[i+1][j+1][1] = Q0[i][j][1];
//             q[i+1][j+1][2] = Q0[i][j][2];
//             q[i+1][j+1][3] = Q0[i][j][3];
//         }
//     }

//     for (int i=0; i<ny; i++){
//         for(int k =0; k<4; k++){
//             q[i][0][k] = q[i][1][k];
//             q[i][nx-1][k] = q[i][nx-2][k];
//             q[0][i][k] = q[1][i][k];
//             q[ny-1][i][k] = q[ny-2][i][k];
//         }
//     }
//     ofstream init("rho_ic.txt");
//     for (int i = 1; i < nx-1 ; i++) {
//         for (int j = 1; j < nx - 1; j++) {
//             init << q[i][j][0] << " ";
//             // cout << q[i][j][0] << " ";
//         }
//         init << "\n";
//     }
//     init.close();

//     ofstream init_p("p_ic.txt");
//     for (int i = 1; i < nx-1 ; i++) {
//         for (int j = 1; j < nx - 1; j++) {
//             init_p << (gamma-1)*q[i][j][0]*((q[i][j][3]/q[i][j][0])-0.5*( pow(q[i][j][1]/q[i][j][0],2) + pow(q[i][j][2]/q[i][j][0],2) ) ) << " ";
//             // cout << q[i][j][0] << " ";
//         }
//         init_p << "\n";
//     }
//     init_p.close();
//     double*** qs;
//     qs = new double**[nx];
//     for (int i=0; i<ny; i++) {
//         qs[i] = new double*[nx];
//     }
//     for (int i=0; i<ny; i++) {
//         for (int j = 0; j < nx; j++) {
//             qs[i][j] = new double[4];
//         }
//     }

//     double*** res;
//     res = new double**[nx];
//     for (int i=0; i<ny; i++) {
//         res[i] = new double*[nx];
//     }
//     for (int i=0; i<ny; i++) {
//         for (int j = 0; j < nx; j++) {
//             res[i][j] = new double[4];
//         }
//     }

//     GLFW_handler glfw_handler = GLFW_handler(800,800);
//     int glfw_handler_return = glfw_handler.init_screen();
//     if (glfw_handler_return == -1){
//         std::cout << "Screen initialization failed" << std::endl;
//         return -1;
//     }

//     GLFWwindow* window = glfw_handler.window;
//     glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

//     Shader shader = Shader("src/include/shaders/vertex_cfd_mesh.glsl","src/include/shaders/frag_cfd_mesh.glsl");

//     unsigned int vao,vbo;

//     glGenVertexArrays(1, &vao);
//     glGenBuffers(1, &vbo);

// // Integrate using RK-2
//     double t = 0;
//     while(t<tEnd && !glfwWindowShouldClose(window)) {

//         vector<float> vertices;
//         vector<float> densities;
//         vector<float> temperature;
//         double max_quant = -1e8;
//         double max_vel = -1e8;
//         double max_u = -1e8;
//         for (int y = 0; y <=nx-1 ; y++) {
//             for (int x = 0; x <= nx-1; x++) {
//                 // for (int c = 0; c <= 3; c++) {
//                     // qs[y][x][c] = q[y][x][c] - dt*res[y][x][c];
//                     vertices.push_back(xc[x]);
//                     vertices.push_back(yc[y]);
//                     densities.push_back(q[y][x][0]);
//                     // temperature.push_back((mesh->temp-50)/70);
//                     // temperature.push_back(mesh->vel_squared*10);
//                     // temperature.push_back(mesh->prim[1]*50);
//                     // temperature.push_back(mesh->U[3]);
//                     double p = (gamma-1)*(q[y][x][3] - (q[y][x][1]*q[y][x][1] + q[y][x][2]*q[y][x][2])/2/q[y][x][0]);
//                     temperature.push_back(p-2);

//                 // }
//                 // if (q[y][x][0]>max_quant){
//                 //     max_quant = q[y][x][0];
//                 // }
//                 max_quant = std::max(max_quant,q[y][x][0]);
//                 double vel = (q[y][x][1]/q[y][x][0])*(q[y][x][1]/q[y][x][0])+(q[y][x][2]/q[y][x][0])+(q[y][x][2]/q[y][x][0]);
//                 vel = sqrt(vel);
//                 max_vel = std::max(max_vel,vel);
//                 max_u = max(max_u,sqrt(gamma*p/q[y][x][0])+vel);
//             }
//         }
//         // dt = CFL*dx/max_u;

//         cout << max_quant << endl;

//         glBindVertexArray(vao);

//         glBindBuffer(GL_ARRAY_BUFFER, vbo);
//         glBufferData(GL_ARRAY_BUFFER,  vertices.size() * sizeof(float), vertices.data(), GL_DYNAMIC_DRAW);
//         glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2* sizeof(float), (void*)0);
//         glEnableVertexAttribArray(0);

//         GLuint densityBuffer;
//         glGenBuffers(1, &densityBuffer);
//         glBindBuffer(GL_ARRAY_BUFFER, densityBuffer);
//         glBufferData(GL_ARRAY_BUFFER, densities.size() * sizeof(float), densities.data(), GL_DYNAMIC_DRAW);
//         glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
//         glEnableVertexAttribArray(1);

//         GLuint temperatureBuffer;
//         glGenBuffers(1, &temperatureBuffer);
//         glBindBuffer(GL_ARRAY_BUFFER, temperatureBuffer);
//         glBufferData(GL_ARRAY_BUFFER, temperature.size() * sizeof(float), temperature.data(), GL_DYNAMIC_DRAW);
//         glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
//         glEnableVertexAttribArray(2);

//         shader.use();
//         glClear(GL_COLOR_BUFFER_BIT);
//         glDrawArrays(GL_POINTS, 0, nx * ny);

//         // Swap buffers and poll events
//         glfwSwapBuffers(window);
//         glfwPollEvents();

//         res = FLUXSOLVER(q, nx, ny, dx, dy, gamma);
//         for (int y = 0; y <=nx-1 ; y++) {
//             for (int x = 0; x <= nx-1; x++) {
//                 for (int c = 0; c <= 3; c++) {
//                     qs[y][x][c] = q[y][x][c] - dt*res[y][x][c];

//                 }
//             }
//         }
//         // reset bcs
//         for (int y = 0; y <= ny-1; y++) {
//             for (int c = 0; c <= 3; c++) {
//                 qs[y][0][c] = qs[y][1][c];
//                 qs[y][nx - 1][c] = qs[y][nx - 2][c];
//             }
//         }
//         for (int x = 0; x <= nx-1; x++) {
//             for (int c = 0; c <= 3; c++) {
//                 qs[0][x][c] = qs[1][x][c];
//                 qs[ny - 1][x][c] = qs[ny - 2][x][c];
//             }
//         }

//         res = FLUXSOLVER(qs, nx, ny, dx, dy, gamma);
//         for (int y = 0; y <= nx-1 ; y++) {
//             for (int x = 0; x <= nx-1; x++) {
//                 for (int c = 0; c <= 3; c++) {
//                     q[y][x][c] = ((q[y][x][c] + qs[y][x][c]- dt*res[y][x][c]))/2;
//                 }
//             }
//         }

//         for (int y = 0; y <= ny-1; y++) {
//             for (int c = 0; c <=3; c++) {
//                 q[y][0][c] = q[y][1][c];
//                 q[y][nx - 1][c] = q[y][nx - 2][c];
//             }
//         }
//         for (int x= 0; x <= nx-1; x++) {
//             for (int c = 0; c <=3; c++) {
//                 q[0][x][c] = q[1][x][c];
//                 q[ny - 1][x][c] = q[ny - 2][x][c];
//             }
//         }

//         if (t + dt > tEnd) {
//             dt = tEnd - t;
//         }
//         t = t + dt;
//         printf("%lf \n",t);
//     }

//     ofstream r("rho.txt");
//     for (int i = 1; i < nx-1 ; i++) {
//         for (int j = 1; j < nx - 1; j++) {
//             r << q[i][j][0] << " ";
//             // cout << q[i][j][0] << " ";
//         }
//         r << "\n";
//     }
//     r.close();
//     ofstream p("p.txt");
//     for (int i = 1; i < nx-1 ; i++) {
//         for (int j = 1; j < nx - 1; j++) {
//             p << (gamma-1)*q[i][j][0]*((q[i][j][3]/q[i][j][0])-0.5*( pow(q[i][j][1]/q[i][j][0],2) + pow(q[i][j][2]/q[i][j][0],2) ) ) << " ";
//             // cout << q[i][j][0] << " ";
//         }
//         p << "\n";
//     }
//     p.close();
//     return 0;
// }

int main(int arcv, char** argv){
    double CFL = 0.40;                  // CFL
    double tEnd = 100.0;                  // Final time
    int nx = 240;                       // # of points in x
    int ny = 240;                       // # number of points in y
    // nx = 3;
    // ny = 3;
    nx = 128;
    ny = 128;
    nx = 64;
    ny = 64;
    ny = 64;
    // nx = 4;
    // ny = 4;
    // nx = 256;
    // ny = 256;
    // nx = 512;
    // ny = 512;
    double gamma = 1.4;                 // Heat Capacity Ratio
    double Lx = 1;                      // x boundary (left)
    double Ly = 1;                      // y boundary (top)
    // Ly = 6;
    double dx = Lx/nx;                  // x spatial step
    double dy = Ly/ny;                  // y spatial step
    double dt = 0.0011;                 // t time step
    // double dt = 0.00011;                 // t time step
    // double xc[nx];                      // x array
    // double yc[ny];                      // y array

    int ghost_cell = 2;

    // CFDSimulation simulation = CFDSimulation(new KelvinHelholtzFull(),nx,ny,ghost_cell,dx,dy,gamma);
    // CFDSimulation simulation = CFDSimulation(new ShockTube(),nx,ny,ghost_cell,dx,dy,gamma);
    CFDSimulation simulation = CFDSimulation(new GravPressureFull(),nx,ny,ghost_cell,dx,dy,gamma);
    // CFDSimulation simulation = CFDSimulation(new ConvectionCurrent(),nx,ny,ghost_cell,dx,dy,gamma);
    // CFDSimulation simulation = CFDSimulation(new ShearStress(),nx,ny,ghost_cell,dx,dy,gamma);

    float*** q = simulation.UPtr();

    GLFW_handler glfw_handler = GLFW_handler(800,800);
    int glfw_handler_return = glfw_handler.init_screen();
    if (glfw_handler_return == -1){
        std::cout << "Screen initialization failed" << std::endl;
        return -1;
    }

    GLFWwindow* window = glfw_handler.window;
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    Shader shader = Shader("src/include/shaders/vertex_cfd_mesh.glsl","src/include/shaders/frag_cfd_mesh.glsl");

    unsigned int vao,vbo;
    GLuint densityBuffer;
    GLuint temperatureBuffer;

    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glGenBuffers(1, &densityBuffer);
    glGenBuffers(1, &temperatureBuffer);

    float last_time = glfwGetTime();
    float last_update_time = glfwGetTime();
    float accumulated_time = 0.f;

    std::chrono::steady_clock::time_point lastTime;
    int frameCount = 0;
    float fps = 0.0f;
    float u,v;

    // Integrate using RK-2
    vector<float> vertices;
    vector<float> densities;
    vector<float> temperature;
    vertices.resize(2*(nx+2*ghost_cell)*(ny+2*ghost_cell),0);
    densities.resize((nx+2*ghost_cell)*(ny+2*ghost_cell),0);
    temperature.resize((nx+2*ghost_cell)*(ny+2*ghost_cell),0);
    double t = 0;
    float min_T = -1e20, max_T = 1e20;
    float min_T_buffer = 0, max_T_buffer = 1e20;
    float T;
    float scale, offset;
    while(t<tEnd && !glfwWindowShouldClose(window)) {
        // float min_phi = 1e20;
        // float max_phi = -1e20;
        scale = 1/(max_T-min_T);

        for (int y = 0; y < ny+2*ghost_cell; y++) {
            for (int x = 0; x < nx+2*ghost_cell; x++) {
                // for (int c = 0; c <= 3; c++) {
                    // qs[y][x][c] = q[y][x][c] - dt*res[y][x][c];
                    // vertices.push_back(simulation.xc[x]);
                    // vertices.push_back(simulation.yc[y]);
                    // densities.push_back(q[y][x][0]);
                    // temperature.push_back(q[y][x][0]-1);

                    vertices[2*(nx+2*ghost_cell)*y+2*x] = simulation.xc[x];
                    vertices[2*(nx+2*ghost_cell)*y+2*x+1] = simulation.yc[y];
                    densities[(nx+2*ghost_cell)*y+x] = q[y][x][0];

                    double p = (gamma-1)*(q[y][x][3] - (q[y][x][1]*q[y][x][1] + q[y][x][2]*q[y][x][2])/2/q[y][x][0]);
                    // temperature.push_back(p);
                    u = q[y][x][1]/q[y][x][0];
                    v = q[y][x][2]/q[y][x][0];

                    T = q[y][x][3]/q[y][x][0]-(u*u+v*v)/2;

                    temperature[(nx+2*ghost_cell)*y+x] = (T-min_T)*scale;
                    max_T_buffer = max(max_T_buffer,T);
                    min_T_buffer = min(min_T_buffer,T);

                    // plot u velocity
                    
                    // temperature[(nx+2*ghost_cell)*y+x] = (u-min_T)*scale;
                    // max_T_buffer = max(max_T_buffer,u);
                    // min_T_buffer = min(min_T_buffer,u);

            }
        }
        max_T = max_T_buffer;
        min_T = min_T_buffer;
        max_T_buffer = -1e20;
        min_T_buffer = 1e20;

        // cout << max_T << endl;

        glBindVertexArray(vao);

        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER,  vertices.size() * sizeof(float), vertices.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2* sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);


        glBindBuffer(GL_ARRAY_BUFFER, densityBuffer);
        glBufferData(GL_ARRAY_BUFFER, densities.size() * sizeof(float), densities.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
        glEnableVertexAttribArray(1);

        glBindBuffer(GL_ARRAY_BUFFER, temperatureBuffer);
        glBufferData(GL_ARRAY_BUFFER, temperature.size() * sizeof(float), temperature.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
        glEnableVertexAttribArray(2);

        shader.use();
        glClear(GL_COLOR_BUFFER_BIT);
        glDrawArrays(GL_POINTS, 0, (nx+2*ghost_cell) * (ny+2*ghost_cell));



        // cout << max_quant << endl;
        // Update FPS once per second
        float current_time = glfwGetTime();
        float dt_screen = current_time - last_time;
        last_time = current_time;
        accumulated_time += dt_screen;
        
        float update_dt = current_time - last_update_time;
        last_update_time = current_time;
        accumulated_time = 0.f;


        frameCount++;
        auto currentTime = std::chrono::steady_clock::now();
        float timeInterval = std::chrono::duration<float>(currentTime - lastTime).count();

        if (timeInterval >= 1.0f) {
            fps = frameCount / timeInterval;
            lastTime = currentTime;
            frameCount = 0;

            // Update the window title with the current FPS
            std::ostringstream title;
            title << "FPS: " << std::fixed << std::setprecision(2) << fps;
            glfwSetWindowTitle(window, title.str().c_str());
        }


        // Swap buffers and poll events
        glfwSwapBuffers(window);
        glfwPollEvents();
        
        simulation.update();
    }

}