#ifndef MESH_H
#define MESH_H

#include <glm/glm.hpp>
#include <GLFW/glfw3.h>
#include <particles.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <thread>
#include <future>
#include <thread_pool.h> // Hypothetical or third-party thread pool header
#include <Eigen/Sparse>
#include <Eigen/Dense>
// #include <initialization.h>
// #include <algorithm>
// #include <initialization.h>

using namespace std;
using namespace Eigen;

class CFDInitializer;

class CFDMesh
{
public:
    CFDMesh(int ind, glm::vec2 center, float dx, Vector4f U);
    // CFDMesh(int ind, glm::vec2 center, float dx, vector <float> U);
    void update_mesh();

    glm::vec2 center;
    // cell paramters
    // Vector2f dx;
    // cartesianal square meshes
    float dx;
    float vol;
    int ind;

    // cell center 
    // glm::vec2 rho_v;
    float temp;
    // float R = 1;
    // float mu=0.1;
    float mu = 0.;
    float k = 0.1;
    // float p;
    // float v;
    // float e;
    float vel_squared;
    // float cv=1;
    float gamma = 1.4;

    // face paramters stored right as 0, up as 1, left as 2, down as 3
    vector <CFDMesh*> neighbor;
    // std::vector <CFDMesh*> neighbor;
    // vector <Vector4f> face_vals;
    // vector <Vector4f> flux;
    // vector <VectorXf> face_vals_prim;
    // vector <Eigen::Matrix<float, 5, 1>> face_vals_prim;

    // Matrix4f face_vals;
    // Matrix4f flux;
    // MatrixXf face_vals_prim;


    void cons2prim();
    void prim2cons();
    

    // VectorXf prim;
    Eigen::Matrix<float, 5, 1> prim;
    Vector4f U;
    Vector4f dUdt;

    // vector <float> U{4};
    // vector <float> dUdt{4};


    void computeGradient();
    void computeFaceVals();

    void computeFlux();
    // void computeFlux2();
    void computeFluxCentralDifference();
    void computeSource();
    // void computeDE();
    // void computeFluxMUSCL();
    void computeFluxHLLE();
    // void computeFluxSecondOrder();
    void minmodSlope();
    void computeHLLFlux();

    // Vector4f U1;
    // Vector4f k1,k2;

    float fx=0,fy=0;
    void calcSlope();
private:

    // Vector4f U_L,U_R;
    float minmod(float a, float b, float c);
    Vector4f J; // source
    Vector4f flux;

    // vector <float> J{4}; // source
    // vector <float> flux{4};

    // Vector4f slope_L,slope_R,slope_C;

    // 0,0 dudx, 1,0 dudy, 0,1 dvdx, 1,1 dvdy
    // Matrix2f gradient;
    // vector <Matrix2f> face_gradient;


    void computeConvFlux();
    void computeViscFlux();


    Vector4f computeFluxForFace(const Vector4f& U, float p, float vel, float dx, int face_orientation);
    Vector4f computeGradient(int i) {return Vector4f::Zero();};
    Vector4f slope_x, slope_y;
    void HLLE(Vector4f& flux, const Vector4f& qL, const Vector4f& qR, int nx, int ny);
    void Rusanov(Vector4f& flux, const Vector4f& qL, const Vector4f& qR, int nx, int ny);
    void bodyForce();
    void reconstructStates(Vector4f& U_L, Vector4f& U_R, const Vector4f& grad) const;

    // vector <float> computeFluxForFace(const vector<float> & U, float p, float vel, float dx, int face_orientation);
    // vector <float> slope_x{4}, slope_y{4};

    // void HLLE(vector<float> & flux, const vector<float> & qL, const vector<float>& qR, int nx, int ny);
    // void Rusanov(vector<float>& flux, const vector<float>& qL, const vector<float>& qR, int nx, int ny);
    // void bodyForce();


    // // Vector4f computeRiemannFlux(const Vector4f& U_L, const Vector4f& U_R, int face_orientation) const;
    // void reconstructStates(vector <float> & U_L, vector <float> & U_R, const vector<float> & grad) const;

};


class CFDMeshJKL
{
public:
    CFDMeshJKL(int ind, glm::vec2 center, float dx, Vector4f U);
    // CFDMesh(int ind, glm::vec2 center, float dx, vector <float> U);
    void update_mesh();

    glm::vec2 center;
    // cell paramters
    // Vector2f dx;
    // cartesianal square meshes
    float dx;
    float vol;
    int ind;

    // cell center 
    // glm::vec2 rho_v;
    float temp;
    // float R = 1;
    // float mu=0.1;
    float mu = 0.;
    float k = 0.1;
    // float p;
    // float v;
    // float e;
    float vel_squared;
    // float cv=1;
    float gamma = 1.4;

    // face paramters stored right as 0, up as 1, left as 2, down as 3
    vector <CFDMesh*> neighbor;
    // std::vector <CFDMesh*> neighbor;
    // vector <Vector4f> face_vals;
    // vector <Vector4f> flux;
    // vector <VectorXf> face_vals_prim;
    // vector <Eigen::Matrix<float, 5, 1>> face_vals_prim;

    // Matrix4f face_vals;
    // Matrix4f flux;
    // MatrixXf face_vals_prim;


    void cons2prim();
    void prim2cons();
    

    // VectorXf prim;
    Eigen::Matrix<float, 5, 1> prim;
    Vector4f U;
    Vector4f dUdt;

    // vector <float> U{4};
    // vector <float> dUdt{4};


    void computeGradient();
    void computeFaceVals();

    void computeFlux();
    // void computeFlux2();
    void computeFluxCentralDifference();
    void computeSource();
    // void computeDE();
    // void computeFluxMUSCL();
    void computeFluxHLLE();
    // void computeFluxSecondOrder();
    void minmodSlope();
    void computeHLLFlux();

    // Vector4f U1;
    // Vector4f k1,k2;

    float fx=0,fy=0;
    void calcSlope();
private:

    // Vector4f U_L,U_R;
    float minmod(float a, float b, float c);
    Vector4f J; // source
    Vector4f flux;

    // vector <float> J{4}; // source
    // vector <float> flux{4};

    // Vector4f slope_L,slope_R,slope_C;

    // 0,0 dudx, 1,0 dudy, 0,1 dvdx, 1,1 dvdy
    // Matrix2f gradient;
    // vector <Matrix2f> face_gradient;


    void computeConvFlux();
    void computeViscFlux();


    Vector4f computeFluxForFace(const Vector4f& U, float p, float vel, float dx, int face_orientation);
    Vector4f computeGradient(int i) {return Vector4f::Zero();};
    Vector4f slope_x, slope_y;
    void HLLE(Vector4f& flux, const Vector4f& qL, const Vector4f& qR, int nx, int ny);
    void Rusanov(Vector4f& flux, const Vector4f& qL, const Vector4f& qR, int nx, int ny);
    void bodyForce();
    void reconstructStates(Vector4f& U_L, Vector4f& U_R, const Vector4f& grad) const;
};

// class CFDMesh
// {
// public:
//     CFDMesh(int ind, glm::vec2 center, float dx, std::vector<float> U);
//     void update_mesh();

//     glm::vec2 center;
//     float dx;
//     float vol;
//     int ind;
//     glm::vec2 rho_v;
//     float temp;
//     float mu = 0.;
//     float k = 0.1;
//     float vel_squared;
//     float gamma = 1.4;

//     std::vector<CFDMesh*> neighbor;
//     std::vector<float> prim; // Replaces Eigen::Matrix<float, 5, 1>
//     std::vector<float> U;    // Replaces Vector4f
//     std::vector<float> dUdt; // Replaces Vector4f
//     std::vector<float> U_temp; 
//     std::vector<float> k1, k2; 

//     void cons2prim();
//     void prim2cons();
//     void computeGradient();
//     void computeFaceVals();
//     void computeFlux();
//     void computeFlux2();
//     void computeFluxCentralDifference();
//     void computeSource();
//     void computeDE();
//     void computeFluxMUSCL();
//     void computeFluxHLLE();
//     void computeFluxSecondOrder();
//     void minmodSlope();
//     void computeHLLFlux();
//     float minmod(float a, float b, float c);
    
//     std::vector<float> slope_x, slope_y;
// };




class NbodyMesh
{
public:
    // list of ptrs to particles within mesh index
    std::vector<Particle*> particles;

    // coord of mesh location 
    std::vector<glm::vec2> corners;
    glm::vec2 center;
    glm::vec2 com;
    float char_len;

    float tot_mass;
    float tot_energy;
    float tot_ang_momentum;
    float tot_momentum;

    float width;
    float height;

    int active_particles;


    NbodyMesh(std::vector<Particle*> particles,float width, float height, glm::vec2 center, NbodyMesh* parent_mesh);
    void mesh_particle_updater(Particle* particle, int direction);

    // updates mesh and particles pointers 
    virtual void update_mesh();

    // ~NbodyMesh();
    std::vector<NbodyMesh*> sub_mesh;
    NbodyMesh* parent_mesh;
    virtual void remove_particle(int ind);
    void propagate_mesh_tree_states();
    void get_mesh_tree_states();
    virtual void update_particles();

protected:

private:
    // for changing mesh definition (barnes hut)
    void update_mesh_states();
    virtual void modify_mesh(){};
    void reset_mesh_tree_states();

};


class BarnesHutMesh: public NbodyMesh{
    public:
    BarnesHutMesh(std::vector<Particle*> particles,float width,float height,glm::vec2 center,std::vector<int> active_ind,int active);
    // void update_particles() override;
    std::vector<BarnesHutMesh*> sub_mesh;
    void update_mesh() override;

    std::vector<int> active_ind;
    void reinitialize(std::vector<Particle*> particles,float width, float height, glm::vec2 center, std::vector<int> active_ind,int active);

    void remove_particle(int ind) override;
    // std::vector<NbodyMesh*> sub_mesh;
    int npart_active;

    glm::vec2 com;
    
    void return_corners(std::vector<std::vector<glm::vec2>>* vertices,std::vector<glm::vec2>* com);

    private:
    void modify_mesh() override;
    void splitMesh();
    void checkActive();
    void deleteTree();;
    // void garbageFunc(){this->npart_active=this->npart_active;};
    // BarnesHutMesh* top_mesh;
    // std::vector<BarnesHutMesh*> bot_mesh;
};


class BarnesHutCFDMesh{
public:
    BarnesHutCFDMesh(std::vector<CFDMesh*> cfd_meshes,float width,float height,glm::vec2 center,std::vector<int> active_ind,int active,BarnesHutCFDMesh* parent_mesh);
    // void update_particles() override;
    std::vector<BarnesHutCFDMesh*> sub_mesh;
    std::vector<CFDMesh*> cfd_meshes;

    // coord of mesh location 
    std::vector<glm::vec2> corners;
    glm::vec2 center;
    glm::vec2 com;
    float char_len;

    float tot_mass;

    float width;
    float height;

    // int active_particles;
    void update_mesh_mass();
    void reset_mesh_tree_states();

    std::vector<int> active_ind;
    // void reinitialize(std::vector<CFDMesh*> cfd_meshes,float width, float height, glm::vec2 center, std::vector<int> active_ind,int active);

    void remove_particle(int ind);
    // std::vector<NbodyMesh*> sub_mesh;
    int npart_active;
    
    void return_corners(std::vector<std::vector<glm::vec2>>* vertices,std::vector<glm::vec2>* com);
    BarnesHutCFDMesh* parent_mesh;
    void initializeMesh();
    void applyGrav(CFDMesh* mesh);
private:
    void splitMesh();
    void checkActive();
    void deleteTree();;

};

// takes mesh (which contains particles) and performs physics on it (leapfrog integration)
// includes handlers for boundary conditions
class Mesh:public NbodyMesh
{
public:
    float width, height;


    Mesh(std::vector<Particle*> particles,int nmesh_hrz,int nmesh_vrt,float width,float height,glm::vec2 center);

    // // updates mesh and particles pointers 
    // void update_mesh();


};


class CFDMeshWrapper{
public:
    CFDMeshWrapper(int nx, int ny, int ghost_cell, float dx);
    void update_mesh();
    vector <CFDMesh*> mesh;
    BarnesHutCFDMesh* grav_mesh;

    // CFDMesh* CFDMesh;
private:
    CFDInitializer* initializer;
    int nx,ny,ghost_cell;

};

#endif