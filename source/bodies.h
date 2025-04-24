#include "n_vector.h"

struct Particle{
    nVec position;
    nVec velocity;
    float mass;
    nVec color;
    float visualRadius;

    Particle(): position(nVec(3)), velocity(nVec(3)), mass(1.0f), color(nVec(3)), visualRadius(1.0f) {};
    Particle(int dimensions): position(nVec(dimensions)), velocity(nVec(dimensions)), mass(1.0f), color(nVec(dimensions)), visualRadius(1.0f) {};

private:
};

struct Cluster{
    nVec mean;
    nVec centerOfMass;
    nVec velocity;
    float totalMass;
    std::vector<&Particle> particles;
    
private:
};

class ParticleSimulation{
    
};