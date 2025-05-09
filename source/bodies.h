#include <vector>
#include <cmath>
#include <bits/stdc++.h>
#include <random>
#include <algorithm>
#include <string> 
#include <fstream>
#include <sstream>
#include <map>
#include <iostream>
#include <format>

/* debugging function */

void clearLinePrint(int row, std::string line);
void clearLinePrint(int row, std::string line, bool newLine);

namespace fs = std::filesystem;

float getRandomFloat(float min, float max);

std::vector<float> generateRandomFloats(float min, float max, size_t count);

struct Particle{
    float x;
    float y;
    float vx;
    float vy;
    float mass;
    int red;
    int green;
    int blue;
    float visualRadius;

    int id;

    Particle(): x(0.0f), y(0.0f), vx(0.0f), vy(0.0f), mass(1.0f), red(255), blue(255), green(255), visualRadius(1.0f) {};

private:
};

struct Cluster{
/*
    mean: _x, _y
        center of the cluster where it is the closest to each point inside the cluster
    cm (centerOfMass):  _x, _y
        center of mass of the cluster, used for calculating
        gravity exerted on other clusters
        also made an equation to determine the acceleration without a specific particle internally
    velocity: _x, _y
        velocity the cluster experiences due to other clusters
    totalMass: 
        total amount of mass in the cluster for intercluster gravitational calculations
        change values when adding or removing a particle
    particles:
        data structure that holds all of the particles within this specific cluster
    inertia:
        calculates the inertia by summing up the distances squared of all particles in the cluster
    addParticle:
        adds a particle to this cluster as a membership
    hasParticle:
        checks if a particle is a member of cluster
    remParticle:
        removes the particle if it is a member of this cluster
*/
    static std::map<int, Particle>* simulationParticles;

    float mean_x;
    float mean_y;
    float cm_x;
    float cm_y;
    float vx;
    float vy;
    float totalMass;
    float radiusSqr;
    std::map<int, float> particle_ID_membership;
    
    float inertia();
//    void addParticle(Particle& refParticle);
    bool hasParticle(Particle& refParticle);
//    void remParticle(Particle& refParticle);
    bool inBounds(Particle& refParticle);

    Cluster(): mean_x(0.0f), mean_y(0.0f), cm_x(0.0f), cm_y(0.0f), vx(0.0f), vy(0.0f), totalMass(1.0f), radiusSqr(0.0f), particle_ID_membership(std::map<int, float>()) {};
private:
};

enum class ClusteringType{
    custom_0 = 0,
    k_means = 1,
    fuzzy_k_means = 2
};

struct SimulationParams{
/*
    particleCount:
        Amount of particles for the simulation to start with
    clusterCount:
        amount of clusters to generate
    clusterStarts:
        amount of different iterations/ intial means to test for convergence in cluster initializations
    steps:
        Amount of frames/steps the simulation will run for, # of images saved
    timeScale:
        Amount multiplied into the velocity per step
    fps:
        Frames per second, or steps per second for the output video
    gravity:
        Gravitational constant for the simulation
    mass_lower:
        lower bound for mass range, for random generation
    mass_upper:
        upper bound for mass range, for random generation
    width:
        width of the simulation projection bounds for the Particle Position
    height:
        height of the simulation projection bounds for the Particle Position
    bounds:
        scaling of the bounds, < 1 in screen, > 1 off screen, = 1 perfect
    resolutionWidth:
        the resolution width of the output video in pixels
    resolutionHeight:
        the resolution height of the output video in pixels
    pixelScale:
        the scale factor for scaling particles from projection to frames
        used when rendering, and saving frames, also for bounds checking
        using an 'effective' radius for the particle bounding box
    maxForce:
        max force allowed, 
        use to calculate max velocity based off mass
*/
    ClusteringType clusterAlgorithm = ClusteringType::k_means;
    int particleCount = 500;
    int clusterCount = 4;
    int clusterStarts = 25;
    int steps = 300;
    float timeScale = 1.0f;
    int fps = 30;
    float gravity = 1.08*std::pow(10,1);
    float mass_lower = 1.0f;
    float mass_upper = 1.0f;
    float width = 1600.0f;
    float height = 900.0f;
    float bounds = 1.0f;
    int resolutionWidth = 1920;
    int resolutionHeight = 1080;
    float pixelScale = 1.2f;
    float maxForce = 1000.0f;
    float softening = 10.0f;
};

class ParticleSimulation{
/*
    particles:
        dictionary of particles that are in the simulation <particle.id, Particle>
        id should match index order created,
    clusters:
        dictionary of clusters that are in the simulation
    parameters:
        simulation parameters, steps, particle/cluster count, look at class comments.

    km_membership:
        vector[particle][cluster], denotes true membership, hard clustering, k means
    fkm_uig: 
        vector[particle][cluster], denotes membership probability, soft clustering, fuzzy k means

    addParticle:
        (int): adds a random [amount] of particles to the simualtion,
                random positions and mass based of parameters
        (float): adds a particle of [mass] in a random positions,
                based off the parameters
        (float, float, float): adds a particle of [mass] in [x, y] position,

    createClusters:
        creates the clusters based of clusterCount in parameters
    optimizeClusters_kmeans:
        optimizes the clusters using the k-means algorithm
    optimizeClusters_FKM:
        optimizes the clusters using fuzzy k-means algorithm

    createDirectories:
        checks and ensures the working directory is setup for rendering/saving
    runSim:
        loop/function that runs the simulation, logic, order, control
    calculateForcesCluster:
        calculates the forces between all clusters,
        setting their velocities
    calculateForcesParticle:
        calculates the forces between all particles INTERNAL to a cluster,
        setting their velocities
    stepCluster:
        steps using timeScale as a factor and each clusters' velocity
    stepParticle:
        steps using timeScale as a factor and each clusters' velocity

    renderStep:
        render all particles visually, might use openGL to use CPU and GPU
    saveFrame:
        saves the rendered frame that used renderStep
    makeFilename:
        checks the working directory, and creates a file name based of the parameters
    compileVideo:
        uses the saved frames to compile a video using FFMPEG
    clearFrames:
        deletes the frames in the working directory
*/
    public:
    std::map<int, Particle> particles;
    std::map<int, Cluster> clusters;
    SimulationParams parameters;

    ParticleSimulation(SimulationParams params): parameters(params){};

    void addParticle(int amountOf);
    void addParticle(float mass);
    void addParticle(float mass, float posi_x, float posi_y);

    void createClusters();
    void optimizeClusters_custom(std::map<int, Cluster>& clusterMap);
    void optimizeClusters_kmeans(std::map<int, Cluster>& clusterMap);
    void optimizeClusters_FKM(std::map<int, Cluster>& clusterMap);

    void createDirectories();
    void runSim();
    void calculateForcesCluster();
    void calculateForcesParticle();
    void stepCluster();
    void stepParticle();
    
    void renderStep();
    void saveFrame(int frameNumber);
    std::string getFilenameID();
    void compileVideo(std::string fileID);
    void saveParameters(std::string fileID);
    void loadParameters(std::string fileID);
    void clearFrames();

    private:
};