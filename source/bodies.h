#include "n_vector.h"

struct Particle{
    nVec position;
    nVec velocity;
    float mass;
    nVec color;
    float visualRadius;

    Particle(): position(nVec(2)), velocity(nVec(2)), mass(1.0f), color(nVec(3)), visualRadius(1.0f) {};
    Particle(int dimensions): position(nVec(dimensions)), velocity(nVec(dimensions)), mass(1.0f), color(nVec(3)), visualRadius(1.0f) {};

private:
};

struct Cluster{
/*
    mean:
        center of the cluster where it is the closest to each point inside the cluster
    centerOfMass:
        center of mass of the cluster, used for calculating
        gravity exerted on other clusters
        also made an equation to determine the acceleration without a specific particle internally
    velocity:
        velocity the cluster experiences due to other clusters
    particles:
        data structure that holds all of the particles within this specific cluster

*/
    nVec mean;
    nVec centerOfMass;
    nVec velocity;
    float totalMass;
    std::vector<Particle&> particles;
    
    Cluster(): position(nVec(2)), centerOfMass(nVec(2)), velocity(nVec(2)), totalMass(1.0f), particles(std::vector<Particle&>()) {};
    Cluster(int dimensions): position(nVec(dimensions)), centerOfMass(nVec(dimensions)), velocity(nVec(dimensions)), totalMass(1.0f), particles(std::vector<Particle&>()) {};
private:
};

struct SimulationParams{
/*
    dimensions:
        Default to 2, want to implement n dimensions, using nVec
    particleCount:
        Amount of particles for the simulation to start with
    clusterCount:
        amount of clusters to generate
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
    int dimensions = 2;
    int particleCount = 500;
    int clusterCount = 4;
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
    float maxForce = 1000.0f
};

class ParticleSimulation{
/*
    particles:
        list of particles that are in the simulation
    clusters:
        list of clusters that are in the simulation
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
        (float, nVec): adds a particle of [mass] in [nvec] position,

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
    std::vector<Particle> particles;
    std::vector<Cluster> clusters;
    SimulationParams parameters;

    std::vector<std::vector<int>> km_membership;
    std::vector<std::vector<float>> fkm_uig;

    void addParticle(int amountOf);
    void addParticle(float mass);
    void addParticle(float mass, nVec position);

    void createClusters();
    void optimizeClusters_kmeans()
    void optimizeClusters_FKM()

    void createDirectories();
    void runSim();
    void calculateForcesCluster();
    void calculateForcesParticle();
    void stepCluster();
    void stepParticle();
    
    void renderStep();
    void saveFrame(int frameNumber);
    std::string makeFilename();
    void compileVideo(std::string filename);
    void clearFrames();

    private:
};