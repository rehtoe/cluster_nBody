#include "source/n_vector.h"
#include "source/bodies.h"
#include <iostream>

int main(){
/*
    nVec pos1(2), pos2(4), pos3;
    pos1[0] = 1;
    pos1[1] = 2;
    
    pos2[0] = 3;
    pos2[1] = 4;
    pos2[2] = 5;
    pos2[3] = 6;


    pos3 = pos1 + pos2;
    std::cout << "+ Vector pos3 unlocked: ";
    for(int i = 0; i < pos3.size(); i++){
        std::cout << pos3[i] << ", ";
    }
    std::cout << std::endl;

    pos1.lock();
    pos3 = pos1 + pos2;
    std::cout << "+ Vector pos3 locked: ";
    for(int i = 0; i < pos3.size(); i++){
        std::cout << pos3[i] << ", ";
    }
    std::cout << std::endl;

    pos3 = pos1 - pos2;
    std::cout << "- Vector pos3: ";
    for(int i = 0; i < pos3.size(); i++){
        std::cout << pos3[i] << ", ";
    }
    std::cout << std::endl;
    
    pos3 = pos1 * 2;
    std::cout << "* Vector pos3: ";
    for(int i = 0; i < pos3.size(); i++){
        std::cout << pos3[i] << ", ";
    }
    std::cout << std::endl;

    pos3 = pos2 / 2;
    std::cout << "/ Vector pos3: ";
    for(int i = 0; i < pos3.size(); i++){
        std::cout << pos3[i] << ", ";
    }
    std::cout << std::endl;
*/
    SimulationParams _params;
    _params.clusterAlgorithm = ClusteringType::k_means;
    _params.particleCount = 500;
    _params.clusterCount = 4;
    _params.clusterStarts = 25;
    _params.steps = 300;
    _params.timeScale = 1.0f;
    _params.fps = 30;
    _params.gravity = 1.08*std::pow(10,1);
    _params.mass_lower = 1.0f;
    _params.mass_upper = 1.0f;
    _params.width = 1600.0f;
    _params.height = 900.0f;
    _params.bounds = 1.0f;
    _params.resolutionWidth = 1920;
    _params.resolutionHeight = 1080;
    _params.pixelScale = 1.2f;
    _params.maxForce = 1000.0f;
    _params.softening = 10.0f;

    ParticleSimulation sim(_params);

    sim.createDirectories();
    sim.addParticle(sim.parameters.particleCount);
    sim.createClusters();
    switch(sim.parameters.clusterAlgorithm){
        default:
        case ClusteringType::custom_0:
            sim.optimizeClusters_custom(sim.clusters); break;
        case ClusteringType::k_means:
            sim.optimizeClusters_kmeans(sim.clusters); break;
        case ClusteringType::fuzzy_k_means:
            sim.optimizeClusters_FKM(sim.clusters); break;
    }

    std::map<int, int> pID_cG;
    for(auto [g, cluster]:sim.clusters){
        for(auto [id, particle]:cluster.particle_ID_membership){
            pID_cG[id] = g;
        }
    }
    int iii = 0;
    for(auto [id, g]:pID_cG){ 
        if(iii++ == 49){ std::cout << g << std::endl; iii = 0; }
        else{ std::cout << g << ", "; }
    }

    return 0;
}