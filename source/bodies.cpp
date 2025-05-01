#include "bodies.h"

/* Random generators */

float getRandomFloat(float min, float max) {
    static std::random_device rd;  // Seed (only once)
    static std::mt19937 gen(rd()); // Mersenne Twister engine
    std::uniform_real_distribution<float> dis(min, max);
    return dis(gen);
}

std::vector<float> generateRandomFloats(float min, float max, size_t count) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(min, max);

    std::vector<float> numbers(count);
    for (auto& num : numbers) {
        num = dis(gen);
    }
    return numbers;
}

/* Clusters */
float Cluster::inertia(){
    float sum = 0.0f;
    for(auto &p:particles){ sum += std::pow(mean_x-p->x,2) + std::pow(mean_y-p->y,2); }
    return sum;
}

/* ParticleSimulation */

void ParticleSimulation::addParticle(int amountOf) {
    std::vector<float> posX = generateRandomFloats(parameters.width - parameters.width*parameters.bounds, parameters.width*parameters.bounds, amountOf);
    std::vector<float> posY = generateRandomFloats(parameters.height - parameters.height*parameters.bounds, parameters.height*parameters.bounds, amountOf);
    std::vector<float> masses = generateRandomFloats(parameters.mass_lower, parameters.mass_upper, amountOf);
    for(int i = 0; i < amountOf; i++){
        Particle ph_part;
        ph_part.x = posX[i];
        ph_part.y = posY[i];
        ph_part.vx = 0;
        ph_part.vy = 0;
        ph_part.mass = masses[i];
        particles.push_back(ph_part);
    }
}
void ParticleSimulation::addParticle(float mass) {
    Particle ph_part;
    ph_part.x = getRandomFloat(parameters.width - parameters.width*parameters.bounds, parameters.width*parameters.bounds);
    ph_part.y = getRandomFloat(parameters.height - parameters.height*parameters.bounds, parameters.height*parameters.bounds);
    ph_part.vx = 0;
    ph_part.vy = 0;
    ph_part.mass = mass;
    particles.push_back(ph_part);
}
void ParticleSimulation::addParticle(float mass, float posi_x, float posi_y) {
    Particle ph_part;
    ph_part.x = posi_x;
    ph_part.y = posi_y;
    ph_part.vx = 0;
    ph_part.vy = 0;
    ph_part.mass = mass;
    particles.push_back(ph_part);
}

void ParticleSimulation::createClusters() {
    std::vector<Cluster> ph_clusters(parameters.clusterCount);
    for(int starts = 0; starts < parameters.clusterStarts; starts++){
        /*  use the random floats function to generate potential indexes
            converts them into ints and puts them into a vector
        */
        std::vector<float> ph_floats = generateRandomFloats(0, particles.size()-1, parameters.clusterCount);
        std::vector<int> ph_ints;
        for(auto fl:ph_floats){ ph_ints.push_back(std::floor(fl)); }
        
        /*  sort the int vector is sorted least to greatest,
            then ensure they are unique
        */
        std::sort(ph_ints.begin(), ph_ints.end());
        for(int i = 1; i < ph_ints.size(); i++){
            if (ph_ints[i] <= ph_ints[i-1]) {
                ph_ints[i] = ph_ints[i-1] + 1;
            }
        }

        /*  setting 1 reference initialized cluster means
            optimizing the clusters per ClusteringType    
        */
        for(int i = 0; i < parameters.clusterCount; i++){
            Cluster ph_cluster;
            ph_cluster.mean_x = particles[ph_ints[i]].x; 
            ph_cluster.mean_y = particles[ph_ints[i]].y; 
            ph_clusters.push_back(ph_cluster);
        }
        if(parameters.clusterAlgorithm == ClusteringType::k_means){
            optimizeClusters_kmeans(ph_clusters);
        } else if(parameters.clusterAlgorithm == ClusteringType::fuzzy_k_means){
            optimizeClusters_FKM(ph_clusters);
        }
        
        /*  setting first converged clusters as clusters
            comparing new iteration to current saved configuration
            to check if it is more optimized
        */
        if(starts == 0){
            for(auto cl:ph_clusters){
                clusters.push_back(cl);
            }
        } else {
            float ph_sum_curr = 0.0f;
            float ph_sum = 0.0f;
            for(int i = 0; i < parameters.clusterCount; i++){
                ph_sum_curr += clusters[i].inertia();
                ph_sum += ph_clusters[i].inertia();
            }
            if(ph_sum < ph_sum_curr){
                /* free up memory of all particles in each cluster before replacing clusters */
                for(auto cl_vec:clusters){ for(auto cl_ptr:cl_vec.particles){ delete cl_ptr; } }
                clusters.clear();
                for(auto cl:ph_clusters){
                    clusters.push_back(cl);
                }
            }
        }

        /* free up memory of all particles in ph_clusters before the next iteration of the loop */
        for(auto cl_vec:ph_clusters){ for(auto cl_ptr:cl_vec.particles){ delete cl_ptr; } }
        ph_clusters.clear();
    }
}
void ParticleSimulation::optimizeClusters_kmeans(std::vector<Cluster>& clusterVec) {}
void ParticleSimulation::optimizeClusters_FKM(std::vector<Cluster>& clusterVec) {}

void ParticleSimulation::createDirectories() {}
void ParticleSimulation::runSim() {}
void ParticleSimulation::calculateForcesCluster() {}
void ParticleSimulation::calculateForcesParticle() {}
void ParticleSimulation::stepCluster() {}
void ParticleSimulation::stepParticle() {}

void ParticleSimulation::renderStep() {}
void ParticleSimulation::saveFrame(int frameNumber) {}
std::string ParticleSimulation::makeFilename() { return ""; }
void ParticleSimulation::compileVideo(std::string filename) {}
void ParticleSimulation::clearFrames() {}