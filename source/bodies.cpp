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
    for(auto &p:particles){
        float dx = mean_x-p->x;
        float dy = mean_y-p->y;
        sum += dx*dx + dy*dy;
    }
    return sum;
}
void Cluster::addParticle(Particle& refParticle){
    /* uses hasParticle to check if reference exists */
    if(hasParticle(refParticle)){ particles.push_back(&refParticle); }
}
bool Cluster::hasParticle(Particle& refParticle){
    /*  uses std::find to get an iterator to the reference in question
        then return it, if no reference is found it returns vector<>.end() */
    auto it = std::find(particles.begin(), particles.end(), &refParticle) != particles.end();
    return it;
}
void Cluster::remParticle(Particle& refParticle){ 
    /*  uses hasParticle to find iterator
        if particle exists shift it to the end and remove it from the array/vector  */
    if(!hasParticle(refParticle)){ return; }
    else{
        particles.erase(std::remove(particles.begin(), particles.end(), &refParticle), particles.end());
    }
}
bool Cluster::inBounds(Particle& refParticle){
    float dx = mean_x-refParticle.x;
    float dy = mean_y-refParticle.y;
    float distSqr = dx*dx + dy*dy;
    return distSqr <= radiusSqr;
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
                /*  free up memory of all particles in each cluster before replacing clusters 
                    Particles in Cluster struct are not owned by the struct, should not delete them.
                    this instead clears the vector to have it ready for the next iteration of the loop*/
    //            for(auto cl_vec:clusters){ for(auto cl_ptr:cl_vec.particles){ delete cl_ptr; } }
                clusters.clear();
                for(auto cl:ph_clusters){
                    clusters.push_back(cl);
                }
            }
        }

        /*  free up memory of all particles in ph_clusters before the next iteration of the loop 
            Particles in Cluster struct are not owned by the struct, should not delete them.
            this instead clears the vector to have it ready for the next iteration of the loop*/
    //    for(auto cl_vec:ph_clusters){ for(auto cl_ptr:cl_vec.particles){ delete cl_ptr; } }
        ph_clusters.clear();
    } // end for loop
}
void ParticleSimulation::optimizeClusters_kmeans(std::vector<Cluster>& clusterVec) {
    /*  has the objective function target change for convergence
        placeholder intial previous value, a current value for when grouping datapoints
        max iterations to run is the integer in the while loop */
    double objective_target = 1.0e-4, objective_prev = 1928374650.0, objective_curr = 0.0;
    int iters = 0;
    while(iters < 100){
        /*  iterates over every particle in the ParticleSimulation
            makes a float array for distances to each cluster(same order)
            keeps track of the current 'i' and the 'smallest' distance indexes */
        #pragma omp parallel for reduction(+:objective_curr)
        for(auto &parti : particles){
            float min_dist = 192837465000.0f, ph_dist = 0.0f;
            int closest_cluster = 0;
            /*  calculates the distance from the particle to EVERY cluster mean
                saves the index of the cluster in which the particle is closes to */
            for (int i = 0; i < clusterVec.size(); i++) {
                float dx_ = clusterVec[i].mean_x - parti.x;
                float dy_ = clusterVec[i].mean_y - parti.y;
                ph_dist = dx_*dx_ + dy_*dy_;
                if (ph_dist < min_dist) {
                    min_dist = ph_dist;
                    closest_cluster = i;
                }
            }
            /*  adds the particle to the cluster mean it is closest to
                calculates and adds the objective value(distance squared) of the current iteration setup */
            clusterVec[closest_cluster].addParticle(parti);
            objective_curr += min_dist;
            clusterVec[closest_cluster].radiusSqr = std::max(clusterVec[closest_cluster].radiusSqr, min_dist);
        }
        /*  changes the means of the clusters based on the particles it contains
            averages all the x values and y values of all particles
            sets cluster's mean to the same cluster */
        for(auto &clstr:clusterVec){
            if (clstr.particles.empty()){ continue; }
            float new_x = 0.0, new_y = 0.0;
            for(auto cl_parti:clstr.particles){ new_x += cl_parti->x; new_y += cl_parti->y; }
            clstr.mean_x = new_x/clstr.particles.size();
            clstr.mean_y = new_y/clstr.particles.size();
        }
        /*  uses calculated objective value to compare it to the previous
            when under the threshhold of objective_target breaks out of the loop
            if it has not passed the theshold then set the previous value to the recent value
            and reset the next value for the next optimization.
            iterations goes up to count how many */
        if(abs(objective_curr-objective_prev) < objective_target){ break; }
        else{ objective_prev = objective_curr; objective_curr = 0.0; }
        iters++;
    }
}
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