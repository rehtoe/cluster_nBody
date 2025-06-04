#include "bodies.h"

/*  debugging function to make it easier 
    function goes to a certain row in the terminal
    clears the printed line at said row
    outputs message passed as 'line' parameter
*/
void clearLinePrint(int row, std::string line){
    std::string clearLineStart = "\033[" + std::to_string(row) + "H\033[2K\r";
    std::cout << clearLineStart << line;
}
void clearLinePrint(int row, std::string line, bool newLine){
    std::string clearLineStart = "\033[" + std::to_string(row) + "H\033[2K\r";
    std::cout << clearLineStart << line;
    if(newLine){ std::cout << std::endl; }
}

/* Random Number Generators */

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

std::map<int, Particle>* Cluster::simulationParticles = nullptr;

float Cluster::inertia(){
    float sum = 0.0f;
    for(auto &[id, member]:particle_ID_membership){
        float dx = mean_x-member*simulationParticles->at(id).x;
        float dy = mean_y-member*simulationParticles->at(id).y;
        sum += dx*dx + dy*dy;
    }
    return sum;
}
bool Cluster::hasParticle(Particle& refParticle){
    /*  uses std::find to get an iterator to the reference in question
        then return it, if no reference is found it returns vector<>.end() */
    return particle_ID_membership.count(refParticle.id);
}
bool Cluster::inBounds(Particle& refParticle){
    float dx = mean_x-refParticle.x;
    float dy = mean_y-refParticle.y;
    float distSqr = dx*dx + dy*dy;
    return distSqr <= radiusSqr;
}

/* ParticleSimulation */

void ParticleSimulation::addParticle(int amountOf) {
    std::vector<float> posX = generateRandomFloats(0, parameters.width, amountOf);
    std::vector<float> posY = generateRandomFloats(0, parameters.height, amountOf);
    std::vector<float> masses = generateRandomFloats(parameters.mass_lower, parameters.mass_upper, amountOf);
    int next_id = getNextID();
    for(int i = 0; i < amountOf; i++){
        Particle ph_part;
        ph_part.x = posX[i];
        ph_part.y = posY[i];
        ph_part.vx = 0;
        ph_part.vy = 0;
        ph_part.mass = masses[i];
        ph_part.id = next_id + i;
        particles[ph_part.id] = ph_part;
    }
}
void ParticleSimulation::addParticle(float mass) {
    int next_id = getNextID();
    Particle ph_part;
    ph_part.x = getRandomFloat(0, parameters.width);
    ph_part.y = getRandomFloat(0, parameters.height);
    ph_part.vx = 0;
    ph_part.vy = 0;
    ph_part.mass = mass;
    ph_part.id = next_id;
    particles[next_id] = ph_part;
}
void ParticleSimulation::addParticle(float mass, float posi_x, float posi_y) {
    int next_id = getNextID();
    Particle ph_part;
    ph_part.x = posi_x;
    ph_part.y = posi_y;
    ph_part.vx = 0;
    ph_part.vy = 0;
    ph_part.mass = mass;
    ph_part.id = next_id;
    particles[next_id] = ph_part;
}

void ParticleSimulation::createClusters() {
    /*  set the static reference to the particles map.
        create placeholder cluster map
        objective values to know which start is more optimized
        do clusterStarts amount of different starts
    */
    Cluster::simulationParticles = &particles;
    std::map<int,Cluster> ph_clusters;
    float objective_current = 0.0, objective_previous = std::numeric_limits<float>::max();
    for(int starts = 0; starts < parameters.clusterStarts; starts++){
        clearLinePrint(2, "Start " + std::to_string(starts + 1), true);
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
            if (ph_ints[i] == ph_ints[i-1]) {
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
            ph_clusters[i] = ph_cluster;
        }
        if(parameters.clusterAlgorithm == ClusteringType::custom_0){ clearLinePrint(1, "Custom", true); optimizeClusters_custom(ph_clusters); }
        else if(parameters.clusterAlgorithm == ClusteringType::k_means){ clearLinePrint(1, "K-Means", true); objective_current = optimizeClusters_kmeans(ph_clusters); }
        else if(parameters.clusterAlgorithm == ClusteringType::fuzzy_k_means){ clearLinePrint(1, "Fuzzy-K-Means", true); optimizeClusters_FKM(ph_clusters); }
        
        /*  setting first converged clusters as clusters
            comparing new iteration to current saved configuration
            to check if it is more optimized
        */
        if(objective_current < objective_previous){
            /*  free up memory of all particles in each cluster before replacing clusters 
            Particles in Cluster struct are not owned by the struct, should not delete them.
            this instead clears the vector to have it ready for the next iteration of the loop*/
            clusters.clear();
            for(auto &[g, cl]:ph_clusters){
                clusters[g] = cl;
            }
            objective_previous = objective_current;
            objective_current = 0;
        }

        /*  free up memory of all particles in ph_clusters before the next iteration of the loop 
            Particles in Cluster struct are not owned by the struct, should not delete them.
            this instead clears the vector to have it ready for the next iteration of the loop*/
        std::cout << "Objective Current: " << objective_current << std::endl;
        std::cout << "Objective Previous: " << objective_previous << std::endl;
        std::cout << ph_clusters.size() << " clusters in ph_clusters map<int, Cluster>\n";
        for(auto &[g, cluster]:ph_clusters){
            std::cout << g << " Cluster, (X,Y) = (" << cluster.mean_x << ","<< cluster.mean_y << ")\n";
            std::cout << g << " Cluster, Particle Membership count: " << cluster.particle_ID_membership.size() << std::endl;
        }
        ph_clusters.clear();
    } // end for loop
}
void ParticleSimulation::optimizeClusters_custom(std::map<int, Cluster>& clusterMap) {
    /*  yet to finish code */

    /*  objective function values
        sets up an iteration variable and loop for optimization
    */
    double objective_target = 1.0e-4, objective_prev = 1928374650.0, objective_curr = 0.0;
    int iters = 0;
    while(iters < parameters.clusterStartIterations){

    }

    /* previou progress/ ideas */
    for(auto &[partID, parti]:particles){
        float objectiveValues[clusters.size()];
        for(int g = 0; g < clusters.size(); g++){
            float objSum = 0.0f;
            float weight_x = 1 - (std::abs(parti.x - clusters[g].mean_x)/(parameters.width*parameters.bounds));
            float weight_y = 1 - (std::abs(parti.y - clusters[g].mean_x)/(parameters.height*parameters.bounds));
            float weight_vx = 1 - (std::abs(parti.vx - clusters[g].vx)/(parameters.width/parameters.timeScale));
            float weight_vy = 1 - (std::abs(parti.vy - clusters[g].vy)/(parameters.height/parameters.timeScale));
            float weight_dist = weight_x*weight_y;
            float weight_velo = weight_vx*weight_vy;
            objSum += weight_dist*(parti.x-clusters[g].mean_x)*(parti.x-clusters[g].mean_x);
            objSum += weight_dist*(parti.y-clusters[g].mean_y)*(parti.y-clusters[g].mean_y);
            objSum += weight_velo*(parti.vx-clusters[g].vx)*(parti.vx-clusters[g].vx);
            objSum += weight_velo*(parti.vy-clusters[g].vy)*(parti.vy-clusters[g].vy);
            objectiveValues[g] = objSum;
        }
        /*  averages memberships to be between 0-1 as a percent.
            each particle membership to each cluster is in float objectiveValues[]
        */
        float objValSum = 0.0f;
        for(auto val:objectiveValues){ objValSum += val; }
        for(auto &val:objectiveValues){ val /= objValSum; }
    }
}
float ParticleSimulation::optimizeClusters_kmeans(std::map<int, Cluster>& clusterMap) {
    /*  the objective function of kmeans optimizing(minimizing) the joint distance function
        _target -> threshold for objective function
        _prev -> previous objective function value, used for threshold checks
        _curr -> current objective function value, used for threshold checks
        iters -> iterations for calculations
    */
    float objective_target = 1.0e-2f;
    float objective_prev = std::numeric_limits<float>::max();
    float objective_curr = 0.0f, objective_value = 0.0f;
    int iters = 0;
    std::map<int, Cluster> ph_clustermap = clusterMap;

    /*  loops optimization for parameters.clusterStartIterations times */
    while (iters < parameters.clusterStartIterations && abs(objective_prev - objective_curr) > objective_target) {
        clearLinePrint(3, "Iteration " + std::to_string(iters + 1), true);
        /*  Clear current memberships, doesnt change cluster means */
        for (auto& [g, cluster] : clusterMap) {
            cluster.particle_ID_membership.clear();
        }

        /*  Assign particles to nearest cluster, uses kmeans joint distance function 
            for each particle-cluster relation */
        for (const auto &[id, particle]:particles) {
            float min_distance = std::numeric_limits<float>::max();
            int best_cluster = -1;

            /*  Find closest cluster, smallest JDF value for this particle */
            for (const auto &[g, cluster]:clusterMap) {
                clearLinePrint(4, "Calculating JDF... Particle " + std::to_string(id) + ", Cluster " + std::to_string(g), true);
                /*  Calculate JDF (position and velocity variables) */
                float dx = particle.x - cluster.mean_x;
                float dy = particle.y - cluster.mean_y;
                float dvx = particle.vx - cluster.vx;
                float dvy = particle.vy - cluster.vy;
                /*  JDF squared for n'th particle and g'th cluster */
                float distance_sqr = dx*dx + dy*dy + dvx*dvx + dvy*dvy;
                // float distance = std::sqrt(distance_sqr); 
                /*  update closest cluster index by comparing JDF values */
                if (distance_sqr < min_distance) {
                    min_distance = distance_sqr;
                    best_cluster = g;
                }
            }

            // Assign particle membership to cluster
            clusterMap[best_cluster].particle_ID_membership[id] = 1.0f;
            objective_curr += min_distance;
        }

        // Update cluster mean position and velocity
        for (auto &[g, cluster]:clusterMap) {
            // Keep current means if no particles assigned
            if (cluster.particle_ID_membership.empty()) { continue; }
            // intialize reference variables to find clusters' mean_x, mean_y, vx, vy
            float sum_x = 0.0f, sum_y = 0.0f;
            float sum_vx = 0.0f, sum_vy = 0.0f;
            int count = 0;
            // sum particles in clusters based off their membership(k means, 0 or 1 so no 'variable*member' need)
            for (const auto& [id, member] : cluster.particle_ID_membership) {
                // Only consider members with membership > 0
                if (member > 0.0f) {
                    sum_x += particles.at(id).x;
                    sum_y += particles.at(id).y;
                    sum_vx += particles.at(id).vx;
                    sum_vy += particles.at(id).vy;
                    count++;
                }
            }

            if (count > 0) {
                cluster.mean_x = sum_x / count;
                cluster.mean_y = sum_y / count;
                cluster.vx = sum_vx / count;
                cluster.vy = sum_vy / count;
            }
        }
        // update objective function values for new means
        objective_curr = 0.0f;
        for(const auto &[id,particle]:particles){
            for(const auto &[g,cluster]:clusterMap){
                float dx = particle.x - cluster.mean_x;
                float dy = particle.y - cluster.mean_y;
                float dvx = particle.vx - cluster.vx;
                float dvy = particle.vy - cluster.vy;
                // distance function squared
                float distance_sqr = dx*dx + dy*dy + dvx*dvx + dvy*dvy;
                objective_curr += distance_sqr;
            }
        }
        /*  if threshold is not crossed
            setup objective placeholder values for next iteration */
        objective_value = objective_curr;
        objective_curr = 0.0;
        objective_prev = objective_value;

        iters++;
    }
    return objective_value;
}
void ParticleSimulation::optimizeClusters_FKM(std::map<int, Cluster>& clusterMap) {}
void ParticleSimulation::calculateClusterMass(std::map<int, Cluster>& clusterMap){
    for(auto &[g,cluster]:clusters){
        cluster.totalMass = 0;
        for (auto &[id, membership]:cluster.particle_ID_membership){
            cluster.totalMass += particles.at(id).mass*membership;
            cluster.cm_x += particles.at(id).x*particles.at(id).mass*membership;
            cluster.cm_y += particles.at(id).y*particles.at(id).mass*membership;
        }
        cluster.cm_x /= cluster.totalMass;
        cluster.cm_y /= cluster.totalMass;
    }
}

float ParticleSimulation::objective_kmeans(std::map<int, Cluster>& clusterMap){
    float objective_value = 0.0;
    for(const auto &[id,particle]:particles){
        for(const auto &[g,cluster]:clusterMap){
            float dx = particle.x - cluster.mean_x;
            float dy = particle.y - cluster.mean_y;
            float dvx = particle.vx - cluster.vx;
            float dvy = particle.vy - cluster.vy;
            // distance function squared
            float distance_sqr = dx*dx + dy*dy + dvx*dvx + dvy*dvy;
            objective_value += distance_sqr;
        }
    }
    return objective_value;
}

void ParticleSimulation::createDirectories() {
    /*  frames:
            subfolder where all the frames rendered are saved
        sims:
            subfolder where all the compiled videos are saved 
    */
    std::filesystem::create_directory("./frames/");
    std::filesystem::create_directory("./sims/");
}
void ParticleSimulation::runSim() {
    /*  Creates the directories in the directory where program was executed from,
        adds all the particles, random positions and mass within parameter ranges
        creates the clusters and tests different starts
            createClusters() also optimizes based off parameters.clusterAlgorithm 
    */
    createDirectories();
    addParticle(parameters.particleCount);
    createClusters();
    /* 
    switch(parameters.clusterAlgorithm){
        default:
        case ClusteringType::custom_0:
            optimizeClusters_custom(clusters); break;
        case ClusteringType::k_means:
            optimizeClusters_kmeans(clusters); break;
        case ClusteringType::fuzzy_k_means:
            optimizeClusters_FKM(clusters); break;
    }   */

    /*  create reference start frame
        loop for running the sim
        calculations, steps, rendering/saving frames
    */
    pythonPlot(-1);
    #pragma omp parallel for schedule(dynamic) 
    for(int i = 0; i < parameters.steps; i++){
        clearLinePrint(1, "Step " + std::to_string(i), true);
        switch(parameters.clusterAlgorithm){
            default:
            case ClusteringType::custom_0:
                optimizeClusters_custom(clusters); break;
            case ClusteringType::k_means:
                optimizeClusters_kmeans(clusters); break;
            case ClusteringType::fuzzy_k_means:
                optimizeClusters_FKM(clusters); break;
        }
        calculateClusterMass(clusters);
        calculateForcesCluster();
        calculateForcesParticle();
        stepCluster();
        stepParticle();
        pythonPlot(i);
        /*renderStep();
        saveFrame(i);*/
    // add a function to optimize mminimally,(not full recalculation)
    }

    /*  creates and gets the next file_id
        compiles the video to fileID.mp4 with fileID
        saves the parameters to fileID_info.txt with the same fileID
        then clears the frames directory to be ready for the next sim
    */
    std::string file_id = getFilenameID();
    compileVideo(file_id);
    saveParameters(file_id);
    //  testing for initial start up, manual call until publish
    //  clearFrames();
}
void ParticleSimulation::calculateForcesCluster() {
    /*  update each cluster's total mass based off membership 
        and change their center of mass */
    for(auto &[g,cluster]:clusters){
        if(cluster.totalMass == 0){
            for (auto &[id, membership]:cluster.particle_ID_membership){
                cluster.totalMass += particles.at(id).mass*membership;
                cluster.cm_x += particles.at(id).x*particles.at(id).mass*membership;
                cluster.cm_y += particles.at(id).y*particles.at(id).mass*membership;
            }
            cluster.cm_x /= cluster.totalMass;
            cluster.cm_y /= cluster.totalMass;
        }
    }
    /*  loop over each cluster to set each of its velocities proportional to forces acting upon it */
    for(auto &[g, cluster]:clusters){
        for(auto &[g2, cluster2]:clusters){
            /*  calculate unique pairs once, 1 2 is same as 2 1 */
            if(g >= g2){ continue; }
            else{ // (g < g2)
                /*  previous velocities before calculation for vf = vo + at */
                float vx01 = cluster.vx;
                float vy01 = cluster.vy;
                float vx02 = cluster2.vx;
                float vy02 = cluster2.vy;
                /*  calculate dx and dy cluster1(curr calc) -> cluster2 
                    calculate distance sqrd */
                float dx = cluster2.cm_x - cluster.cm_x;
                float dy = cluster2.cm_y - cluster.cm_y;
                float dsqrd = dx*dx + dy*dy;
                /*  ratios to know what % of force is x and y */
                float dxratio = dx/dsqrd;
                float dyratio = dy/dsqrd;
                /*  for extremely close centers or particles
                    calculates dist after softening/buffer */
                dsqrd += parameters.softening*parameters.softening;
                float d = std::sqrt(dsqrd);
                /*  force of gravity between both
                    fg = G* m1m2/d^2
                    (Nm^2/kg^2) * (kg^2/m^2) = N = kg*m/s^2 */
                float Fg = (parameters.gravity)*(cluster.totalMass)*(cluster2.totalMass) / (dsqrd);
                /*  calculates % of forces
                    current object is positive/ the origin
                    the other object experiences an equal and opposite(-) force
                    N * (m/m) = N */
                float fx = Fg * dxratio;
                float fy = Fg * dyratio;
                /*  cluster was the current so it was (0,0)
                    cluster 2 has opposite force so - */
                cluster.vx += (fx/cluster.totalMass);
                cluster.vy += (fy/cluster.totalMass);
                cluster2.vx += ((-fx)/cluster2.totalMass);
                cluster2.vy += ((-fy)/cluster2.totalMass);
            }
        }
    }
}
void ParticleSimulation::calculateForcesParticle() {
    /*  update each cluster's total mass based off membership 
        and center of mass */
    for(auto &[g,cluster]:clusters){
        for (auto &[id, membership]:cluster.particle_ID_membership){
            /*  previous values for v0, vf = v0 + at */
            /*  de-biasing clusters by removing the given particle's values from cluster */
            float refTotalMass = cluster.totalMass - (particles.at(id).mass*membership);
            float refCMx = cluster.cm_x - (particles.at(id).mass*particles.at(id).x*membership/cluster.totalMass);
            float refCMy = cluster.cm_y - (particles.at(id).mass*particles.at(id).y*membership/cluster.totalMass);
            refCMx *= cluster.totalMass/refTotalMass;
            refCMy *= cluster.totalMass/refTotalMass;
            
            /*  distances between points particle -> Cluster */
            float dx = refCMx - particles.at(id).x;
            float dy = refCMy - particles.at(id).y;
            float distSqrd = dx*dx + dy*dy;
            /*  for extremely close centers or particles */
            distSqrd += parameters.softening*parameters.softening;
            float dist = std::sqrt(distSqrd);

            /*  force of gravity between particles and de-biased cluster */
            float force = (parameters.gravity)*(refTotalMass)*(particles.at(id).mass) / (distSqrd);
            float fx = force * dx/dist;
            float fy = force * dy/dist;
            particles.at(id).vx += (fx/particles.at(id).mass);
            particles.at(id).vy += (fy/particles.at(id).mass);
        }
    }
}
void ParticleSimulation::stepCluster() {
    for(auto &[g, cluster]:clusters){
        /*  moves center of mass and mean of the cluster a velocity*timestep
            the movement is internal to each cluster
            hence being allowed to move mean and cm 
            multiply by parameters.timeScale since that was NOT done in the calculations */
        float pos_x_change = cluster.vx * parameters.timeScale;
        float pos_y_change = cluster.vy * parameters.timeScale;
        cluster.mean_x += pos_x_change;
        cluster.mean_y += pos_y_change;
        cluster.cm_x += pos_x_change;
        cluster.cm_y += pos_y_change;
        for(auto &[id, membership]:cluster.particle_ID_membership){
            /*  moves all particles by the same amount for the cluster */
            particles.at(id).x += pos_x_change;
            particles.at(id).y += pos_y_change;
        }
    }
}
void ParticleSimulation::stepParticle() {
    for(auto &[g, cluster]:clusters){
        for(auto &[id, membership]:cluster.particle_ID_membership){
            /*  moves all particles by a timestep per its velocity */
            particles.at(id).x = particles.at(id).vx * parameters.timeScale;
            particles.at(id).y = particles.at(id).vy * parameters.timeScale;
        }
    }
}

// void ParticleSimulation::renderStep() {}
// void ParticleSimulation::saveFrame(int frameNumber) {}
void ParticleSimulation::renderStep() {
    float scaleWidth = parameters.resolutionWidth/parameters.width;
    float scaleHeight = parameters.resolutionHeight/parameters.height;
    /*
    for(const auto &[g,Cluster]:clusters){
        for(const auto &[id, member]:Cluster.particle_ID_membership){
            Particle ph_part = particles[id];
            renderer.setPixel(ph_part.x*scaleWidth, ph_part.y*scaleHeight, ph_part.red, ph_part.green, ph_part.blue);
        }
    }
    renderer.render();*/
}

void ParticleSimulation::saveFrame(int frameNumber) {
   
}
std::string ParticleSimulation::getFilenameID() {
    /*  setup an initial #0 simulation 
        for each either [mp4 | txt] checks each one to get highest ID
        sets ph_ID to the next ID to use, returns it.
    */
    int ph_id = 0;
    for (const auto& entry : fs::directory_iterator("./sims/")) {
        if (entry.is_regular_file() && entry.path().extension() == ".txt") {
            std::string s = entry.path().filename();
            size_t start = 0;
            size_t end = s.find('_');
            s = s.substr(start, end);
            if(ph_id <= std::stoi(s)) { ph_id = std::stoi(s)+1; }
        }
    }
    return std::to_string(ph_id); 
}
void ParticleSimulation::compileVideo(std::string fileID) {
    /*  Format: ffmpeg -framerate 30 -start_number 0 -i "pathTo/picture%d.png" -c:v libx264 -pix_fmt yuv420p output.mp4
        -framerate, sets fps
        -start_number, starting number for %d which is numbered index, default is 1
        -i, included files, %d represents a numberin ascending order, cant skip, like 0,1,2, 4
        -c:v, video encoding, libx264 -> H.264 (universal)
        -pix_fmt: format for color encoding, yuv420p, sampled method(not lossless) that is universal to players 
    */
    std::string video_command = (
        "ffmpeg -framerate " + std::to_string(parameters.fps) +
        " -start_number 0 -i ./frames/%d.png -c:v libx264 -pix_fmt yuv420p -y ./sims/" +
        fileID + ".mp4"
    );
    system(video_command.c_str());
}
void ParticleSimulation::saveParameters(std::string fileID) {
    /*  makes a map for the variables
        initializes fileName for the new file
    */
    std::map<std::string, std::string> variable_map;
    std::string fileName = "./sims/" + fileID + "_info.txt";
    /*  maps the SimulationParams to the variable_map
    */
    variable_map["particleCount"] = std::to_string(parameters.particleCount);
    variable_map["clusterCount"] = std::to_string(parameters.clusterCount);
    variable_map["clusterStarts"] = std::to_string(parameters.clusterStarts);
    variable_map["clusterStartIterations"] = std::to_string(parameters.clusterStartIterations);
    variable_map["steps"] = std::to_string(parameters.steps);
    variable_map["timeScale"] = std::to_string(parameters.timeScale);
    variable_map["fps"] = std::to_string(parameters.fps);
    variable_map["gravity"] = std::to_string(parameters.gravity);
    variable_map["mass_lower"] = std::to_string(parameters.mass_lower);
    variable_map["mass_upper"] = std::to_string(parameters.mass_upper);
    variable_map["width"] = std::to_string(parameters.width);
    variable_map["height"] = std::to_string(parameters.height);
    variable_map["bounds"] = std::to_string(parameters.bounds);
    variable_map["resolutionWidth"] = std::to_string(parameters.resolutionWidth);
    variable_map["resolutionHeight"] = std::to_string(parameters.resolutionHeight);
    variable_map["pixelScale"] = std::to_string(parameters.pixelScale);
    variable_map["maxForce"] = std::to_string(parameters.maxForce);
    variable_map["softening"] = std::to_string(parameters.softening);
    /*  creates file @ filename, opens a output file stream
        adds all the mapped variables as "variable_name,value"
        closes file stream
    */
    std::ofstream newFile(fileName);
    for(auto &var_pair:variable_map){
        newFile << var_pair.first + "," + var_pair.second + '\n';
    }
    newFile.close();
}
void ParticleSimulation::loadParameters(std::string fileID) {
    /*  initializes a map usd for the variables
        initializes the fileName to load based off given fileID
    */
    std::map<std::string, std::string> variable_map;
    std::string fileName = "./sims/" + fileID + "_info.txt";
    std::ifstream file(fileName);
    std::string line;
    /*  map out the parameters to make the setting more easily readable
    */
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string key, value;

        if (std::getline(ss, key, ',') && std::getline(ss, value)) {
            variable_map[key] = value;
        }
    }
    /*  setting each parameter to whatever config chosen/parsed previously
    */
    parameters.particleCount = std::stoi(variable_map["particleCount"]);
    parameters.clusterCount = std::stoi(variable_map["clusterCount"]);
    parameters.clusterStarts = std::stoi(variable_map["clusterStarts"]);
    parameters.clusterStartIterations = std::stoi(variable_map["clusterStartIterations"]);
    parameters.steps = std::stoi(variable_map["steps"]);
    parameters.timeScale = std::stof(variable_map["timeScale"]);
    parameters.fps = std::stoi(variable_map["fps"]);
    parameters.gravity = std::stof(variable_map["gravity"]);
    parameters.mass_lower = std::stof(variable_map["mass_lower"]);
    parameters.mass_upper = std::stof(variable_map["mass_upper"]);
    parameters.width = std::stof(variable_map["width"]);
    parameters.height = std::stof(variable_map["height"]);
    parameters.bounds = std::stof(variable_map["bounds"]);
    parameters.resolutionWidth = std::stoi(variable_map["resolutionWidth"]);
    parameters.resolutionHeight = std::stoi(variable_map["resolutionHeight"]);
    parameters.pixelScale = std::stof(variable_map["pixelScale"]);
    parameters.maxForce = std::stof(variable_map["maxForce"]);
    parameters.softening = std::stof(variable_map["softening"]);
}
void ParticleSimulation::clearFrames() {
    /*  Called after 'compileVideo()'
            deletes all rendered frames
            saves space, preps folder for next sim
    */
    for (const auto& entry : std::filesystem::directory_iterator("./frames/")) {
        std::filesystem::remove(entry.path());
    }
}

int ParticleSimulation::getNextID(){
    int ph_id = particles.size();
    while(particles.count(ph_id)){ ph_id++; }
    return ph_id;
}
void ParticleSimulation::pythonPlot(int next_frame_id){
    /* creates a filename
        opens a file stream
        creates imports for matplotlib and the lists. */
    std::string pythonFilename = "clusterPlot_";
    pythonFilename +=  std::to_string(next_frame_id) + ".py";
    std::ofstream py_script(pythonFilename);
    py_script   << "import matplotlib.pyplot as plt\n"
                << "import math\n"
                << "x_list = []\n"
                << "y_list = []\n"
                << "g_list = []\n";
    /* passing the data from C++ to Python lists */
    for(auto [g, cluster]:clusters){
        for(auto [id, index]:cluster.particle_ID_membership){
            py_script   << "x_list.append(" << std::to_string(particles[id].x) << ")\n"
                        << "y_list.append(" << std::to_string(particles[id].y) << ")\n"
                        << "g_list.append(" << std::to_string(g) << ")\n";
        }
        py_script   << "print(len(x_list))\n"
                    << "print(len(y_list))\n"
                    << "print(len(g_list))\n";
    }
    py_script   << "print(x_list[500:])\n"
                << "print(y_list[500:])\n"
                << "print(g_list[500:])\n";
    /* creating the plot in python */
    py_script   << "plt.figure(figsize=(16, 9))\n"
                << "colors = plt.get_cmap('viridis', len(set(g_list)))\n"
                << "plt.scatter(x_list, y_list, c=g_list, cmap=colors, alpha=0.7)\n"
                //<< "plt.colorbar(label='Cluster')\n"  // Optional colorbar
                << "plt.xlabel('X Axis')\n"
                << "plt.ylabel('Y Axis')\n"
                << "plt.xlim(" << 0 << ", " << parameters.width*parameters.bounds << ")\n"
                << "plt.ylim(" << 0 << ", " << parameters.height*parameters.bounds << ")\n"
                << "plt.title('Grouped Scatter Plot')\n"
                << "plt.savefig('frames/"+std::to_string(next_frame_id)+".png')\n"
                << "plt.close()\n"
                << "#plt.show()\n";
    py_script.close();
    /* uses a virtual environment inside the working directory */
    /* doesnt work yet

    std::string pythonVenv= "envPlot";
    std::string venvCommand;
    std::string exitVenv = "deactivate";
    if(!fs::exists(pythonVenv)){
        std::string createVenv = "python -m venv " + pythonVenv;
        std::string enterVenv = "source " + pythonVenv + "/bin/activate";
        std::string upgradePip = "pip install --upgrade pip";
        std::string downloadMPL = "pip install matplotlib";
        venvCommand = createVenv + enterVenv + upgradePip + downloadMPL;
    } else{
        std::string enterVenv = "source " + pythonVenv + "/bin/activate";
        venvCommand = enterVenv;
    }
    */
    /*  remove file if it exists
        run the virtual environment command. */
  //  remove("bin/"+std::to_string(next_frame_id)+".png");
  //  system(venvCommand.c_str());
    
    /* for script: create command, run command */
    std::string command = "python " + pythonFilename;
    system(command.c_str());
    /* leave python virtual environment */
  //  system(exitVenv.c_str());
    /* deleting python file */
    remove(pythonFilename.c_str());
    /* open image made by python file */
  //  system("xdg-open bin/cluster_output.png");

}