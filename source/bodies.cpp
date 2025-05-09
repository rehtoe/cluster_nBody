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
    std::vector<float> posX = generateRandomFloats(0, parameters.width*parameters.bounds, amountOf);
    std::vector<float> posY = generateRandomFloats(0, parameters.height*parameters.bounds, amountOf);
    std::vector<float> masses = generateRandomFloats(parameters.mass_lower, parameters.mass_upper, amountOf);
    int next_id = particles.size();
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
    int next_id = particles.size();
    Particle ph_part;
    ph_part.x = getRandomFloat(0, parameters.width*parameters.bounds);
    ph_part.y = getRandomFloat(0, parameters.height*parameters.bounds);
    ph_part.vx = 0;
    ph_part.vy = 0;
    ph_part.mass = mass;
    ph_part.id = next_id;
    particles[next_id] = ph_part;
}
void ParticleSimulation::addParticle(float mass, float posi_x, float posi_y) {
    int next_id = particles.size();
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
        do clusterStarts amount of different starts
    */
    Cluster::simulationParticles = &particles;
    std::map<int,Cluster> ph_clusters;
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
            ph_clusters[i]= ph_cluster;
        }
        if(parameters.clusterAlgorithm == ClusteringType::custom_0){ clearLinePrint(1, "Custom", true); optimizeClusters_custom(ph_clusters); }
        else if(parameters.clusterAlgorithm == ClusteringType::k_means){ clearLinePrint(1, "K-Means", true); optimizeClusters_kmeans(ph_clusters); }
        else if(parameters.clusterAlgorithm == ClusteringType::fuzzy_k_means){ clearLinePrint(1, "Fuzzy-K-Means", true); optimizeClusters_FKM(ph_clusters); }
        
        /*  setting first converged clusters as clusters
            comparing new iteration to current saved configuration
            to check if it is more optimized
        */
        if(starts == 0){
            for(auto &[g, cl]:ph_clusters){
                clusters[g] = cl;
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
    //            for(auto cl_vec:clusters){ for(auto cl_ptr:cl_vec.particle_ID_membership){ delete cl_ptr; } }
                clusters.clear();
                for(auto &[g, cl]:ph_clusters){
                    clusters[g] = cl;
                }
            }
        }

        /*  free up memory of all particles in ph_clusters before the next iteration of the loop 
            Particles in Cluster struct are not owned by the struct, should not delete them.
            this instead clears the vector to have it ready for the next iteration of the loop*/
    //    for(auto cl_vec:ph_clusters){ for(auto cl_ptr:cl_vec.particle_ID_membership){ delete cl_ptr; } }
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
void ParticleSimulation::optimizeClusters_kmeans(std::map<int, Cluster>& clusterMap) {
    /*  has the objective function target change for convergence
        placeholder intial previous value, a current value for when grouping datapoints
        max iterations to run is the integer in the while loop */
    double objective_target = 1.0e-4, objective_prev = 1928374650.0, objective_curr = 0.0;
    int iters = 0;
    while(iters < parameters.clusterStartIterations){
        clearLinePrint(3, "Iteration " + std::to_string(iters + 1), true);
        /*  iterates over every particle in the ParticleSimulation
            makes a float array for distances to each cluster(same order)
            keeps track of the current 'i' and the 'smallest' distance indexes */
        #pragma omp parallel for reduction(+:objective_curr)
        for(auto &[id,parti] : particles){
            float min_dist = 192837465000.0f, ph_dist = 0.0f;
            int closest_cluster = 0;
            /*  calculates the distance from the particle to EVERY cluster mean
                saves the index of the cluster in which the particle is closes to */
            for (int i = 0; i < clusterMap.size(); i++) {
                clearLinePrint(4, "Calculation: Particle ID " + std::to_string(id) + ", Cluster # " + std::to_string(i), true);
                float dx_ = clusterMap[i].mean_x - parti.x;
                float dy_ = clusterMap[i].mean_y - parti.y;
                ph_dist = dx_*dx_ + dy_*dy_;
                if (ph_dist < min_dist) {
                    min_dist = ph_dist;
                    closest_cluster = i;
                }
            }
            //std::cout << std::endl;
            /*  adds the particle to the cluster mean it is closest to
                calculates and adds the objective value(distance squared) of the current iteration setup */
        //    clusterMap[closest_cluster].addParticle(parti);
            clusterMap[closest_cluster].particle_ID_membership[id] = 1.0f;
            objective_curr += min_dist;
            clusterMap[closest_cluster].radiusSqr = std::max(clusterMap[closest_cluster].radiusSqr, min_dist);
        }
        /*  changes the means of the clusters based on the particles it contains
            averages all the x values and y values of all particles
            sets cluster's mean to the same cluster */
        for(auto &[g, clstr]:clusterMap){
            if (clstr.particle_ID_membership.empty()){ continue; }
            float new_x = 0.0, new_y = 0.0;
            for(auto &[id, member]:clstr.particle_ID_membership){
                clearLinePrint(5,"Update: Cluster " + std::to_string(g) + ", Particle ID " + std::to_string(id), true);
                new_x += member*clstr.simulationParticles->at(id).x;
                new_y += member*clstr.simulationParticles->at(id).y;
            }
            clstr.mean_x = new_x/clstr.particle_ID_membership.size();
            clstr.mean_y = new_y/clstr.particle_ID_membership.size();
        }
        std::cout << std::endl;
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
void ParticleSimulation::optimizeClusters_FKM(std::map<int, Cluster>& clusterMap) {}

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
    switch(parameters.clusterAlgorithm){
        default:
        case ClusteringType::custom_0:
            optimizeClusters_custom(clusters); break;
        case ClusteringType::k_means:
            optimizeClusters_kmeans(clusters); break;
        case ClusteringType::fuzzy_k_means:
            optimizeClusters_FKM(clusters); break;
    }

    /*  loop for running the sim
        calculations, steps, rendering/saving frames
    */
    for(int i = 0; i < parameters.steps; i++){
        calculateForcesCluster();
        calculateForcesParticle();
        stepCluster();
        stepParticle();
        renderStep();
        saveFrame(i);
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
    clearFrames();
}
void ParticleSimulation::calculateForcesCluster() {}
void ParticleSimulation::calculateForcesParticle() {}
void ParticleSimulation::stepCluster() {}
void ParticleSimulation::stepParticle() {}

void ParticleSimulation::renderStep() {}
void ParticleSimulation::saveFrame(int frameNumber) {}
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