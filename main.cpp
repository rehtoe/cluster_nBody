#include "source/bodies.h"
#include <iostream>

int main(){

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



    return 0;
}