#include "n_vector.h"

float nVec::length() const{
    float sum = 0.0f;
    for(auto val:values)
        sum += val*val;
    return std::sqrt(sum);
}

nVec nVec::normalized() const{ return nVec(values)/size(); }

int nVec::size() const{
    return values.size();
}

bool nVec::isLocked(){ return locked; }

void nVec::lock(){ locked = true; }

void nVec::unlock(){ locked = false; }

nVec nVec::operator+(nVec& other){
    int phSize = lock_state_size(other);
    std::vector<float> phVec(phSize,0);
    for(int i = 0; i < phSize; i++)
        if(i<size()) { phVec[i] = values[i] + other[i];}
        else{ phVec[i] = other[i]; }
        
    return nVec(phVec);
}

nVec nVec::operator-(nVec& other){
    int phSize = lock_state_size(other);
    std::vector<float> phVec(phSize,0);
    for(int i = 0; i < phSize; i++)
        if(i<size()) { phVec[i] = values[i] - other[i];}
        else{ phVec[i] = other[i]; }
        
    return nVec(phVec);
}

nVec nVec::operator*(float scalar){
    std::vector<float> phVec = values;
    for(int i = 0; i < size(); i++)
        phVec[i] *= scalar;
    return nVec(phVec);
}

nVec nVec::operator/(float scalar){
    if(scalar == 0)
        return nVec(values);
    std::vector<float> phVec = values;
    for(int i = 0; i < size(); i++)
        phVec[i] /= scalar;
    return nVec(phVec);
}

float& nVec::operator[](int index){
    if(index < values.size())
        return values[index];
    return dummy;
}


int nVec::lock_state_size(nVec change){
    if(isLocked())
        return size();
    return std::max(size(), change.size());
}