#include <vector>
#include <cmath>
#include <bits/stdc++.h>

struct nVec {
    std::vector<float> values;

    nVec() : values(std::vector<float>(2,0)) {}
    nVec(int dimensions) : values(std::vector<float>(dimensions, 0)) {}
    nVec(std::vector<float> dimensions) : values(dimensions) {}
    
    float length() const;
    nVec normalized() const;
    int size() const;
    bool isLocked();
    void lock();
    void unlock();

    nVec operator+(nVec& other);
    nVec operator-(nVec& other);
    nVec operator*(float scalar);
    nVec operator/(float scalar);
    float& operator[](int index);

private:
    inline static float dummy = 0.0f;
    bool locked = false;
    int lock_state_size(nVec change);
};