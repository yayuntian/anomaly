#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <climits>
#include <cmath>
#include <cfloat>


// Minimum representable floating-point number.
#define NA FLT_MIN


class AnomalyzerConf
{
public:
    AnomalyzerConf(float UpperBound, int ActiveSize, int NSeasons, std::vector<std::string> Methods);
    ~AnomalyzerConf();

    void validateConf();
    void printConf();

    bool Delay;
    float Sensitivity;
    float UpperBound;
    float LowerBound;
    int ActiveSize;
    int referenceSize;
    int NSeasons;
    int PermCount;
    std::vector<std::string> Methods;   
    
};



class Anomalyzer
{
public:

    Anomalyzer(const AnomalyzerConf& conf);
    Anomalyzer(const AnomalyzerConf& conf, const std::vector<float>& data);

    ~Anomalyzer();


    float eval();
    float push(float f);


    const AnomalyzerConf* Conf;
    std::vector<float> Data;

};


typedef float (*Algorithm)(const std::vector<float>&, const AnomalyzerConf&);
typedef bool (*compare)(float, float);


extern std::unordered_map<std::string, Algorithm> Algorithms;


float BootstrapKsTest(const std::vector<float>& data, AnomalyzerConf& conf);
float MagnitudeTest(const std::vector<float>& data, AnomalyzerConf& conf);
float DiffTest(const std::vector<float>& data, AnomalyzerConf& conf);
float RankTest(const std::vector<float>& data, AnomalyzerConf& conf);
float ReverseRankTest(const std::vector<float>& data, AnomalyzerConf& conf);
float CDFTest(const std::vector<float>& data, AnomalyzerConf& conf);
float FenceTest(const std::vector<float>& data, AnomalyzerConf& conf);


static inline bool greaterThan(float x, float y)
{
    if (x > y) {
        return true;
    } else {
        return false;
    }
}


static inline bool lessThan(float x, float y)
{
    if (x < y) {
        return true;
    } else {
        return false;
    }
}


// Sum returns the sum of the vector.
static inline float sum(const std::vector<float>& data)
{
    return std::accumulate(data.begin(), data.end(), 0.0);
}


    template <class T>
int getArrayLen(T& array)
{ 
    return (sizeof(array) / sizeof(array[0]));
}


    template <class T1, class T2>
bool exists(T1 v, T2& array)
{
    auto it = std::find(array.begin(), array.end(), v);
    if (it != array.end()) {
        return true;
    } else {
        return false;
    }
}


    template <class T>
void printVector(T& array)
{
    std::cout << "Vectors: [ ";
    for (auto it : array) {
        std::cout << it << " ";
    }
    std::cout <<"]"<< std::endl;
}


    template <class T>
void printMap(T& myMap)
{
    for(auto elem : myMap) {
        std::cout << elem.first << " " 
            //<< elem.second.first << " " 
            << elem.second << std::endl;
    }
}


// Returns a value within a given window (xmin and xmax).
static inline float cap(float x, float min, float max)
{
    return std::max(std::min(x, max), min);
}

// Mean returns the mean of the vector.
static inline float mean(const std::vector<float>& array)
{
    float s = sum(array);
    float n = (float) array.size();

    return s / n;
}
