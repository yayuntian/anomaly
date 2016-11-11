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

typedef struct {
    bool Delay;
    float Sensitivity;
    float UpperBound;
    float LowerBound;
    int ActiveSize;
    int referenceSize;
    int NSeasons;
    int PermCount;
    std::vector<std::string> Methods;
} AnomalyzerConf;


typedef struct {
    AnomalyzerConf Conf;
    std::vector<float> Data;
} Anomalyzer;


typedef float (*Algorithm)(std::vector<float>, AnomalyzerConf&);

typedef bool (*compare)(float, float);


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
static inline float sum(std::vector<float> x)
{
    return std::accumulate(x.begin(), x.end(), 0.0);
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


extern std::unordered_map<std::string, Algorithm> Algorithms;


float BootstrapKsTest(std::vector<float> data, AnomalyzerConf& conf);
float MagnitudeTest(std::vector<float> data, AnomalyzerConf& conf);
float DiffTest(std::vector<float> data, AnomalyzerConf& conf);
float RankTest(std::vector<float> data, AnomalyzerConf& conf);
float ReverseRankTest(std::vector<float> data, AnomalyzerConf& conf);
float CDFTest(std::vector<float> data, AnomalyzerConf& conf);
float FenceTest(std::vector<float> data, AnomalyzerConf& conf);
std::vector<float> split(std::vector<float> x, int start, int end);


float mean(const std::vector<float>& array);
float eval(Anomalyzer& anomaly);
