#include <iostream>
#include <cctype>
#include <cstdlib>
#include <list>
#include <iterator>
#include <queue>
#include <set>

#include <algorithms.h>

using namespace std;

float variance(const std::vector<float>& x);
float sd(std::vector<float>& array);
float weightedSum(std::vector<float>& probs, std::vector<float>& weights);
float weightedMean(std::vector<float> probs, std::vector<float> weights);
std::vector<float> extractReference(const std::vector<float>& data,
        int refSize, int activeSize, int minRefSize);
std::vector<float> extractActive(const std::vector<float>& data,
        int refSize, int activeSize, int minRefSize);
float weightExp(float x, float base);
std::vector<float> Rank(const std::vector<float>& data);

std::vector<float> split(const std::vector<float>& data, int start, int end);


void variance_test(std::vector<float>& x)
{
    cout << variance(x) << endl;
}

void mean_test(std::vector<float>& x)
{
    auto ret = mean(x);
    cout << ret << endl;
}

void split_test(const std::vector<float>& x)
{
    auto y = split(x, 3, 5);

    y[0] = 10;
    printVector(x);
    printVector(y);
}

int main()
{
    std::vector<float> data = {0.01, 1.3, 0.5, 0.2, 10.51, 0.31, 0.1, 0.4, 0.11, 5.2, 0.12, 0.34};

    // "magnitude", "diff", "highrank", "lowrank", "fence", "ks", "cdf"
    AnomalyzerConf conf = AnomalyzerConf(5, 2, 4, {"magnitude", "diff", "highrank", "lowrank", "fence", "ks", "cdf"});
    Anomalyzer anomaly = Anomalyzer(conf, data);
    float ret = anomaly.eval();

    cout << "Anomalous Probability: " << ret << endl;

    return 0; 
}