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
Anomalyzer NewAnomalyzer(AnomalyzerConf& conf, const std::vector<float>& data);
float eval(Anomalyzer& anomaly);
std::vector<float> Rank(const std::vector<float>& data);
float push(Anomalyzer& anomaly, float data);


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
    AnomalyzerConf conf;

    conf.Delay = true;
    conf.Sensitivity = 0.1;
    conf.UpperBound = 5;
    conf.LowerBound = NA;
    conf.ActiveSize = 2;
    conf.NSeasons = 4;
    conf.PermCount = 500;
    // "magnitude", "diff", "highrank", "lowrank", "fence", "ks", "cdf"
    conf.Methods = {"magnitude", "diff", "highrank", "lowrank", "fence", "ks", "cdf"};

    std::vector<float> data = {0.01, 1.3, 0.5, 0.2, 10.51, 0.31, 0.1, 0.4, 0.11, 5.2, 0.12, 0.34};

#if 1
    Anomalyzer anomaly = NewAnomalyzer(conf, data);
    //float ret = eval(anomaly);
    float ret = push(anomaly, 12);

    cout << "Anomalous Probability: " << ret << endl;
#endif
#if 0
    split_test(x);
    mean_test(x);
    variance_test(x);

    float ret = weightedMean(x, y);
    cout << "weightedSum test: "<< ret << endl;

    auto m = extractReference(x, 4, 2, 4);
    auto n = extractActive(x, 4, 2, 4);

    printVector(m);
    printVector(n);

    auto v = Rank(data);
    printVector(v);
#endif
    return 0; 
}
