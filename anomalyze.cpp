#include <iostream>
#include <cctype>
#include <cstdlib>
#include <list>
#include <iterator>
#include <queue>
#include <set>

#include <algorithms.h>

using namespace std;

/*
 * Use essentially similar weights.  However, if either the magnitude
 * or fence methods have high probabilities, upweight them significantly.
 */
float getWeight(string name, float prob)
{
    float weight = 0.5;

    std::vector<string> dynamicWeights = {"magnitude", "fence"};

    /* If either the magnitude and fence methods don't have any
       probability to contribute, we don't want to hear about it.
       If they do, we upweight them substantially.*/

    if (exists(name, dynamicWeights)) {
        if (prob > 0.8) {
            weight = 5.0;
        }
    }

    return weight;
}


// Variance caclulates the variance of the vector
float variance(const std::vector<float>& data)
{
    float m = mean(data);

    int n = data.size();
    if (n < 2) {
        return 0.0;
    }

    float ss = 0.0;

    for (int i = 0; i < n; i++) {
        ss += pow(data[i] - m, 2.0);
    }

    return ss / (n - 1);
}


// Sd calculates the standard deviation of the vector
float sd(const std::vector<float>& array)
{
    return sqrt(variance(array));
}


/*
 * weightedSum returns the weighted sum of the vector.  This is really only useful in
 * calculating the weighted mean.
 */
float weightedSum(const std::vector<float>& probs, const std::vector<float>& weights)
{
    if (probs.size() != weights.size()) {
        cout << "Length of weights unequal to vector length" << endl;
        return NA;
    }

    float ws = 0.0;
    int size = (int) probs.size();

    for (int i = 0; i < size; i++) {
        ws += probs[i] * weights[i];
    }

    return ws;
}


// WeightedMean returns the weighted mean of the vector for a given vector of weights.
float weightedMean(std::vector<float> probs, std::vector<float> weights)
{
    float ws = weightedSum(probs, weights);
    float sw = sum(weights);

    // if all the weights are zero return NA
    if (ws == NA || !isgreater(sw, 0.0)) {
        return NA;
    }

    return ws / sw;
}





AnomalyzerConf:: AnomalyzerConf(float UpperBound, int ActiveSize,
    int NSeasons, std::vector<std::string> Methods)
{
        UpperBound = UpperBound;
        ActiveSize = ActiveSize;
        NSeasons = NSeasons;
        Methods = Methods;

        Delay = true;
        PermCount = 500;
        Sensitivity = 0.1;
        LowerBound = NA;

        validateConf();
}


AnomalyzerConf::~AnomalyzerConf()
{
    
}

void AnomalyzerConf::printConf()
{
    cout << "Delay:" << Delay << endl;
    cout << "Sensitivity:" << Sensitivity << endl;
    cout << "UpperBound:" << UpperBound << endl;
    cout << "LowerBound:" << LowerBound << endl;
    cout << "ActiveSize:" << ActiveSize << endl;
    cout << "referenceSize:" << referenceSize << endl;
    cout << "NSeasons:" << NSeasons << endl;
    cout << "PermCount:" << PermCount << endl;

    printVector(Methods);
}


void AnomalyzerConf::validateConf()
{
    // if supplied, make sure the detection methods are supported
    std::vector<string> supportedMethods = {"magnitude", "diff", "highrank", "lowrank", "fence", "ks", "cdf"};
    std::vector<string> minimumMethods = {"magnitude", "ks"};

    if (Methods.empty()) {
        Methods = minimumMethods;
    } else {
        for (auto it : Methods) {
            if (!exists(it, supportedMethods)) {
                std::cout << "Unsupported detection method: " << it << endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    // if number of seasons are not specified, default it to 4
    if (NSeasons == 0) {
        NSeasons = 4;
    }
    // if delay is not specified, default to false. this means calculations
    // of anomalousness will be returned as soon as we can

    // make reference window some multiple of the active window size
    referenceSize = NSeasons * ActiveSize;

    // window sizes must be positive ints
    if (ActiveSize < 1) {
        printf("Active window size must be at least of size 1\n");
        exit(EXIT_FAILURE);
    }

    if (referenceSize < 4) {
        printf("The combination of active window (%d) and nseasons (%d) yields"
                " a reference window that is too small for analysis."
                "  Please increase one or both.\n", ActiveSize, NSeasons);
        exit(EXIT_FAILURE);
    }


    // validation for the fence test
    if (exists("fence", Methods)) {
        if (UpperBound == LowerBound) {
            cout << "Fence test included with identical bounds on the fences" << endl;
            exit(EXIT_FAILURE);
        } else if (UpperBound < LowerBound) {
            printf("UpperBound (%0.2f) was lower than the LowerBound (%0.2f)\n", UpperBound, LowerBound);
            exit(EXIT_FAILURE);
        }
    }

    // validation for the permutation tests
    if (exists("highrank", Methods) || exists("lowrank", Methods) ||
            exists("ks", Methods) || exists("diff", Methods)) {
        if (PermCount == 0) {
            PermCount = 500;
        }
    }

    if (exists("magnitude", Methods)) {
        if (Sensitivity == 0.0) {
            Sensitivity = 0.1;
        }
    }

    this->printConf();
}


Anomalyzer::Anomalyzer(const AnomalyzerConf& conf)
{

    Conf = &conf;
    Data = std::vector<float>(1000);
}


Anomalyzer::Anomalyzer(const AnomalyzerConf& conf, const std::vector<float>& data)
{

    Conf = &conf;
    Data = data;
}


Anomalyzer::~Anomalyzer()
{

}


/*
 * Return the weighted average of all statistical tests
 * for anomaly detection, which yields the probability that
 * the currently observed behavior is anomalous.
 */
float Anomalyzer::eval()
{
    int threshold = Conf->referenceSize + Conf->ActiveSize;

    if (Conf->Delay && (int) Data.size() < threshold) {
        return 0.0;
    }

    std::unordered_map<std::string, float> probmap;
    for (auto method : Conf->Methods) {
        auto algorithm = Algorithms[method];

        float prob = cap(algorithm(Data, *Conf), 0, 1);
        if (prob != NA) {
            // if highrank and lowrank methods exist then only listen to the max of either
            if (method == "highrank" || method == "lowrank") {
                auto it = probmap.find("rank");
                if (it == probmap.end()) {
                    probmap["rank"] = 0;
                }
                probmap["rank"] = std::max(probmap["rank"], prob);
            } else {
                probmap[method] = prob;
            }
        }
    }

    std::vector<float> probs, weights;
    for (auto it : probmap) {
        string method = it.first;
        float prob = it.second;

        //cout << "Method: " << method << ", Prob: " << prob << endl;
        if (method == "magnitude" && prob < Conf->Sensitivity) {
            return 0.0;
        }

        probs.push_back(prob);
        weights.push_back(getWeight(method, prob));
    }

    /* ignore the error since we force the length 
       of probs and the weights to be equal */
    float weighted = weightedMean(probs, weights);

    return (weighted == NA) ? 0.0 : weighted;
}


float Anomalyzer::push(float f)
{
    // add the new point to the data
    Data.push_back(f);

    // evaluate the anomalous probability
    return eval();
}