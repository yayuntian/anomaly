#include <iostream>
#include <cctype>
#include <cstdlib>
#include <list>
#include <iterator>
#include <queue>
#include <set>

#include <algorithms.h>


void printConf(AnomalyzerConf conf) {
    cout << "Delay:" << conf.Delay << endl;
    cout << "Sensitivity:" << conf.Sensitivity << endl;
    cout << "UpperBound:" << conf.UpperBound << endl;
    cout << "LowerBound:" << conf.LowerBound << endl;
    cout << "ActiveSize:" << conf.ActiveSize << endl;
    cout << "referenceSize:" << conf.referenceSize << endl;
    cout << "NSeasons:" << conf.NSeasons << endl;
    cout << "PermCount:" << conf.PermCount << endl;

    printVector(conf.Methods);
}


void validateConf(AnomalyzerConf& conf) {
    // if supplied, make sure the detection methods are supported
    std::vector<string> supportedMethods = {"magnitude", "diff", "highrank", "lowrank", "fence", "ks", "cdf"};
    std::vector<string> minimumMethods = {"magnitude", "ks"};

    if (conf.Methods.empty()) {
        conf.Methods = minimumMethods;
    } else {
        for (auto it : conf.Methods) {
            if (!exists(it, supportedMethods)) {
                std::cout << "Unsupported detection method: " << it << endl;
                exit(1);
            }
        }
    }

    // if number of seasons are not specified, default it to 4
    if (conf.NSeasons == 0) {
        conf.NSeasons = 4;
    }
    // if delay is not specified, default to false. this means calculations
    // of anomalousness will be returned as soon as we can

    // make reference window some multiple of the active window size
    conf.referenceSize = conf.NSeasons * conf.ActiveSize;

    // window sizes must be positive ints
    if (conf.ActiveSize < 1) {
        printf("Active window size must be at least of size 1\n");
        exit(1);
    }

    if (conf.referenceSize < 4) {
        printf("The combination of active window (%d) and nseasons (%d) yields"
                " a reference window that is too small for analysis."
                "  Please increase one or both.\n", conf.ActiveSize, conf.NSeasons);
        exit(1);
    }


    // validation for the fence test
	if (exists("fence", conf.Methods)) {
        if (conf.UpperBound == conf.LowerBound) {
            cout << "Fence test included with identical bounds on the fences" << endl;
            exit(1);
        } else {
            if (conf.UpperBound < conf.LowerBound) {
                printf("UpperBound (%0.2f) was lower than the LowerBound (%0.2f)\n", conf.UpperBound, conf.LowerBound);
                exit(1);
            }
        }
	}

    // validation for the permutation tests
    if (exists("highrank", conf.Methods) || exists("lowrank", conf.Methods) ||
		exists("ks", conf.Methods) || exists("diff", conf.Methods)) {
        if (conf.PermCount == 0) {
            conf.PermCount = 500;
        }
    }

    if (exists("magnitude", conf.Methods)) {
        if (conf.Sensitivity == 0.0) {
            conf.Sensitivity = 0.1;
        }
    }

    printConf(conf);
}



Anomalyzer NewAnomalyzer(AnomalyzerConf& conf, std::vector<float> data)
{

    validateConf(conf);

    Anomalyzer anomaly;

    anomaly.Conf = conf;
    anomaly.Data = data;

    return anomaly;
}



void update(Anomalyzer& anomaly, std::vector<float> data)
{
    // add new elememnts to the vector
    anomaly.Data.insert(anomaly.Data.end(), data.begin(), data.end());

    printVector(anomaly.Data);

    int offset = anomaly.Data.size() - (anomaly.Conf.ActiveSize + anomaly.Conf.referenceSize);
    if (offset < 0) {
        offset = 0;
    }

    anomaly.Data.erase(anomaly.Data.begin(), anomaly.Data.begin() + offset);

    printVector(anomaly.Data);
}


float push(Anomalyzer& anomaly, float data)
{
    // add the new point to the data
    anomaly.Data.push_back(data);

    // evaluate the anomalous probability
    return eval(anomaly);
}


// Use essentially similar weights.  However, if either the magnitude
// or fence methods have high probabilities, upweight them significantly.
float getWeight(Anomalyzer& anomaly, string name, float prob)
{
    float weight = 0.5;

    std::vector<string> dynamicWeights = {"magnitude", "fence"};

    // If either the magnitude and fence methods don't have any
    // probability to contribute, we don't want to hear about it.
    // If they do, we upweight them substantially.

    if (exists(name, dynamicWeights)) {
        if (prob > 0.8) {
            weight = 5.0;
        }
    }

    return weight;
}


// Mean returns the mean of the vector.
float mean(const std::vector<float> array)
{
    float s = sum(array);
    float n = (float) array.size();

    return s / n;
}


// Variance caclulates the variance of the vector
float variance(std::vector<float>& x)
{
    float m = mean(x);

    int n = x.size();
    if (n < 2) return 0.0;

    float ss = 0.0;
    for (float v : x) {
        ss += pow(v - m, 2.0);
    }

    return ss / (n - 1);
}


// Sd calculates the standard deviation of the vector
float sd(std::vector<float>& array)
{
    return sqrt(variance(array));
}


// weightedSum returns the weighted sum of the vector.  This is really only useful in
// calculating the weighted mean.
float weightedSum(std::vector<float>& probs, std::vector<float>& weights)
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


// Return the weighted average of all statistical tests
// for anomaly detection, which yields the probability that
// the currently observed behavior is anomalous.
float eval(Anomalyzer& anomaly)
{
    int threshold = anomaly.Conf.referenceSize + anomaly.Conf.ActiveSize;

    if (anomaly.Conf.Delay && (int) anomaly.Data.size() < threshold) {
        return 0.0;
    }

    std::unordered_map<std::string, float> probmap;
    for (auto method : anomaly.Conf.Methods) {
        auto algorithm = Algorithms[method];

        float prob = cap(algorithm(anomaly.Data, anomaly.Conf), 0, 1);
        if (prob != NA) {
            // if highrank and lowrank methods exist then only listen to
            // the max of either
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

		cout << "Method: " << method << ", Prob: " << prob << endl;
        if (method == "magnitude" && prob < anomaly.Conf.Sensitivity) {
            return 0.0;
        }

        probs.push_back(prob);
        weights.push_back(getWeight(anomaly, method, prob));
    }

    // ignore the error since we force the length of probs
    // and the weights to be equal
    float weighted = weightedMean(probs, weights);

	return (weighted == NA) ? 0.0 : weighted;
}


// Get the results and weights of each test. Useful for debugging
void evalByTest(Anomalyzer& anomaly)
{
    std::unordered_map<std::string, float> probmap;
    for (auto method : anomaly.Conf.Methods) {

        auto algorithm = Algorithms[method];
        float prob = cap(algorithm(anomaly.Data, anomaly.Conf), 0, 1);

        if (prob != NA) {
            // if highrank and lowrank methods exist then only listen to
            // the max of either
            if (method == "highrank" || method == "lowrank") {
                if (isnan(probmap["rank"])) {
                    probmap["rank"] = 0;
                }
                probmap["rank"] = std::max(probmap["rank"], prob);
            } else {
                probmap[method] = prob;
            }
        }
    }

    std::unordered_map<std::string, float> weightmap;
    for (auto it : probmap) {
        string method = it.first;
        float prob = it.second;

        cout << it.first << it.second << endl;
        weightmap[method] = getWeight(anomaly, method, prob);
    }

    printMap(probmap);
    printMap(weightmap);
}
