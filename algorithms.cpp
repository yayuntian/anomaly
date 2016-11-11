#include <iostream>
#include <ctime>
#include <cstdlib>

#include <algorithms.h>

using namespace std;

std::unordered_map<std::string, Algorithm> Algorithms = {
    {"magnitude", MagnitudeTest},
    {"diff", DiffTest},
    {"highrank", RankTest},	
    {"lowrank", ReverseRankTest},
    {"cdf", CDFTest},
    {"fence", FenceTest},
    {"ks", BootstrapKsTest}
};


std::vector<float> split(std::vector<float> data, int start, int end)
{
    int n = end - start;
    if (n <= 0) {
        return std::vector<float>();
    }

    std::vector<float> y(n);
    for (int i = 0; i < n; i++) {
        y[i] = data[start++];
    }
    return y;
}


/*
 * Return a vector slice for the active window and reference window.
 * Some tests require different minimum thresholds for sizes of reference windows.
 * This can be specified in the minRefSize parameter. If size isn't important, use -1
 */
std::vector<float> extractReference(const std::vector<float>& data,
        int refSize, int activeSize, int minRefSize)
{
    int n = data.size();

    activeSize = std::min(activeSize, n);
    refSize = std::min(refSize, n - activeSize);

    /* make sure the reference size is at least as big as the active size
       note that this penalty might be overly severe for some tests */
    if (refSize < minRefSize) {
        std::cout << "Reference size must be at least as big as active size" << std::endl;
        return std::vector<float>();
    }

    // return reference windows
    return split(data, n - activeSize - refSize, n - activeSize);

}

// same as extractReference
std::vector<float> extractActive(const std::vector<float>& data,
        int refSize, int activeSize, int minRefSize)
{
    int n = data.size();

    activeSize = std::min(activeSize, n);
    refSize = std::min(refSize, n - activeSize);

    /* make sure the reference size is at least as big as the active size
       note that this penalty might be overly severe for some tests */
    if (refSize < minRefSize) {
        std::cout << "Reference size must be at least as big as active size" << std::endl;
        return std::vector<float>();
    }

    // return active windows
    return split(data, n-activeSize, n);
}


/*
 * This is a function will sharply scale values between 0 and 1 such that
 * smaller values are weighted more towards 0. A larger base value means a
 * more horshoe type function.
 */
float weightExp(float num, float base)
{
    return (pow(base, num) - 1) / (base - 1);
}


/*
 * This function can be used to test whether or not data is getting close to a
 * specified upper or lower bound.
 */
float FenceTest(std::vector<float> data, AnomalyzerConf& conf)
{
    // we don't really care about a reference window for this one
    auto active = extractActive(data, conf.referenceSize, conf.ActiveSize, -1);

    float m = mean(active);
    float distance = 0.0;
    if (conf.LowerBound == NA) {
        // we only care about distance from the upper bound
        distance = m / conf.UpperBound;
    } else {
        // we care about both bounds, so measure distance from midpoint
        float bound = (conf.UpperBound - conf.LowerBound) / 2;
        float mid = conf.LowerBound + bound;

        distance = (abs(m - mid)) / bound;
    }

    return weightExp(cap(distance, 0, 1), 10);
}


// Diff returns a vector of length (n - 1) of the differences in the input vector
std::vector<float> Diff(std::vector<float>& data)
{
    int n = data.size();
    if (n < 2) {
        return std::vector<float>();
    }

    std::vector<float> d(n - 1);
    for (int i = 1; i < n; i++) {
        d[i - 1] = abs(data[i] - data[i - 1]);
    }

    return d;
}


// RelDiff returns a vector of the relative differences of the input vector
std::vector<float> RelDiff(const std::vector<float>& data)
{
    int n = data.size();
    if (n < 2) {
        return std::vector<float>();
    }

    std::vector<float> d(n - 1);
    for (int i = 1; i < n; i++) {
        d[i - 1] = std::abs((data[i] - data[i - 1])/ data[i]);
    }

    return d;
}


// Rank returns a vector of the ranked values of the input vector.
std::vector<float> Rank(const std::vector<float>& data)
{
    std::vector<float> y = data;
    std::sort(y.begin(), y.end());

    // equivalent to a minimum rank (tie) method
    int rank = 0;
    int size = data.size();
    std::vector<float> ranks(size, -1);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (y[i] == data[j] && ranks[j] == -1) {
                ranks[j] = (float) rank;
            }

        }
        rank++;
    }

    return ranks;
}


// Order returns a vector of untied ranks of the input vector.
std::vector<float> Order(std::vector<float> data)
{
    std::vector<float> y = data;
    std::sort(y.begin(), y.end());

    int rank = 0;
    int size = y.size();
    std::vector<float> order(size, -1);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (y[i] == data[j] && order[j] == -1) {
                order[j] = (float) rank;
                rank++;
                break;
            }
        }
    }

    return order;
}


/*
 * Generates permutations of reference and active window values to determine
 * whether or not data is anomalous. The number of permutations desired has
 * been set to 500 but can be increased for more precision.
 */
float DiffTest(std::vector<float> data, AnomalyzerConf& conf)
{
    // Find the differences between neighboring elements and rank those differences.
    auto rd = RelDiff(data);
    auto ranks = Rank(rd);

    /*
     * The indexing runs to length-1 because after applying .Diff(), We have
     * decreased the length of out vector by 1.
     */
    auto active = extractActive(ranks, conf.referenceSize - 1, conf.ActiveSize, conf.ActiveSize);
    if (active.empty()) {
        return NA;
    }
    /*
     * Consider the sum of the ranks across the active data. This is the sum that
     * we will compare our permutations to.
     */
    float activeSum = sum(active);
    int significant = 0;

    /*
     * Permute the active and reference data and compute the sums across the tail
     * (from the length of the reference data to the full length).
     */
    for (int i = 0; i < conf.PermCount; i++) {
        std::random_shuffle(data.begin(), data.end());
        auto rd = RelDiff(data);

        auto permRanks = Rank(rd);
        auto permActive = extractActive(permRanks, conf.referenceSize - 1, conf.ActiveSize, conf.ActiveSize);

        /*
         * If we find a sum that is less than the initial sum across the active data,
         * this implies our initial sum might be uncharacteristically high. We increment
         * our count.
         */
        if (sum(permActive) < activeSum) {
            significant++;
        }
    }
    // We return the percentage of the number of iterations where we found our initial sum to be high.
    return (float) significant / (float) conf.PermCount;
}



/*
 * Generates permutations of reference and active window values to determine
 * whether or not data is anomalous. The number of permutations desired defaults
 * to 500 but can be increased for more precision. A comparison function above
 * can be specified to create Rank and ReverseRank tests.
 */
float RankTest(std::vector<float> data, AnomalyzerConf& conf, compare comparison)
{
    // Rank the elements of a vector
    auto ranks = Rank(data);

    auto active = extractActive(ranks, conf.referenceSize, conf.ActiveSize, conf.ActiveSize);
    if (active.empty()) {
        return NA;
    }

    /* Consider the sum of the ranks across the active data. This is the sum that
       we will compare our permutations to. */
    float activeSum = sum(active);
    int significant = 0;

    /* Permute the active and reference data and compute the sums across the tail
       (from the length of the reference data to the full length).*/
    for (int i = 0; i < conf.PermCount; i++) {
        std::random_shuffle(data.begin(), data.end());
        auto permRanks = Rank(data);
        auto permActive = extractActive(permRanks, conf.referenceSize, conf.ActiveSize, conf.ActiveSize);

        /* If we find a sum that is less than the initial sum across the active data,
           this implies our initial sum might be uncharacteristically high. We increment our count.*/
        float permSum = sum(permActive);
        if (comparison(permSum, activeSum)) {
            significant++;
        }
    }
    // We return the percentage of the number of iterations where we found our initial sum to be high.
    return (float) significant / (float) conf.PermCount;
}

// high rank
float RankTest(std::vector<float> data, AnomalyzerConf& conf)
{
    return RankTest(data, conf, lessThan);
}


// low rank
float ReverseRankTest(std::vector<float> data, AnomalyzerConf& conf)
{
    return RankTest(data, conf, greaterThan);
}


/*
 * Ecdf returns the empirical cumulative distribution function.  The ECDF function
 * will return the percentile of a given value relative to the vector.
 */
float Ecdf(std::vector<float> data, float q)
{
    std::vector<float> y = data;
    std::sort(y.begin(), y.end());

    int n = data.size();

    for (int i = 0; i < n; i++) {
        if (q < y[i]) {
            return i / (float) n;
        }
    }
    return 1.0;
}


// Generates the cumulative distribution function using the difference in the means for the data.
float CDFTest(std::vector<float> data, AnomalyzerConf& conf)
{
    auto diffs = Diff(data);

    auto reference = extractReference(diffs, conf.referenceSize - 1, conf.ActiveSize, conf.ActiveSize);
    auto active = extractActive(diffs, conf.referenceSize - 1, conf.ActiveSize, conf.ActiveSize);
    if (reference.empty() || active.empty()) {
        return NA;
    }

    // Difference between the active and reference means.
    float activeDiff = mean(active) - mean(reference);

    // Apply the empirical distribution function to that difference.
    float percentile = Ecdf(reference, activeDiff);

    // Scale so max probability is in tails and prob at 0.5 is 0.
    return (2 * std::abs(0.5 - percentile));
}



/*
 * Generates the percent difference between the means of the reference and active
 * data. Returns a value scaled such that it lies between 0 and 1.
 */
float MagnitudeTest(std::vector<float> data, AnomalyzerConf& conf)
{
    auto reference = extractReference(data, conf.referenceSize, conf.ActiveSize, 1);
    auto active = extractActive(data, conf.referenceSize, conf.ActiveSize, 1);
    if (reference.empty() || active.empty()) {
        return NA;
    }

    float activeMean = mean(active);
    float refMean = mean(reference);

    // If the baseline is 0, then the magnitude should be Inf, but we'll round to 1.
    if (refMean == 0) {
        return (activeMean == 0) ? 0 : 1;
    }

    float pdiff = abs(activeMean - refMean) / refMean;
    return pdiff;
}


// A helper function for KS that rescales a vector to the desired length npoints.
std::vector<float> interpolate(float min, float max, int npoints)
{
    std::vector<float> interp(npoints);
    float step = (max - min) / (float) (npoints - 1);
    interp[0] = min;

    for (int i = 1; i < npoints; i++) {
        interp[i] = interp[i - 1] + step;
    }

    return interp;
}


// Calculate a Kolmogorov-Smirnov test statistic.
float KsStat(std::vector<float> data, AnomalyzerConf conf)
{
    auto reference = extractReference(data, conf.referenceSize, conf.ActiveSize, 1);
    auto active = extractActive(data, conf.referenceSize, conf.ActiveSize, 1);
    if (reference.empty() || active.empty()) {
        return NA;
    }

    int n1 = reference.size();
    int n2 = active.size();

    if (n1 % n2 != 0) {
        return NA;
    }

    /* We want the reference and active vectors to have the same length n, so we
       consider the min and max for each and interpolated the points between.*/

    float min = std::min(*std::min_element(reference.begin(), reference.end()),
            *std::min_element(active.begin(), active.end()));

    float max = std::max(*std::max_element(reference.begin(), reference.end()),
            *std::max_element(active.begin(), active.end()));


    std::vector<float> interpolated = interpolate(min, max, n1+n2);

    // Then we apply the distribution function over the interpolated data.
    int size = interpolated.size();
    std::vector<float> activeDist(size);
    std::vector<float> refDist(size);

    for (int i = 0; i < size; i++) {
        activeDist.push_back(Ecdf(active, interpolated[i]));
        refDist.push_back(Ecdf(reference, interpolated[i]));
    }

    // Find the maximum displacement between both distributions.
    float d = 0.0;
    for (int i = 0; i < n1+n2; i++) {
        d = std::max(d, abs(activeDist[i] - refDist[i]));
    }

    return d;
}

/*
 * Calculates the Kolmogorov-Smirnov test over active and reference windows
 * and compares that value to KS test scores obtained
 * after permuting all elements in the set.
 */
float BootstrapKsTest(std::vector<float> data, AnomalyzerConf& conf)
{
    float dist = KsStat(data, conf);
    if (dist == NA) {
        return NA;
    }

    int significant = 0;
    for (int i = 0; i < conf.PermCount; i++) {
        std::random_shuffle(data.begin(), data.end());
        float permDist = KsStat(data, conf);

        if (permDist < dist) {
            significant++;
        }
    }
    return significant / (float) conf.PermCount;
}
