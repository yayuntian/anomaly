# anomaly
基于时间序列的异常检测， 参照skyline、anomalyzer使用c++实现

参照： https://github.com/lytics/anomalyzer


# Example

```
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

    Anomalyzer anomaly = NewAnomalyzer(conf, data);
    //float ret = eval(anomaly);
    float ret = push(anomaly, 12);
    cout << "Anomalous Probability: " << ret << endl;

    return 0; 
}
```