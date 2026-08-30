### Multiple Comparison Synergy P-Values Adjustment

SynergyLMM provides time-resolved estimates of drug combination effects across multiple time points, which raises the question whether p-values should be adjusted for multiple testing. Users can choose to apply common correction methods, such as Bonferroni or false discovery
rate (FDR), to control for type I error. However, we do not recommend systematic adjustment of p-values.

Instead, we recommend reporting the nominal p-values with confidence intervals. If stricter control of type I error is desired, users can lower the significance threshold (e.g., $\alpha = 0.01$). Based on the study design or users’ interpretation of single or multiple p-values, users can decide if an appropriate p-value adjustment is needed.

**When to adjust?**

- Adjust if the goal is to reduce false positives before advancing to other confirmatory studies or clinical trials.
- Avoid stringent adjustment in exploratory studies, where detection sensitivity can be more important than specificity.

Available options for adjustment include those in [stats::p.adjust](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust):
- `none`: Default. No multiple comparison p-value adjustment.
- `BH`: [Benjamini & Hochberg (1995)](https://academic.oup.com/jrsssb/article/57/1/289/7035855).
- `fdr`: False discovery rate. Alias and equivalent for 'BH'.
- `holm`: [Holm (1979)](https://www.jstor.org/stable/4615733).
- `hochberg`: [Hochberg (1988)](https://academic.oup.com/biomet/article/75/4/800/423177).
- `hommel`: [Hommel (1988)](https://academic.oup.com/biomet/article/75/2/383/292949).
- `bonferroni`: Bonferroni correction (p-values are multiplied by the number of comparisons).
- `BY`: [Benjamini & Yekutieli (2001)](https://projecteuclid.org/journals/annals-of-statistics/volume-29/issue-4/The-control-of-the-false-discovery-rate-in-multiple-testing/10.1214/aos/1013699998.full).