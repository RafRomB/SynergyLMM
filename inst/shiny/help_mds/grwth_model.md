## Growth Models

Currently, two growth kinetics models are available: exponential growth, and Gompertz growth.

**Exponential**

SynergyLMM will fit a linear mixed-effect model (LMM) assuming that the tumor growth follows an exponential kinetics. Any departure from this assumption can be tested based on the results in **Model Diagnostics** and **Model Performance** tabs.

The model formula for the LMM following the **exponential tumor growth** is: $$\log RTV_{i}(t) = \beta_{T_i} \cdot t + b_i \cdot t + \varepsilon_{i} (t).$$

-   $\log RTV_{i}(t)$ denotes the value of the logarithm of the relative tumor volume measured for subject $i$ at time $t$.

-   $\beta_{T_i}$ represents the fixed effects for each treatment $T_i$, where $T_i \in \{Control, DrugA, DrugB, Combination\}$ in the case of two-drugs combination experiments, or $T_i \in \{Control, DrugA, DrugB, DrugC, Combination\}$ in the case of three-drugs combination experiments, and indicates the tumor-specific growth rate for each treatment group.

-   $b_i \cdot t$ corresponds to the subject-specific random slope that takes into account the longitudinal nature of the data, where $b_i \sim \mathcal{N}(0,\sigma^2_b)$ is the random effect for subject $i$.

-   $\varepsilon_{i}(t) \sim \mathcal{N}(0,\sigma^2)$ is the residual error term.

**Gompertz**

SynergyLMM will fit a non-linear mixed effect (NLME) model assuming the tumor growth follows a Gompertz growth kinetics. Any departure from this assumption can be tested based on the results in **Model Diagnostics** and **Model Performance** tabs.

The model formula for the non-linear mixed-effect model following the **Gompertz tumor growth** is: $$\log RTV_{i}(t) = \frac{r_{0_{T_i}}+b_{0_i}}{\rho_{T_i}+b_{1_i}}(1-e^{-(\rho_{T_i}+b_{1_i})\cdot t})+\varepsilon_{i}(t).$$

-   $\log RTV_{i}(t)$ denotes the value of the logarithm of the relative tumor volume measured for subject $i$ at time $t$.
-   $T_i \in \{Control, DrugA, DrugB, Combination\}$ in the case of two-drugs combination experiments, or $T_i \in \{Control, DrugA, DrugB, DrugC, Combination\}$ in the case of three-drugs combination experiments, indicates the treatment groupt of individual $i$.
-   $r_{0_{T_i}}$ is the fixed effect for the initial growth rate for treatment group $T_i$.
-   $\rho_{T_i}$ is the fixed effect for the constant accounting for the reduction in the tumor growth rate for treatment group $T_i$.
-   $b_{0_i} \sim \mathcal{N}(0,\sigma^2_{{r_0}})$ is the random effect on $r_0$ for individual $i$.
-   $b_{1_i} \sim \mathcal{N}(0,\sigma^2_{{\rho}})$ is the random effect on $\rho$ for individual $i$.
-   $\varepsilon_{i}(t) \sim \mathcal{N}(0,\sigma^2)$ is the residual error term.
