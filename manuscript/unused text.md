Unused text:
========================================================

**********************************************
Additional benefits of a merger (**To include or not?**)
-----------------------------------

*Covariate effects*  
Another issue with the EMS strategy that can be emeliorated by a merger with occupancy modeling is the way in which covariate effects are explored. As discussed above, many researchers extract the first ordination axis score for each sampled site and associate these scores to any number of site-specific covariates using multiple univariate correlation tests. This approach is flawed due to inflation of the overall type I error rate caused by conducting multiple tests. Additionally, regardless of the type of statistical test, using the ordination scores from a single axis could be leaving out valuable information, and these ordination scores are extracted from a matrix that is most likely flawed by detection error in the first place. Occupancy modeling helps to remedy these issues by estimating species-specific covariate effects on occurrence probability based on the raw data, not ordination scores. Then, using Bayesian model selection and model averaging, dominant covariates and their effects can be estimated without the need for multiple tests. 

The influence of confounding covariate effects on detection and occurrence probabilities can also be explored within the occupancy modeling framework. For example, if a single covariate, such as vegetation density, has positive effects on species' occupancy but negative effects on species' detectability, how might this influence our ability to detect or accurately assign metacommunity structure? The simulation method described above could be used to explore this question in a future study by altering the correlation between the covariate's effect on each species' detection and occurence probabilities and observing how metacommunity structure changes.  


*Temporal dynamics*    

With longitudinal metacommunity data becoming more and more available, an interesting question that arises is how does metacommunity structure vary across time points and what are the driving mechanisms for these changes. For example, @Eros2014 explored how metacommunity structure of stream fishes changed over time and how this structure was influenced by environmental covariates. Dynamic occupancy modeling has a few complementary features that could aid in a more cohesive approach to studying such questions.

Occupancy modeling could be used to distinguish if changes in metacommunity structure are more related to changing covariate values at sites over time, whether the dominant structuring covariate changes identity over time (e.g. from altitude to pH), or if metacommunity structure changes due to idiosyncratic changes in how species respond to the same covariates over time. Additionally, @Dorazio2010 showed how multi-species occupancy modeling can be used to estimate species- and time-specific covariates on detection, occurrence, colonization and persistence probabilities. Understanding colonization and persistence (i.e. extinction) dynamics could aid in inferring more specific mechanisms for metacommunity structural changes.  


