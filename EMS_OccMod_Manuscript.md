*Improving inference of metacommunity structure using multi-species occupancy models*
========================================================

#### Authors
Joseph R. Mihaljevic$^1$$^*$, Maxwell B. Joseph$^1$, Pieter T.J. Johnson$^1$

$^1$ *Department of Ecology and Evolutionary Biology, University of Colorado, Boulder, CO, 80309*

$^*$ *Corresponding Author*

Abstract
------------------
Blah blah blah

Introduction
-------------------

Broadly speaking, the metacommunity concept seeks to understand how spatial patterns of community composition emerge as a product of both local (e.g. competition) and regional (e.g. dispersal) dyanmics [@Leibold2004; @Chase2005; @Holyoak2005]. To date, metacommunity research can be separated into two nearly distinct categories: mechanism-based approaches and pattern-based approaches. The mechanism-based approach employs mechanistic modeling or controlled experiments to generate and test hypotheses related to, for example, the distribution of species across local sites [@Holyoak2005; @Cottenie2005; @Urban2009; @Logue2011; @Pillai2011; @Economo2011; @Carrara2012]. The pattern-based approach uses an inverse process, attempting to relate empirically observed patterns of species occurrences across a landscape to structuring mechanisms [@Leibold2002; @Presley2010; @Presley2010a; @Lopez-González2012; @Henriques-Silva2013; @Meynard2013]. While both approaches have led to substantial advances in our understanding of community structure and its underpinnings, methodological improvements can be developed to better integrate these often disparate lines of research. Improving metacommunity analytical tools could therefore lead to a more complete understanding of community patterns and sturcturing mechanisms.  

The pattern-based metacommunity approach, often referred to as the 'elements of metacommunity structure' (EMS) paradigm, relies on observed occurrences of species among 'patches' of habitat in the metacommunity (e.g. field sampling sites). This data is typically compiled into a species-by-site incidence matrix and statistical procedures are run to determine if the metacommunity exhibits any of twelve unique structures [@Leibold2002; @Presley2010]. Based on the specific structure observed, inferences can be made as to how the metacommunity assembles along dominant environmental axes. However, the efficacy of this pattern-to-process approach depends heavily upon the quality of the observed data. For instance, problems with species detection could lead to innaccurate incidence matrices and, therefore, innaccurate assessments of metacommunity structure. Fortunately, decades of advancements in occupancy modeling have led to powerful methods that deal with problems such as species detectability. Integrating occupancy modeling with the pattern-based EMS paradigm should improve our ability to assess metacommunity structure and the associated inferences about structuring processes.

Occupancy modeling relies on repeated sampling surveys to distinguish between the probability of a species occuring at a site and the probability of a species being detected at a site in which it occurs [**citations**]. These models can also use observed occurrences to simultaneously estimate the effects of covariates on detection and occurrence probabilities. Recently, dynamic multi-species occupancy models have been developed for use with longitudinal datasets to better estimate aspects of metacommunity structure, such as alpha, beta and gamma diversity metrics, as well as species-, site- and time-specific covariate effects on occurrence probabilities [@Dorazio2010; @Burton2012]. In their article, Dorazio et al. (2010) also specifically identify a closer union of their modeling framework with metacommunity theory as a priority in the field. While the integration of occupancy modeling with metacommunity theory is nascent, there is enormous potential to improve the study of metacommunity dynamics.

In this article, we will highlight some of the difficulties with current EMS methods and emphasize how integrating assessment of EMS with occupancy modeling will help to resolve these issues. We will also illustrate unique advantages of occupancy modeling that will complement the EMS paradigm. First, we will begin by briefly reviewing the methods involved in both the analysis of EMS and multi-species occupancy models. Next, we will verbally outline various problems or advantages and couple these explanations with illustrative figures generated with simluated datasets. Our motivation for this article is to highlight the breadth of utilities gained by integrating these two approaches, rather than to explore any specific question or topic in great detail. We have also made our modeling and simulation code publicly available and fully annotated so that these methods can be appropriately and broadly utilized. We believe our approach will lead to a more complete exploration of metacommunity structure and its underpinning mechanisms from occurrence datasets.

Elements of metacommunity structure
-------------------------------------

Multi-species occupancy models
-------------------------------------



