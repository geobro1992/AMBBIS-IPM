# AMBBIS-IPM
Supporting code for the manuscript: Chandler, H. C., G. C. Brooks, Y. Jiao, and C. A. Haas. Predicting long-term population viability for an imperiled salamander under future climate changes.

Population viability analyses (PVAs) represent a key component of many recovery plans for threatened and endangered species. Demography links the processes that affect individuals to population-level patterns, and hence projections constructed from demographic data are the most common tools for PVAs. We constructed a size-structured integral projection model (IPM) for the United States federally endangered Reticulated Flatwoods Salamander, _Ambystoma bishopi_, to evaluate demographic influences on population growth and predict the efficacy of future management actions. Climate change has been shown to have a variety of effects on amphibians breeding in ephemeral wetlands, altering both demographic processes and underlying environmental conditions. These changes can threaten the long-term viability of populations and limit the effectiveness of ongoing management programs. To examine the response of flatwoods salamander populations to climate change, we incorporated climate projections into the stochastic Integral Projection Model. We examined 6 climate change scenarios, combining 3 General Circulation Models (GCM) under two emission scenarios to estimate extinction probability from 2030–2100. 


## contents
### scripts
**RFS_IPM_kernels_and_parameters.R** - R code to create the IPM kernels\
**RFS_IPM_Projections.R** - R code to simulate population dynamics through time under various climate scenarios\
**RFS_IPM_Figures.R** - R code to generate prediction and sensitivity plots

### data
**P4_climate_projections.csv** - climate data for Pond 4 referenced in the mansucript\
**P5_climate_projections.csv** - climate data for Pond 5 referenced in the mansucript\
These data files include values used to identify hydrology (breeding potential) as: 0) unsuitable, 1) suitable with small metamorphs, or 3) suitable with large metamorphs. There are also values determining whether or not there were potential negative interactions between hydrology and phenology in a given breeding season (1 = yes, 0 = no). Values are shown separately for the two study sites.
