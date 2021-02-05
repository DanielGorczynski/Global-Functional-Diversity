# Global-Functional-Diversity README
A variety of factors can affect the biodiversity of tropical mammal communities,
but their relative importance and directionality remain uncertain.
Previous global investigations of mammal functional diversity have relied
on range maps instead of observational data to determine community
composition. We test the effects of species pools, habitat heterogeneity, primary
productivity and human disturbance on the functional diversity
(dispersion and richness) of mammal communities using the largest standardized
tropical forest camera trap monitoring system, the Tropical Ecology
Assessment and Monitoring (TEAM) Network. We use occupancy values
derived from the camera trap data to calculate occupancy-weighted functional
diversity and use Bayesian generalized linear regression to
determine the effects of multiple predictors. Mammal community functional
dispersion increased with primary productivity, while functional richness
decreased with human-induced local extinctions and was significantly
lower in Madagascar than other tropical regions. The significant positive
relationship between functional dispersion and productivity was evident
only when functional dispersion was weighted by species’ occupancies.
Thus, observational data from standardized monitoring can reveal the
drivers of mammal communities in ways that are not readily apparent
from range map-based studies. The positive association between
occupancy-weighted functional dispersion of tropical forest mammal communities
and primary productivity suggests that unique functional traits
may be more beneficial in more productive ecosystems and may allow
species to persist at higher abundances.

# Description of all code and data files included in the repository
## Code
Global FD Code ProcB.R - All code used in the manuscript, including calculation of functional diversity metrics, Bayesian regressions, and figures

## Data
All TEAM Occupancy.csv- Data frame containing median occupancy estimates for all species monitored at TEAM sites over the years that surveys were conducted at each site. Formatted for use in Global FD Code ProcB.R

Full trait list bi 1 kg.csv- Data frame containing all trait data of all species found in communities or species pools of the TEAM sites. Body mass and litter size are continous variables, diet and activity period are categorical binary variables, and social group and substrate use are binary. Formatted for use in Global FD Code ProcB.R

Hull.Color.csv- Small data frame giving the R color codes used for to color the convex hulls of the trait space figures as part of figure 2 in the manuscript. 

Site Coords.csv- Small data frame containing the latitude and longitude of each of the TEAM sites, which are used in making the map as a part of figure 2 in the manuscript

Site Trait and Occ data- folder containing data frames listing the communities and species pools at each of the TEAM sites. These lists were used to construct the trait spaces in Figure 2. Formatted for use in Global FD Code ProcB.R

Species1.csv- Data frame listing the species pools of each of the TEAM sites. Formatted for use in Global FD Code ProcB.R

TEAM Site Funct bi.csv- Data frame containing site level variables for use in the Bayesian regression component of the analysis. The data frame includes the different functional diversity response variables and the historical, environmental and anthropogenic predictor variables.


### License: NA

