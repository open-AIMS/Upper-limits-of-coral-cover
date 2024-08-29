# Upper limits of coral communities
This repository contains the data and code to reproduce the results of the manuscript "Spatial variation in upper limits of coral cover on the Great Barrier Reef".

The scripts were run in R version 4.3.2 using the following packages:

- brms (version 2.20.4)
- patchwork (version 1.2.0)
- viridis (version 0.6.4)
- dplyr (version 1.1.4)
- ggplot2 (version 3.5.1)

To run the code, the working directory must have three folders: *Data* (which should *coral_cover.csv* and *environmental_variables.csv*), *Code* (which should contain *analyses.R* and *functions.R*), and *Outputs* (where models and plots will be saved).

    ## Data files ##
   - coral_cover.csv: contains coral cover (Scleractinian corals) estimates from the manta tow surveys. Each row corresponds to one tow.
     + YEAR_CODE: financial year when the survey was conducted (e.g., 201415 is for 2014-2015)
     + REEF_NAME: reef
     + ZONE: wave exposure zone (*flank1* and *flank2* have intermediate wave exposure, *front* has high wave exposure, and *back* is sheltered)
     + HC_MIN: lower end of the hard coral cover interval
     + HC_MAX: upper end of the hard coral cover interval
     + Central_LAT: latitude at the center of the tow
     + Central_LON: longitude at the center of the tow
     + HARD_COVER: midpoint of the coral cover interval
       
   - environmental_variables.csv: contains the environmental variables for each site (wave exposure zone within reef)
     + Reef: reef
     + Zone: wave exposure zone (*flank1* and *flank2* have intermediate wave exposure, *front* has high wave exposure, and *back* is sheltered)
     + Latitude: site latitude
     + Longitude: site longitude
     + median_ubed90: 90<sup>th</sup>  quantile of horizontal water velocity at bed (m s<sup>-1</sup>) (summarised as the median value for each site)
     + per_suitable: proportion of hards substrate available
     + temp_median: median temperature (degrees C)
     + Secchi_median: median Secchi depth (m)


     ## Code ##
     - analyses.R:
       + estimates upper limits from time series coral cover data
       + fits models
       + generates figures shown in the manuscript
     - functions.R:
       + includes the required functions to run the analyses


