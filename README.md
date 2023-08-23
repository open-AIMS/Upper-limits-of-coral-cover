# Upper limits of coral communities
This repository containds the data and code to reproduce the results of the manuscript "Upper limits in the cover of coral communities across spatial scales".

The scripts were run in R version 4.0.5 using the following packages:

- brms (version 2.18.0)
- patchwork (version 1.1.1)
- viridis (version 0.6.2)
- dplyr (version 1.1.1)
- ggplot2 (version 3.4.2)
- reshape2 (version 1.4.4)

    # Data files
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
     + median_ubed90: 90<sup>th</sup>  quantile of bottom stress (summarised as the median value for each site)
     + per_suitable: proportion of hards substrate available
     + PAR8_median: median daily integrated photosynthetically active radiation (PAR) at 8m 
     + temp_median: median temperature

   - realised_upper_limits.csv: contains the estimates of realised upper limits for each site (wave exposure zone within reef)
     + Reef: reef
     + Zone: wave exposure zone (*flank1* and *flank2* have intermediate wave exposure, *front* has high wave exposure, and *back* is sheltered)
     + carrying: realised upper limit estimate (in percentage)
       
   - rds files: contain the fitted models, can be called in brms (using the file function) to avoid rerunning the models.


     # Scripts
     - Get_realised_upper_limits.R: summarises the manta tow data by site (wave exposure zone within reef) for each year, interpolates for missing years, calculates the 5yr-running mean, and indentifies instances when coral cover estabilises around its maximum level. This script generates the file realised_upper_limits.csv.
     - Fit_models.R: fits models to predict 1) realised upper limits as a function of potential upper limits, and 2) realised upper limits as a function of temperature, benthic stress, and light (when the percentage of hard substrate is at least 80%). This scripts generates the figures and results shown in the manuscritp.
