# Functions



get_stable_cover = function(data, running_mean, limit, interpolate) {
  
  # Function to estimate upper limits. 
  # - data contains time-series of coral cover 
  # - running_mean is the running_mean window (e.g., 5)
  # - limit is the proportion of the maximum running mean that the coral cover
  # must stabilise within (e.g., 0.8 is 80%)
  # - interpolate (TRUE or FALSE) determines whether coral cover for unsampled 
  # years is interpolated
  
  reefs <- unique(data[ ,grepl("reef", colnames(data), ignore.case = TRUE)])
  
  # Summarise coral cover by site and year
  sdat <- data %>% group_by(REEF_NAME, YEAR_CODE, ZONE) %>%
    summarise("cover" = median(HARD_COVER, na.rm = TRUE),
              "sd" = sd(HARD_COVER, na.rm = TRUE))
  
  sdat <- as.data.frame(sdat)
  
  # Create a column for sampling year
  sdat$Year <- NA
  for (i in 1:nrow(sdat)) {
    sdat$Year[i] <- as.numeric(substr(as.character(sdat$YEAR_CODE[i]), start = 1,
                                      stop = 4)) + 0.5
  }
  
  
  
  # Data frame to save the running mean values
  running <- data.frame("Reef" = factor(), "Zone" = factor(),
                        "Year" = numeric(), "Cover" = numeric(),
                        "Running" = numeric(), "True_cover" = factor(), 
                        "stable" = factor(),
                        "n_years" = numeric(), 
                        "max_cover" = numeric(), "high_cover" = factor())
  
  sdat$Reef_zone <- paste(sdat$REEF_NAME, sdat$ZONE)
  
  # Create a row for each year, even if it wasn't sampled and interpolate values
  # if required
  for(rz in unique(sdat$Reef_zone)) {
    
    sub_dat <- sdat[sdat$Reef_zone == rz, ]
    sub_dat <- sub_dat[order(sub_dat$Year), ]
    
    years <- seq(min(sub_dat$Year), max(sub_dat$Year), by = 1)
    years_m <- sort(sub_dat$Year)
    
    run <- data.frame("Reef" = factor(), "Zone" = factor(),
                      "Year" = numeric(), "Cover" = numeric(),
                      "Running" = numeric())
    
    for (y in years) {
      cover <- NULL
      cover <-  sub_dat[which(years_m == y), ]$cover
      
      if (length(cover) == 0) {
        cover <- NA
      }
      
      run1 <- data.frame("Reef" = sub_dat$REEF_NAME[1],
                         "Zone" = sub_dat$ZONE[1],
                         "Year" = y,
                         "Cover" = cover,
                         "Running" = NA)
      run <- rbind(run, run1)
      
    }
    
    run$True_cover = ifelse(is.na(run$Cover) == TRUE, "no", "yes")
    
    if (interpolate ) {
      # Linearly interpolate missing years
      covers <- is.na(run$Cover)
      true_covers <- which(covers == FALSE)
      
      for (i in 1:nrow(run)) {
        
        if (is.na(run$Cover[i]) == TRUE) {
          prev_cover <- max(true_covers[true_covers < i ])
          next_cover <- min(true_covers[true_covers > i])
          
          delta_year <- (run$Year[next_cover] - run$Year[prev_cover])
          
          run$Cover[i] <- run$Cover[prev_cover] + 
            (run$Cover[next_cover] - run$Cover[prev_cover]) /
            (delta_year) * (run$Year[i] - run$Year[prev_cover])
          
        }
        
      }
    }
    
    
    run$stable <- NA
    run$Running <- NA
    run$n_years <- NA
    
    interval =  running_mean/2 - 0.5
    
    if (nrow(run) > (running_mean - 1)) {
      
      min_int = interval + 1
      max_int = max(min_int, nrow(run) - interval)
      for (i in min_int:max_int) {
        run$Running[i] <- mean(run$Cover[(i-interval) : (i+interval)], na.rm = TRUE)
        IsNA =NULL; IsNA = length(which(run$True_cover[(i-interval):(i+interval)] == "no"))
        run$n_years[i] <- running_mean - IsNA
        
        
      }
    }
    
    # check whether yearly cover is close to the maximum running mean
    run$max_cover <- max(run$Running, na.rm = TRUE)
    run$max_cover <- ifelse(run$max_cover == -Inf, NA, run$max_cover)
    run$high_cover <- ifelse(run$Running / run$max_cover > limit, "yes", "no")
    run$stable <- NA
    
    # check whether coral cover is stable 
    if (nrow(run) >= running_mean) {
      min_int = interval + 1
      max_int = max(min_int, nrow(run) - interval)
      for (i in min_int:max_int) {
        stable = NULL
        stable = run$high_cover[(i-interval):(i+interval)]
        stable = stable[complete.cases(stable)]
        run$stable[i] <- ifelse("no" %in% stable, "no",
                                ifelse(length(stable) < 5, "no", "yes"))
        run$stable[i] <- ifelse(run$stable[i] == "yes" & run$n_years[i] > interval, "yes", "no")
      }
    }
    
    running <- rbind(running, run)
    
    
  }
  
  # Check if any stable covers occur at the beginning or end of the sampling period
  running$start = NA
  running$end = NA
  running$year_num = NA
  running$order = NA
  
  
  
  
  running <- running[complete.cases(running$stable), ]
  
  run1 <- running[running$stable == "yes" & running$n_years > interval, ] %>%
    group_by(Reef,  Zone) %>%
    summarise("carrying" = max(Running))
  
  
  run1 <- as.data.frame(run1)
  run1 = run1[complete.cases(run1$carrying),]
  run1$year_num = NA
  
  
  for (i in 1:nrow(run1)) {
    subd = NULL
    subd = running[running$Reef == run1$Reef[i] &
                     running$Zone == run1$Zone[i], ]
    
    Year_num = NULL; year_num = NULL
    year_num = running[running$Reef == run1$Reef[i] &
                         running$Zone == run1$Zone[i] &
                         running$Running == run1$carrying[i] &
                         running$stable == "yes", ]$year_num
    Year_num = ifelse(year_num %in% c(interval, (nrow(subd) - (interval)):nrow(subd)), "yes", "no")  
    
    if (length(Year_num) == 1) {
      run1$year_num[i] = Year_num
    } else {
      if (length(which(Year_num == "yes")) > 0) {
        run1$year_num[i] = "yes"
      } else{
        run1$year_num[i] = "no"
      }
    }
    
    
  }
  
  
  
  sdat2 <- run1
  
  # Get rid of sites with inconsistent swimming paths
  sdat2 <- sdat2[!(sdat2$Reef == "15077S" & sdat2$Zone == "Back"),  ]
  sdat2 <- sdat2[!sdat2$Reef == "HELSDON REEF", ]
  sdat2 <- sdat2[!(sdat2$Reef == "ST CRISPIN REEF" & sdat2$Zone %in% c("Flank1", "Front")), ]
  sdat2 <- sdat2[!(sdat2$Reef == "THETFORD REEF" & sdat2$Zone == "Flank1"),  ]
  sdat2 <- sdat2[!(sdat2$Reef == "CHICKEN REEF" & sdat2$Zone == "Flank1"),  ]
  sdat2 <- sdat2[!(sdat2$Reef == "DIP REEF" & sdat2$Zone %in% c("Flank1", "Back", "Front")), ]
  sdat2 <- sdat2[!(sdat2$Reef == "RIB REEF" & sdat2$Zone %in% c("Flank1", "Back", "Flank1")), ]
  sdat2 <- sdat2[!(sdat2$Reef == "ERSKINE REEF" & sdat2$Zone %in% c("Front", "Flank1")), ]
  sdat2 <- sdat2[!(sdat2$Reef == "FARQUHARSON REEF (NO 1)" & sdat2$Zone == "Flank1"),  ]
  sdat2 <- sdat2[!sdat2$Reef == "21302S", ]
  sdat2 <- sdat2[!(sdat2$Reef == "22084S" & sdat2$Zone == "Front"),  ]
  sdat2 <- sdat2[!(sdat2$Reef == "ONE TREE REEF" & sdat2$Zone == "Flank1"),  ]
  sdat2 <- sdat2[!(sdat2$Reef == "20354S" & sdat2$Zone == "Flank1"),  ]
  sdat2 <- sdat2[!(sdat2$Reef == "MOORE REEF" & sdat2$Zone == "Flank1"),  ]
  sdat2 <- sdat2[!(sdat2$Reef == "OSBORNE REEF" & sdat2$Zone %in% c("Back", "Flank1")), ]
  sdat2 <- sdat2[!sdat2$Reef == "ASHMORE BANKS (2)", ]
  
  
  # sdat2 contains the upper limits
  # running contains the running means and observed coral covers for each year in
  # all sites
  
  return(list("carrying" = sdat2, "running" = running))
  
}



get_data_quantile = function(data) {
  
  # Function to extract the 95th quantile of the coral cover time-series
  # data for sites with at least 10 data points. It also converts the 
  # environmental variables as z-scores.
  # The input data is the coral cover time-series.
  # It returns two data frames: one with the 95th coral cover quantiles and
  # one with the means and standard deviation values of the environmental data.
  
  dat1 = data %>% group_by(REEF_NAME, ZONE, YEAR_CODE) %>%
    summarise("cover" = mean(HARD_COVER, na.rm = TRUE),
              "latitude" = mean(Central_LAT, na.rm = TRUE),
              "longitude" = mean(Central_LON, na.rm = TRUE))
  
  sites10 = dat1 %>% group_by(REEF_NAME, ZONE) %>%
    summarise("n" = n())
  sites10 = sites10[sites10$n > 9, ]
  
  
  dat = as.data.frame(sites10)
  dat$quantile = NA
  
  
  for ( i in 1:nrow(dat)) {
    dat$quantile[i] = quantile(dat1[dat1$REEF_NAME == dat$REEF_NAME[i] &
                                      dat1$ZONE == dat$ZONE[i], ]$cover, 0.95, na.rm = TRUE) /100
    
    
  }
  
  colnames(dat)[1:2] = c("Reef", "Zone")
  
  
  
  
  # Read environmental data
  env <- read.csv("Data/environmental_data.csv", sep = ",")
  dat <- merge(dat, env, by = c("Reef", "Zone"), all = TRUE)
  
  # Eliminate sites where not all environmental variables are available
  dat <- dat[complete.cases(dat[ , c("quantile",
                                     "temp_median",
                                     "Secchi_median", 
                                     "per_suitable",
                                     "Zone", "Latitude", "Longitude") ]), ]
  
  
  # Normalise environmental variables
  dat$st_temp_median <- (dat$temp_median - mean(dat$temp_median, na.rm = TRUE)) /
    sd(dat$temp_median, na.rm = TRUE)
  dat$st_Secchi <- (dat$Secchi_median  - mean(dat$Secchi_median, na.rm = TRUE)) /
    (sd(dat$Secchi_median, na.rm = TRUE))
  dat$st_habitat <- (dat$per_suitable - mean(dat$per_suitable, na.rm = TRUE)) /
    sd(dat$per_suitable, na.rm = TRUE)
  
  mean_vals <- data.frame("variable" = c("temperature", "Secchi", 
                                         "habitat"),
                          "mean" = c(mean(dat$temp_median),
                                     mean(dat$Secchi_median),
                                     mean(dat$per_suitable)),
                          "sd" = c(sd(dat$temp_median),
                                   sd(dat$Secchi_median),
                                   sd(dat$per_suitable)))
  
  return(list("dat" = dat, "mean_vals" = mean_vals))
}




get_quantile = function(data, model ) {
  
  dat = get_data_quantile(data)$dat
  file_name = paste(path_out, paste(model, "quantile", sep = "_"), sep = "/")
  
  if (paste(paste(model, "quantile", sep = "_"), "rds", sep = ".") %in% list.files(path_out)) {
    m = readRDS(paste(file_name, "rds", sep = "."))
  }  else {
    if (model == "m1") {
      prior <- get_prior(quantile ~ st_habitat + 
                           st_temp_median +
                           st_Secchi +
                           (1| Reef), 
                         data = dat,
                         family = "beta")
      
      m <- brm(quantile ~ st_habitat + 
                 st_temp_median +
                 st_Secchi + 
                 (1| Reef), 
               data = dat,
               family = "beta", chains = 3, cores = 3, 
               prior = prior, iter = 20000, thin = 5,
               control = list(adapt_delta = 0.99, max_treedepth = 17),
               save_pars = save_pars(all = TRUE), file = file_name)
      
      cond = conditional_effects(m, re_formula = NA, robust = TRUE)
      
      
      cond$st_temp_median$Temperature = convert_scales("temperature", cond$st_temp_median$st_temp_median,
                                                       mean_vals)
      
      
      Figure_quantile <- ggplot() +
        geom_ribbon(data = cond$st_habitat, aes( x = habitat * 100, ymin = lower__* 100, ymax = upper__ * 100),  alpha = 0.1) +
        geom_line(data = cond$st_habitat, aes( x = habitat * 100, y = estimate__ * 100), size = 1) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "red", alpha = 0.5, size = 1.4) +
        geom_point(data = dat, aes( x = per_suitable * 100, y = (quantile* 100)), fill = "grey", col = "black",
                   pch = 21,  size = 2) +
        xlab("Hard substrate availability (%)") +
        ylab("Upper limit (%) (95th percentile)") +
        theme_bw() +
        scale_x_continuous(limits = c(0, 100)) +
        scale_y_continuous(limits = c(0, 100)) +
        theme(text = element_text(size = 15 )) 
      ggsave(paste(path_out, "quantile_temp.png", sep = "/"), Figure_quantile, width = 15, height = 15, units = "cm", dpi = 300)
      
      
    }
    
    if (model == "m2") {
      prior <- get_prior(quantile ~ 
                           st_temp_median +
                           st_Secchi +
                           (1| Reef), 
                         data = dat,
                         family = "beta")
      
      m <- brm(quantile ~  
                 st_temp_median +
                 st_Secchi + 
                 (1| Reef), 
               data = dat,
               family = "beta", chains = 3, cores = 3, 
               prior = prior, iter = 20000, thin = 5,
               control = list(adapt_delta = 0.99, max_treedepth = 17),
               save_pars = save_pars(all = TRUE), file = file_name)
      
    }
    if (model == "m3") {
      prior <- get_prior(quantile ~ st_habitat + 
                           st_Secchi +
                           (1| Reef), 
                         data = dat,
                         family = "beta")
      
      m <- brm(quantile ~ st_habitat + 
                 st_Secchi + 
                 (1| Reef), 
               data = dat,
               family = "beta", chains = 3, cores = 3, 
               prior = prior, iter = 20000, thin = 5,
               control = list(adapt_delta = 0.99, max_treedepth = 17),
               save_pars = save_pars(all = TRUE), file = file_name )
      
    }
    if (model == "m4") {
      prior <- get_prior(quantile ~ st_habitat + 
                           st_temp_median +
                           (1| Reef), 
                         data = dat,
                         family = "beta")
      
      m <- brm(quantile ~ st_habitat + 
                 st_temp_median +
                 (1| Reef), 
               data = dat,
               family = "beta", chains = 3, cores = 3, 
               prior = prior, iter = 20000, thin = 5,
               control = list(adapt_delta = 0.99, max_treedepth = 17),
               save_pars = save_pars(all = TRUE), file = file_name)
    }
    if (model == "m5") {
      prior <- get_prior(quantile ~ st_habitat + 
                           (1| Reef), 
                         data = dat,
                         family = "beta")
      
      m <- brm(quantile ~ st_habitat + 
                 (1| Reef), 
               data = dat,
               family = "beta", chains = 3, cores = 3, 
               prior = prior, iter = 20000, thin = 5,
               control = list(adapt_delta = 0.99, max_treedepth = 17),
               save_pars = save_pars(all = TRUE), file = file_name)
    }
    
    if (model == "m6") {
      prior <- get_prior(quantile ~ 
                           st_Secchi +
                           (1| Reef), 
                         data = dat,
                         family = "beta")
      
      m <- brm(quantile ~ 
                 st_Secchi + 
                 (1| Reef), 
               data = dat,
               family = "beta", chains = 3, cores = 3, 
               prior = prior, iter = 20000, thin = 5,
               control = list(adapt_delta = 0.99, max_treedepth = 17),
               save_pars = save_pars(all = TRUE), file = file_name)
      
    }
    if (model == "m7") {
      prior <- get_prior(quantile ~ 
                           st_temp_median +
                           (1| Reef), 
                         data = dat,
                         family = "beta")
      
      m <- brm(quantile ~ 
                 st_temp_median +
                 (1| Reef), 
               data = dat,
               family = "beta", chains = 3, cores = 3, 
               prior = prior, iter = 20000, thin = 5,
               control = list(adapt_delta = 0.99, max_treedepth = 17),
               save_pars = save_pars(all = TRUE), file = file_name)
      
    }
  }
  
  df = data.frame("model" = model, "running" = "", "percentage" = "",
                  "interpolation" = "", "temperature" = NA, "temp_ci" = NA,
                  "habitat" = NA, "habitat_ci" = NA, "Secchi" = NA, 
                  "Secchi_ci" = NA, "n_sites" = NA, "n_reefs" = NA, "type" = "quantile")
  post <- as.data.frame(as_draws_df(m)[ , 2:4])
  
  if(length(which(grepl("temp", colnames(post))) > 0)) {
    df$temperature = quantile(post[, which(grepl("temp", colnames(post)))], 0.5)
    df$temp_ci = ifelse(quantile(post[, which(grepl("temp", colnames(post)))], 0.025) < 0 &
                          quantile(post[, which(grepl("temp", colnames(post)))], 0.975) > 0 ,"n", "y")
  }
  
  if(length(which(grepl("habitat", colnames(post))) > 0)) {
    df$habitat = quantile(post[, which(grepl("habitat", colnames(post)))], 0.5)
    df$habitat_ci = ifelse(quantile(post[, which(grepl("habitat", colnames(post)))], 0.025) < 0 &
                             quantile(post[, which(grepl("habitat", colnames(post)))], 0.975) > 0 ,"n", "y")
  }
  
  
  if(length(which(grepl("Secchi", colnames(post))) > 0)) {
    df$Secchi = quantile(post[, which(grepl("Secchi", colnames(post)))], 0.5)
    df$Secchi_ci = ifelse(quantile(post[, which(grepl("Secchi", colnames(post)))], 0.025) < 0 &
                            quantile(post[, which(grepl("Secchi", colnames(post)))], 0.975) > 0 ,"n", "y")
  }
  df$n_sites = nrow(m$data)
  df$n_reefs = length(unique(m$data$Reef))
  return(df)
}


convert_scales = function(variable, value, mean_vals) {
  
  # Function to convert from normalised to natural scale
  # - variable can be: "habitat" (hard substrate), "temperature", or "Secchi"
  # -value is the normalised value that needs to be converted back to raw
  # scale
  # - mean_vals is the data frame that contains the means and standard deviations
  # for the three environmental variables
  # The function returns a value on its raw scale.
  mean = mean_vals[mean_vals$variable == variable, ]$mean
  sd = mean_vals[mean_vals$variable == variable, ]$sd
  
  return(value*sd + mean)
  
}
