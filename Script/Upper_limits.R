rm(list = c())

path = "C:/Users/malvarez/OneDrive - Australian Institute of Marine Science/Documents/Recovery/Manuscript/Global Ecology and Biogeography/Resubmission/To Submit/Data"
path_out = "C:/Users/malvarez/OneDrive - Australian Institute of Marine Science/Documents/Recovery/Manuscript/Global Ecology and Biogeography/Resubmission/To Submit/Outputs"

set.seed(9)

library(ggplot2)
library(brms)
library(dplyr)
library(viridis)
library(SpatialEpi)
library(reshape2)
library(egg)
library(patchwork)
library(maps)
library(sf)
library(spdep)
library(ggstance)
library(deldir)
library(sp)
library(purrr)


# Read data from the manta tow surveys
manta <- read.csv(paste(path, "coral_cover.csv", sep = "/"))



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
  
  for (reef in unique(running$Reef)) {
    sub_d = NULL
    sub_d = running[running$Reef == reef, ]
    
    for (zone in unique(sub_d$Zone)) {
      sub_d = running[running$Reef == reef, ]
      sub_d = sub_d[sub_d$Zone == zone, ]
      nrow_subd = nrow(sub_d)
      
      running[running$Reef == reef & running$Zone == zone, ]$year_num = 1:nrow(running[running$Reef == reef & running$Zone == zone, ])
      
      Order = order(running[running$Reef == reef & running$Zone == zone, ]$Running)
      running[running$Reef == reef & running$Zone == zone, ]$order = Order
      
      if (nrow_subd > (running_mean + 2)){ 
        
        start_condition = 0
        start_condition = sub_d$stable[interval +1] == "yes" & 
          sub_d$Running[1:interval] > sub_d$Running[interval + 1]
        start_condition = ifelse(start_condition == TRUE, 1, 0)
        start_condition = sum(start_condition)/length(start_condition)
        
        if (start_condition == 1) {
          running[running$Reef == reef & running$Zone == zone, ]$start = "yes"
        }
        
        
        end_condition = 0
        end_condition = sub_d$stable[nrow_subd - interval] == "yes" & 
          sub_d$Running[(nrow_subd-interval +1):nrow_subd] > sub_d$Running[nrow_subd - interval]
        
        end_condition = ifelse(end_condition == TRUE, 1, 0)
        end_condition = sum(end_condition)/length(end_condition)
        
        if (end_condition == 1) {
          
          running[running$Reef == reef & running$Zone == zone, ]$end = "yes"
        }
        
        
      }
      
    }
  }
  
  
  
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


# Get upper limits with a running mean of five years and with interpolation for
# missing years
dat = get_stable_cover(manta, 5, 0.80, TRUE)[[1]]


# Plot upper limits if needed (change next line to TRUE)
plot_ul = FALSE

if (plot_ul) {
  df = get_stable_cover(manta, 5, 0.80, TRUE)[[2]]
  
  
  pdf(paste(path, "Upper limit plots.pdf", sep = "/"))
  for (reef in unique(df$Reef)) {
    
    
    p = NULL
    p = ggplot() +
      geom_point(data = df[df$Reef == reef & df$True_cover == "yes", ], aes(x = Year, y = Cover), pch = 21, fill = "grey", col = "black") +
      theme_bw() +
      facet_wrap(~Zone, ncol = 2) +
      geom_line(data = df[df$Reef == reef , ],
                aes(x = Year, y = Running)) +
      geom_hline(data = dat[dat$Reef == reef, ], aes(yintercept = carrying), linetype = "dashed",
                 col = "red") +
      ggtitle(reef)
    
    print(p)
    print(reef)
    
  }
  
  dev.off()
}


# Read environmental data
env <- read.csv(paste(path, "environmental_data.csv", sep = "/"))
dat <- merge(dat, env, by = c("Reef", "Zone"), all = TRUE)
dat <- dat[!is.na(dat$carrying), ]
dat$carrying_prop <- dat$carrying / 100


# Save all upper limit estimates in a separate data frame 
dat1 <- dat


# Check correlation between hard substrate availability and water velocity
cor.test(dat1$per_suitable, dat1$median_ubed90, method = "spearman")


ggplot(dat, aes(x = median_ubed90, y = per_suitable * 100)) +
  geom_point(pch = 21, fill = "grey", size = 3) +
  theme_bw() +
  xlab(expression(paste("Water velocity at the bed (m", s^-1, ", 9", 0^"th", " quantile)"))) +
  ylab("Hard substrate available (%)") +
  theme(text = element_text(size = 15))


# Eliminate sites where not all environmental variables are available
dat <- dat[complete.cases(dat[ , c(
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

# Save mean and standard deviation values for each variable
mean_vals <- data.frame("variable" = c("temperature", "Secchi", 
                                       "habitat"),
                        "mean" = c(mean(dat$temp_median),
                                   mean(dat$Secchi_median),
                                   mean(dat$per_suitable)),
                        "sd" = c(sd(dat$temp_median),
                                 sd(dat$Secchi_median),
                                 sd(dat$per_suitable)))


# Get distances between reefs
mat <- matrix(NA, nrow = nrow(dat), ncol = 2)
mat[ ,1] <- dat$Longitude
mat[ ,2] <- dat$Latitude

xy <- latlong2grid(mat)
dat$x <- xy$x
dat$y <- xy$y

data <- dat[ , c("Reef", "Zone", "x", "y")]
data <- data[!duplicated(data), ]


data <- data %>% group_by(Reef) %>%
  summarise("x" = mean(x, na.rm = TRUE),
            "y" = mean(y, na.rm = TRUE))

data <- as.data.frame(data)
reefs <- unique(data[order(data$y), ]$Reef)

spat <- matrix(NA, nrow = length(reefs), ncol = length(reefs))
colnames(spat) <- reefs
row.names(spat) <- reefs
spat2 <- spat

for ( i in 1: nrow(spat)) {
  
  d1 <- data[data$Reef == rownames(spat)[i], ]
  
  x1 <- d1$x ; y1 <- d1$y
  
  for (j in 1: ncol(spat)) {
    
    d2 <- data[data$Reef == colnames(spat)[j], ]
    
    x2 <- d2$x ; y2 <- d2$y
    
    spat[i,j] <- sqrt((x1 - x2) ^ 2 + (y1 - y2) ^ 2) 
    spat2[i,j] <-  sqrt((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
    
    
    spat[i,j] <- ifelse(spat[i,j] > 0,  1/(spat[i,j])^ 0.5, 1 )
    
    if(spat[i,j] < 0 ) {print(i); print(j)}
    
  }
}


melted <- melt(spat)
colnames(melted)[1:2] <- c("X1", "X2")
ggplot(melted) +
  geom_tile(aes(x = X1, y = X2, fill = value)) +
  scale_fill_viridis(option = "B") 


dat$Reef_spat <- dat$Reef

dat$Reef_zone <- paste(dat$Reef, dat$Zone, sep = "-")


# Plot spatial patterns in upper limits
dat$Zone <- ifelse(dat$Zone == "Flank1", "Flank",
                   ifelse(dat$Zone == "Flank2", "Flank",
                          dat$Zone))


world_map <- map_data("world")
Aus <- world_map %>% filter(region == "Australia")

sum_dat1 = dat1 %>% group_by(Reef) %>%
  summarise("Latitude" = mean(Latitude), "Longitude" = mean(Longitude),
            "carrying_prop" = mean(carrying_prop), "n" = n())
sum_dat1 = as.data.frame(sum_dat1)

sizes = c(3,4,5,6)
names(sizes) = c("1", "2", "3", "4")

Figure2a <- 
  ggplot() +
  geom_polygon(data = Aus, aes(x = long, y = lat, group = group), fill = "grey",
               col = "black") +
  coord_map(xlim = c(142, 153.5), ylim = c(-25, -9)) +
  geom_point(data = sum_dat1, aes(x = Longitude, y = Latitude,
                                  fill =  carrying_prop * 100, size = as.factor(n)),
             pch = 21, col = "black", alpha = 0.5) +
  scale_fill_viridis(limits = c(0, 80),
                     breaks = seq(0, 100, by = 20),
                     labels = paste0(seq(0, 100, by = 20), "%"),
                     option = "B",
                     name = "") +
  scale_size_manual(values = sizes, name = "Sites per reef") +
  theme_bw() +
  ylab("") +
  xlab("") +
  scale_y_continuous(breaks = seq(-25, -10, by = 5),
                     labels = c("25°S", "20°S", "15°S", "10°S")) +
  scale_x_continuous(breaks = seq(144, 153, by = 3),
                     labels = c("144°E", "147°E", "150°E", "153°E")) +
  theme_bw() +
  theme(text = element_text(size = 14 ), plot.title = element_text(size = 16)) +
  ggtitle("a- Among reef variation")

dat1$Zone_grouped <- ifelse(dat1$Zone %in% c("Flank1" ,"Flank2"), "Flank", dat1$Zone)

Figure2b <- ggplot() +
  geom_jitter(data = dat1, aes(y = carrying_prop * 100, x = Zone_grouped), height = 0, width = 0.2,
              pch = 21, fill = "grey", alpha = 0.4) +
  geom_boxplot(data = dat1, aes(y = carrying_prop * 100, x = Zone_grouped), fill = NA) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), 
                     limits = c(0, 80)) + 
  scale_x_discrete(breaks = c("Back", "Flank", "Front"),
                   labels = c("Sheltered", 
                              "Intermediate",
                              "Exposed")) +
  xlab("Wave exposure") +
  ylab("Upper limit (%)") +
  theme_bw() +
  theme(text = element_text(size = 14 ), plot.title = element_text(size = 16)) +
  ggtitle("b - Within reef variation")


Figure2 <- Figure2a + Figure2b


# Plot within reef variation
lat_reefs = dat1[order(dat1$Latitude, decreasing = TRUE), ]$Reef
lat_reefs = lat_reefs[!duplicated(lat_reefs)]

FigureS5a = ggplot(dat1, aes(x = Reef, y = carrying*100)) +
  geom_boxplot() +
  scale_x_discrete(limits = lat_reefs) +
  ylab("Upper limits (%)") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("a")

range_sum = dat1 %>% group_by(Reef) %>%
  summarise("n" = n(),
            "Range" = max(carrying) - min(carrying))
range_sum = as.data.frame(range_sum)


FigureS5b = ggplot(range_sum, aes(x = n, y = Range)) +
  geom_boxplot(aes(group = n)) +
  theme_bw() +
  xlab("Number of sites with upper limit estimates per reef") +
  ylab("Range of upper limit \n estimates within reef (%)") +
  ggtitle("b")


FigureS5 <- FigureS5a / FigureS5b

# Make sure within-group covariance matrices are positive definite
spat_new = diag(diag(spat + 0.000001))
colnames(spat_new) = colnames(spat)
row.names(spat_new) = row.names(spat)

# Fit models
prior <- get_prior(carrying_prop ~ st_habitat *
                     st_temp_median *
                     st_Secchi +
                     (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat_new),
                   family = "beta")

m1 <- brm(carrying_prop ~ st_habitat * st_temp_median *
            st_Secchi +
            (1| gr(Reef_spat, cov = spat)),
          data = dat,
          data2 = list(spat = spat_new),
          family = "beta", chains = 3, cores = 3, 
          prior = prior, iter = 10000, thin = 5,
          control = list(adapt_delta = 0.99, max_treedepth = 17),
          save_pars = save_pars(all = TRUE), file = paste(path_out, "m1.rds", sep = "/"))
m1 <- add_criterion(m1, "loo", moment.match = TRUE)



prior <- get_prior(carrying_prop ~ st_habitat * st_temp_median +
                     st_habitat * st_Secchi +
                     st_Secchi * st_temp_median +
                     (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat_new),
                   family = "beta")

m2 <- brm(carrying_prop ~ st_habitat * st_temp_median +
            st_habitat * st_Secchi +
            st_Secchi * st_temp_median +
            (1| gr(Reef_spat, cov = spat)),
          data = dat,
          data2 = list(spat = spat_new),
          family = "beta", chains = 3, cores = 3, 
          prior = prior, iter = 10000, thin = 5,
          control = list(adapt_delta = 0.99, max_treedepth = 17),
          save_pars = save_pars(all = TRUE), file = paste(path_out, "m2.rds", sep = "/"))
m2 <- add_criterion(m2, "loo", moment.match = TRUE)




prior <- get_prior(carrying_prop ~ st_habitat * st_temp_median +
                     st_habitat * st_Secchi +
                     (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat_new),
                   family = "beta")

m3 <- brm(carrying_prop ~ st_habitat * st_temp_median +
            st_habitat * st_Secchi +
            (1| gr(Reef_spat, cov = spat)),
          data = dat,
          data2 = list(spat = spat_new),
          family = "beta", chains = 3, cores = 3, 
          prior = prior, iter = 10000, thin = 5,
          control = list(adapt_delta = 0.99, max_treedepth = 17),
          save_pars = save_pars(all = TRUE), file = paste(path_out, "m3.rds", sep = "/"))
m3 <- add_criterion(m3, "loo", moment.match = TRUE)




prior <- get_prior(carrying_prop ~ st_habitat * st_temp_median +
                     st_Secchi * st_temp_median +
                     (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat_new),
                   family = "beta")

m4 <- brm(carrying_prop ~ st_habitat * st_temp_median +
            st_Secchi * st_temp_median +
            (1| gr(Reef_spat, cov = spat)),
          data = dat,
          data2 = list(spat = spat_new),
          family = "beta", chains = 3, cores = 3, 
          prior = prior, iter = 10000, thin = 5,
          control = list(adapt_delta = 0.99, max_treedepth = 17),
          save_pars = save_pars(all = TRUE), file = paste(path_out, "m4.rds", sep = "/"))
m4 <- add_criterion(m4, "loo", moment.match = TRUE)



prior <- get_prior(carrying_prop ~ 
                     st_habitat * st_Secchi +
                     st_Secchi * st_temp_median +
                     (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat_new),
                   family = "beta")

m5 <- brm(carrying_prop ~ 
            st_habitat * st_Secchi +
            st_Secchi * st_temp_median +
            (1| gr(Reef_spat, cov = spat)),
          data = dat,
          data2 = list(spat = spat_new),
          family = "beta", chains = 3, cores = 3, 
          prior = prior, iter = 10000, thin = 5,
          control = list(adapt_delta = 0.99, max_treedepth = 17),
          save_pars = save_pars(all = TRUE), file = paste(path_out, "m5.rds", sep = "/"))
m5 <- add_criterion(m5, "loo", moment.match = TRUE)


prior <- get_prior(carrying_prop ~ st_habitat * st_temp_median +
                     st_Secchi +
                     (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat_new),
                   family = "beta")

m6 <- brm(carrying_prop ~ st_habitat * st_temp_median +
            st_Secchi +
            (1| gr(Reef_spat, cov = spat)),
          data = dat,
          data2 = list(spat = spat_new),
          family = "beta", chains = 3, cores = 3, 
          prior = prior, iter = 10000, thin = 5,
          control = list(adapt_delta = 0.99, max_treedepth = 17),
          save_pars = save_pars(all = TRUE), file = paste(path_out, "m6.rds", sep = "/"))
m6 <- add_criterion(m6, "loo", moment.match = TRUE)





prior <- get_prior(carrying_prop ~  st_temp_median +
                     st_habitat * st_Secchi +
                     (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat_new),
                   family = "beta")

m7 <- brm(carrying_prop ~ st_temp_median +
            st_habitat * st_Secchi +
            (1| gr(Reef_spat, cov = spat)),
          data = dat,
          data2 = list(spat = spat_new),
          family = "beta", chains = 3, cores = 3, 
          prior = prior, iter = 10000, thin = 5,
          control = list(adapt_delta = 0.99, max_treedepth = 17),
          save_pars = save_pars(all = TRUE),file = paste(path_out, "m7.rds", sep = "/"))
m7 <- add_criterion(m7, "loo", moment.match = TRUE)


prior <- get_prior(carrying_prop ~ st_habitat +
                     st_Secchi * st_temp_median +
                     (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat_new),
                   family = "beta")

m8 <- brm(carrying_prop ~ st_habitat +
            st_Secchi * st_temp_median +
            (1| gr(Reef_spat, cov = spat)),
          data = dat,
          data2 = list(spat = spat_new),
          family = "beta", chains = 3, cores = 3, 
          prior = prior, iter = 10000, thin = 5,
          control = list(adapt_delta = 0.99, max_treedepth = 17),
          save_pars = save_pars(all = TRUE), file = paste(path_out, "m8.rds", sep = "/"))
m8 <- add_criterion(m8, "loo", moment.match = TRUE)



prior <- get_prior(carrying_prop ~ st_habitat + st_temp_median +
                     st_Secchi +
                     (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat_new),
                   family = "beta")

m9 <- brm(carrying_prop ~ st_habitat + st_temp_median +
            st_Secchi +
            (1| gr(Reef_spat, cov = spat)),
          data = dat,
          data2 = list(spat = spat_new),
          family = "beta", chains = 3, cores = 3, 
          prior = prior, iter = 10000, thin = 5,
          control = list(adapt_delta = 0.99, max_treedepth = 17),
          save_pars = save_pars(all = TRUE), file = paste(path_out, "m9.rds", sep = "/"))
m9 <- add_criterion(m9, "loo", moment.match = TRUE)

loo_compare(m1, m2, m3, m4, m5, m6, m7, m8, m9)




prior <- get_prior(carrying_prop ~  st_temp_median +
                     st_habitat * st_Secchi +
                     (1| Reef),
                   data = dat,
                   family = "beta")

m7b <- brm(carrying_prop ~ st_temp_median +
             st_habitat * st_Secchi +
             (1| Reef),
           data = dat,
           family = "beta", chains = 3, cores = 3, 
           prior = prior, iter = 10000, thin = 5,
           control = list(adapt_delta = 0.99, max_treedepth = 17),
           save_pars = save_pars(all = TRUE), file = paste(path_out, "m7b.rds", sep = "/"))
m7b <- add_criterion(m7b, "loo", moment.match = TRUE)

loo_compare(m7, m7b)


## Example 1
example_preds <- data.frame("temperature" = c(25, 26)
)
example_preds$st_temp_median <- (example_preds$temperature - mean_vals[mean_vals$variable == "temperature", ]$mean)/
  mean_vals[mean_vals$variable == "temperature", ]$sd
example_preds$st_Secchi <- 0
example_preds$st_habitat <- 0

example_preds <- cbind(example_preds, as.data.frame(fitted(m7, example_preds, 
                                                           re_formula = NA)))

## Example 2
example_preds <- expand.grid("habitat" = c(0.5, 0.6),
                             "Secchi" = c(15, 25) )

example_preds$st_temp_median <- 0
example_preds$st_Secchi <- (example_preds$Secchi - mean_vals[mean_vals$variable == "Secchi", ]$mean)/
  mean_vals[mean_vals$variable == "Secchi", ]$sd
example_preds$st_habitat <- (example_preds$habitat - mean_vals[mean_vals$variable == "habitat", ]$mean)/
  mean_vals[mean_vals$variable == "habitat", ]$sd

example_preds <- cbind(example_preds, as.data.frame(fitted(m7, example_preds, 
                                                           re_formula = NA)))

# Plot model predictions
dat$rand_m1 <- as.data.frame(fitted(m7, dat))$Estimate - 
  as.data.frame(fitted(m7, dat, re_formula = NA))$Estimate  


cond = conditional_effects(m7, re_formula = NA, robust = TRUE)


convert_scales = function(variable, value, mean_vals) {
  
  # Function to convert from normalised to natural scale
  mean = mean_vals[mean_vals$variable == variable, ]$mean
  sd = mean_vals[mean_vals$variable == variable, ]$sd
  
  return(value*sd + mean)
  
}

cond$st_temp_median$Temperature = convert_scales("temperature", cond$st_temp_median$st_temp_median,
                                                 mean_vals)

Figure3b <- ggplot() +
  geom_point(data = dat, aes( x = temp_median, y = (carrying_prop - rand_m1) * 100), col = "black",
             pch = 21, fill = "grey", size = 2) +
  geom_ribbon(data = cond$st_temp_median, aes( x = Temperature, ymin = lower__* 100, ymax = upper__ * 100), alpha = 0.2) +
  geom_line(data = cond$st_temp_median, aes( x = Temperature, y = estimate__ * 100), size = 1.5) +
  xlab("Temperature (°C)") +
  ylab("Upper limit (%)") +
  theme_bw() +
  theme(text = element_text(size = 14 )) +
  ggtitle("b")


cond$`st_habitat:st_Secchi`$Secchi = convert_scales("Secchi", cond$`st_habitat:st_Secchi`$st_Secchi,
                                                    mean_vals)
cond$`st_habitat:st_Secchi`$habitat = convert_scales("habitat", cond$`st_habitat:st_Secchi`$st_habitat,
                                                     mean_vals)

Figure3c <- ggplot() +
  geom_ribbon(data = cond$'st_habitat:st_Secchi', aes( x = habitat * 100, ymin = lower__* 100, ymax = upper__ * 100, group = Secchi, fill = Secchi), alpha = 0.1) +
  geom_line(data = cond$'st_habitat:st_Secchi', aes( x = habitat * 100, y = estimate__ * 100, group = Secchi, col = Secchi), size = 1) +
  geom_point(data = dat, aes( x = per_suitable * 100, y = (carrying_prop - rand_m1) * 100, fill = Secchi_median), col = "black",
             pch = 21,  size = 2) +
  xlab("Hard substrate availability (%)") +
  ylab("Upper limit (%)") +
  scale_fill_viridis(limits = c(12, 28), name = "Secchi depth (m)") + 
  scale_colour_viridis(limits = c(12, 28), name = "Secchi depth (m)") + 
  theme_bw() +
  theme(text = element_text(size = 14 )) +
  ggtitle("c")



post <- as.data.frame(as_draws_df(m7)[ , 2:5])
colnames(post) <- c("temp", "habitat", "Secchi", "habitat_Secchi")


Figure3a <- ggplot() +
  geom_violin(data = post, aes(x = 0, y = temp), fill = "grey" , alpha = 0.2) +
  geom_linerange(data = post, aes(x = 0, ymin = quantile(temp, 0.025), ymax = quantile(temp, 0.975)), 
                 size = 1.1) +
  geom_point(data = post, aes(x = 0, y = quantile(temp, 0.5)), 
             size = 5) +
  
  geom_violin(data = post, aes(x = 1, y = habitat), fill = "grey" , alpha = 0.2) +
  geom_linerange(data = post, aes(x = 1, ymin = quantile(habitat, 0.025), ymax = quantile(habitat, 0.975)), 
                 size = 1.1) +
  geom_point(data = post, aes(x = 1, y = quantile(habitat, 0.5)), 
             size = 5) +
  
  geom_violin(data = post, aes(x = 2, y = Secchi), fill = "grey", alpha = 0.2) +
  geom_linerange(data = post, aes(x = 2, ymin = quantile(Secchi, 0.025), ymax = quantile(Secchi, 0.975)), 
                 size = 1.1) +
  geom_point(data = post, aes(x = 2, y = quantile(Secchi, 0.5)), 
             size = 5) +
  geom_violin(data = post, aes(x = 3, y = habitat_Secchi), fill = "grey", alpha = 0.2) +
  geom_linerange(data = post, aes(x = 3, ymin = quantile(habitat_Secchi, 0.025), ymax = quantile(habitat_Secchi, 0.975)), 
                 size = 1.1) +
  geom_point(data = post, aes(x = 3, y = quantile(habitat_Secchi, 0.5)), 
             size = 5) +
  
  geom_hline(yintercept = 0, linetype ="dashed", size = 1, col = "red") +
  xlab("") +
  theme_bw() +
  scale_x_continuous(breaks = 0:3, labels = c("Temperature", "Hard substrate","Secchi depth", "Hard substrate * Secchi depth"))+
  theme(text = element_text(size = 14 ), axis.text.x = element_text(size = 14)) +
  ggtitle("a") +
  ylab("Coefficient estimates")


Figure3 <- Figure3a /  (Figure3b + Figure3c + ylab("")) 


# Plot recovery with different upper limits

logistic_eq <- function(K, N, R) {
  # Logistic equation in discrete time
  Nt1 <- N + R * (N * (K - N) / K)
  return(list("N" = Nt1))
  
}


# Values to explore
carrying_prop <- seq(0.01, 1, by = 0.01)
coral_cover <- seq(0.02, 1, by = 0.01)
R <- c(0.25, 0.5)


# Create a table to save how many years it takes to reach a certain
# coral cover from 1% 
df3 <- expand.grid("carrying_prop" = numeric(),
                   "coral_cover" = numeric(),
                   "R" = numeric(),
                   "years_to" = integer()) 

Time <- 1:100

for (cc in carrying_prop) {
  
  for (r in R) {
  
    
    N <- 0.01
    
    for (t in Time[2:length(Time)]) {
      N[t] <- as.numeric(logistic_eq(K = cc, N = N[t - 1], R = r))
    }
    
    traj <- data.frame( "Time" = Time - 1, "N" = N)
    
    years <- c()
    
    for (y in 1:length(coral_cover)) {
      years_above <- which(round(traj$N, 2) >= coral_cover[y])
      years[y] <- ifelse(length(years_above) > 0 , min(years_above),
                         NA)
    }
    
    df3_sub <- data.frame("carrying_prop" = cc,
                          "coral_cover" = coral_cover,
                          "R" = r,
                          "years_to" = years) 
    
    df3 <- rbind(df3, df3_sub)
    
    
  }
  print(which(carrying_prop == cc))
}




Figure4a <- ggplot()+
  geom_tile(data = df3[df3$R == R[1], ], aes(x = carrying_prop * 100, y = coral_cover * 100, fill = years_to) ) +
  theme_classic() +
  scale_fill_viridis(option = "C",
                     name = "Years",  na.value = "black", limits = c(min(df3$years_to, na.rm = TRUE), max(df3$years_to, na.rm = TRUE))) +
  ylab("Coral cover (%)") +
  xlab("Upper limit (%)") +
  theme(text = element_text(size = 13),
        plot.title = element_text(size=14)) +
  ggtitle(expression(paste("a - ", italic(R), " = 0.25"))) +
  theme(strip.background =element_rect(fill="white")) +
  scale_x_continuous(limits = c(0,100))



Figure4b <- ggplot()+
  geom_tile(data = df3[df3$R == R[2], ], aes(x = carrying_prop * 100, y = coral_cover * 100, fill = years_to) ) +
  theme_classic() +
  scale_fill_viridis(option = "C",
                     name = "Years",  na.value = "black", limits = c(min(df3$years_to, na.rm = TRUE), max(df3$years_to, na.rm = TRUE))) +
  ylab("Coral cover (%)") +
  xlab("Upper limit (%)") +
  theme(text = element_text(size = 13),
        plot.title = element_text(size=14)) +
  ggtitle(expression(paste("b - ", italic(R), " = 0.5"))) +
  theme(strip.background =element_rect(fill="white")) +
  scale_x_continuous(limits = c(0,100))


Figure4 <- (Figure4a + xlab("")) / (Figure4b + theme(legend.position = "none"))
