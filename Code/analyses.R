# This code runs the analyses to estimate the upper limits of coral cover 
# across environmental gradients in temperature, water clarity, and hard
# substrate availability. It reproduces the analyses and figures in the manuscript:
# "Spatial variation in upper limits of coral cover on the Great Barrier Reef".
# Global Ecology and Biogeography.


rm(list = c())

set.seed(9)

# R version 4.3.2 (2023-10-31 ucrt)
library(ggplot2)   # ggplot2_3.5.1 
library(brms)      # brms_2.20.4
library(dplyr)     # dplyr_1.1.4    
library(viridis)   # viridis_0.6.4 
library(egg)       # egg_0.4.5 
library(patchwork) # patchwork_1.2.0
library(maps)      # maps_3.4.2 
library(spdep)     # spdep_1.3-1 
library(ggstance)  # ggstance_0.3.6
library(deldir)    # deldir_2.0-2
library(sp)        # sp_2.1-2
library(purrr)     # purrr_1.0.2 
library(ggridges)  # ggridges_0.5.6

# Folder to save outputs
path_out = "Outputs"

# Read data from the manta tow surveys
manta <- read.csv("Data/coral_cover.csv", sep = ",")

# load required functions
source("Code/functions.R")


################################################################################
#  How do hard substrate availability, temperature, and  water clarity relate  #
#  to upper limits - sensitivity analyses (start)                              #
################################################################################


# Create a data frame to save coefficient estimates from the different models
# (with different explanatory variables and different upper limit definitions)
coef_df = expand.grid("model" = paste("m", 1:7, sep = ""),
                      "running" = seq(5, 11, by = 2),
                      "percentage" = c(80, 85, 90),
                      "interpolation" = c(TRUE),
                      "temperature" = NA,
                      "temp_ci" = NA,
                      "habitat" = NA,
                      "habitat_ci" = NA,
                      "Secchi" = NA,
                      "Secchi_ci" = NA,
                      "n_sites" = NA,
                      "n_reefs" = NA)

# Fill in the data frame with the coefficients
for ( i in 1:nrow(coef_df)) {
  
  # Get upper limits with a running mean of five years and with interpolation for
  # missing years
  
  # Generate a file name for each model
  file_name = paste(path_out, paste(coef_df[i, 1:4], collapse = "_"), sep = "/")
  
  
  # If the model has already run, load the saved rds file. Otherwise, run the model
  if (paste(paste(coef_df[i, 1:4], collapse = "_"), "rds", sep = ".") %in% list.files(path_out)) {
    m = readRDS(paste(file_name, "rds", sep = "."))
  } else {
    
    # Extract upper limits according to the running-mean window and the 
    # threshold for the difference between the running mean and the maximum
    # coral cover
    dat = get_stable_cover(manta, coef_df$running[i], coef_df$percentage[i]/100,  
                           coef_df$interpolation[i])[[1]]
    
    # Read environmental data
    env <- read.csv("Data/environmental_data.csv", sep = ",")
    dat <- merge(dat, env, by = c("Reef", "Zone"), all = TRUE)
    dat <- dat[!is.na(dat$carrying), ]
    dat$carrying_prop <- dat$carrying / 100
    
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
    
    
    # models' explanatory variables:
    # m1 = hard substrate + temperature + Secchi depth
    # m2 = temperature + Secchi depth
    # m3 = hard substrate + Secchi depth
    # m4 = hard substrate + temperature 
    # m5 = hard substrate 
    # m6 = Secchi depth
    # m7 = temperature
    
    if (coef_df$model[i] == "m1") {
      prior <- get_prior(carrying_prop ~ st_habitat + 
                           st_temp_median +
                           st_Secchi +
                           (1| Reef), 
                         data = dat,
                         family = "beta")
      
      m <- brm(carrying_prop ~ st_habitat + 
                 st_temp_median +
                 st_Secchi + 
                 (1| Reef), 
               data = dat,
               family = "beta", chains = 3, cores = 3, 
               prior = prior, iter = 20000, thin = 5,
               control = list(adapt_delta = 0.99, max_treedepth = 17),
               save_pars = save_pars(all = TRUE), file = file_name)
    }
    
    if (coef_df$model[i] == "m2") {
      prior <- get_prior(carrying_prop ~ 
                           st_temp_median +
                           st_Secchi +
                           (1| Reef), 
                         data = dat,
                         family = "beta")
      
      m <- brm(carrying_prop ~  
                 st_temp_median +
                 st_Secchi + 
                 (1| Reef), 
               data = dat,
               family = "beta", chains = 3, cores = 3, 
               prior = prior, iter = 20000, thin = 5,
               control = list(adapt_delta = 0.99, max_treedepth = 17),
               save_pars = save_pars(all = TRUE), file = file_name)
      
    }
    if (coef_df$model[i] == "m3") {
      prior <- get_prior(carrying_prop ~ st_habitat + 
                           st_Secchi +
                           (1| Reef), 
                         data = dat,
                         family = "beta")
      
      m <- brm(carrying_prop ~ st_habitat + 
                 st_Secchi + 
                 (1| Reef), 
               data = dat,
               family = "beta", chains = 3, cores = 3, 
               prior = prior, iter = 20000, thin = 5,
               control = list(adapt_delta = 0.99, max_treedepth = 17),
               save_pars = save_pars(all = TRUE), file = file_name )
      
    }
    if (coef_df$model[i] == "m4") {
      prior <- get_prior(carrying_prop ~ st_habitat + 
                           st_temp_median +
                           (1| Reef), 
                         data = dat,
                         family = "beta")
      
      m <- brm(carrying_prop ~ st_habitat + 
                 st_temp_median +
                 (1| Reef), 
               data = dat,
               family = "beta", chains = 3, cores = 3, 
               prior = prior, iter = 20000, thin = 5,
               control = list(adapt_delta = 0.99, max_treedepth = 17),
               save_pars = save_pars(all = TRUE), file = file_name)
    }
    if (coef_df$model[i] == "m5") {
      prior <- get_prior(carrying_prop ~ st_habitat + 
                           (1| Reef), 
                         data = dat,
                         family = "beta")
      
      m <- brm(carrying_prop ~ st_habitat + 
                 (1| Reef), 
               data = dat,
               family = "beta", chains = 3, cores = 3, 
               prior = prior, iter = 20000, thin = 5,
               control = list(adapt_delta = 0.99, max_treedepth = 17),
               save_pars = save_pars(all = TRUE), file = file_name)
    }
    
    if (coef_df$model[i] == "m6") {
      prior <- get_prior(carrying_prop ~ 
                           st_Secchi +
                           (1| Reef), 
                         data = dat,
                         family = "beta")
      
      m <- brm(carrying_prop ~ 
                 st_Secchi + 
                 (1| Reef), 
               data = dat,
               family = "beta", chains = 3, cores = 3, 
               prior = prior, iter = 20000, thin = 5,
               control = list(adapt_delta = 0.99, max_treedepth = 17),
               save_pars = save_pars(all = TRUE), file = file_name)
      
    }
    if (coef_df$model[i] == "m7") {
      prior <- get_prior(carrying_prop ~ 
                           st_temp_median +
                           (1| Reef), 
                         data = dat,
                         family = "beta")
      
      m <- brm(carrying_prop ~ 
                 st_temp_median +
                 (1| Reef), 
               data = dat,
               family = "beta", chains = 3, cores = 3, 
               prior = prior, iter = 20000, thin = 5,
               control = list(adapt_delta = 0.99, max_treedepth = 17),
               save_pars = save_pars(all = TRUE), file = file_name)
      
    }
  }
  
  
  
  # Extract coefficient estimates from the fitted models
  post <- as.data.frame(as_draws_df(m)[ , 2:4])
  
  if(length(which(grepl("temp", colnames(post))) > 0)) {
    coef_df[i, ]$temperature = quantile(post[, which(grepl("temp", colnames(post)))], 0.5)
    coef_df[i, ]$temp_ci = ifelse(quantile(post[, which(grepl("temp", colnames(post)))], 0.025) < 0 &
                                    quantile(post[, which(grepl("temp", colnames(post)))], 0.975) > 0 ,"n", "y")
  }
  
  if(length(which(grepl("habitat", colnames(post))) > 0)) {
    coef_df[i, ]$habitat = quantile(post[, which(grepl("habitat", colnames(post)))], 0.5)
    coef_df[i, ]$habitat_ci = ifelse(quantile(post[, which(grepl("habitat", colnames(post)))], 0.025) < 0 &
                                       quantile(post[, which(grepl("habitat", colnames(post)))], 0.975) > 0 ,"n", "y")
  }
  
  
  if(length(which(grepl("Secchi", colnames(post))) > 0)) {
    coef_df[i, ]$Secchi = quantile(post[, which(grepl("Secchi", colnames(post)))], 0.5)
    coef_df[i, ]$Secchi_ci = ifelse(quantile(post[, which(grepl("Secchi", colnames(post)))], 0.025) < 0 &
                                      quantile(post[, which(grepl("Secchi", colnames(post)))], 0.975) > 0 ,"n", "y")
  }
  
  coef_df[i, ]$n_sites = nrow(m$data)
  coef_df[i, ]$n_reefs = length(unique(m$data$Reef))
  print(i/nrow(coef_df))
  
}


# Add information on the definition of upper limits used (RM = running mean window, 
# threshold = thresholds for the difference between the running mean and the maximum coral cover)
coef_df$type = paste(paste("RM", coef_df$running, sep = " "), 
                     paste("threshold ", coef_df$percentage, "%", sep = ""), sep = "-")



coef_df = coef_df[complete.cases(coef_df$model), ]


# Plot the coefficient estimates for the three environmental variables obtained
# across all models.

temperature_plot = ggplot(coef_df[!coef_df$interpolation == FALSE, ], aes(x = type,
                                                                          y = model, fill = temperature)) +
  geom_tile(col = "black") +
  geom_text(data =coef_df[!coef_df$interpolation == FALSE& coef_df$temp_ci == "y", ],
            aes(x = type, y = model, label = "*")
  ) +
  scale_fill_gradient2(high = "blue", low = "red", mid = "white", midpoint = 0, limits = c(-1, 1),
                       na.value = "#CCCCCC", name = "") +
  scale_y_discrete(limits = c("m7", "m6", "m5", "m4", "m3", "m2", "m1"),
                   labels = c(rep("", 7))) +
  scale_x_discrete(limits = c("RM 5-threshold 80%" ,
                              "RM 5-threshold 85%", "RM 5-threshold 90%",
                              "RM 7-threshold 80%" ,
                              "RM 7-threshold 85%", "RM 7-threshold 90%",
                              "RM 9-threshold 80%" ,
                              "RM 9-threshold 85%", "RM 9-threshold 90%",
                              "RM 11-threshold 80%" ,
                              "RM 11-threshold 85%", "RM 11-threshold 90%"),
                   labels = c( "5-80%", "5-85%", "5-90%", "7-80%",
                              "7-85%", "7-90%", "9-80%", "9-85%", "9-90%",
                              "11-80%", "11-85%", "11-90%")) +
  xlab("") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 15)) +
  ylab("") +
  ggtitle("e - temperature")

habitat_plot = ggplot(coef_df[!coef_df$interpolation == FALSE, ], aes(x = type,
                                                                      y = model, fill = habitat)) +
  geom_tile(col = "black") +
  geom_text(data =coef_df[!coef_df$interpolation == FALSE& coef_df$habitat_ci == "y", ],
            aes(x = type, y = model, label = "*")
  ) +
  scale_fill_gradient2(high = "blue", low = "red", mid = "white", midpoint = 0, limits = c(-1, 1),
                       na.value = "#CCCCCC", name = "") +
  scale_y_discrete(limits = c("m7", "m6", "m5", "m4", "m3", "m2", "m1"),
                   labels = c("temperature", "Secchi", "hard substrate", 
                              "hard substrate + temperature", 
                              "hard substrate + Secchi", "temperature + Secchi",
                              "hard substrate + temperature + Secchi")) +
  scale_x_discrete(limits = c( "RM 5-threshold 80%" ,
                              "RM 5-threshold 85%", "RM 5-threshold 90%",
                              "RM 7-threshold 80%" ,
                              "RM 7-threshold 85%", "RM 7-threshold 90%",
                              "RM 9-threshold 80%" ,
                              "RM 9-threshold 85%", "RM 9-threshold 90%",
                              "RM 11-threshold 80%" ,
                              "RM 11-threshold 85%", "RM 11-threshold 90%"),
                   labels = c( "5-80%", "5-85%", "5-90%", "7-80%",
                              "7-85%", "7-90%", "9-80%", "9-85%", "9-90%",
                              "11-80%", "11-85%", "11-90%")) +
  xlab("") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 15)) +
  ylab("") +
  ggtitle("d - hard substrate availability")

Secchi_plot =
  ggplot(coef_df[!coef_df$interpolation == FALSE, ], aes(x = type,
                                                         y = model, fill = Secchi)) +
  geom_tile(col = "black") +
  geom_text(data =coef_df[!coef_df$interpolation == FALSE& coef_df$Secchi_ci == "y", ],
            aes(x = type, y = model, label = "*")
  ) +
  scale_fill_gradient2(high = "blue", low = "red", mid = "white", midpoint = 0, limits = c(-1, 1),
                       na.value = "#CCCCCC", name = "") +
  scale_y_discrete(limits = c("m7", "m6", "m5", "m4", "m3", "m2", "m1"),
                   labels = c(rep("", 7))) +
  scale_x_discrete(limits = c( "RM 5-threshold 80%" ,
                              "RM 5-threshold 85%", "RM 5-threshold 90%",
                              "RM 7-threshold 80%" ,
                              "RM 7-threshold 85%", "RM 7-threshold 90%",
                              "RM 9-threshold 80%" ,
                              "RM 9-threshold 85%", "RM 9-threshold 90%",
                              "RM 11-threshold 80%" ,
                              "RM 11-threshold 85%", "RM 11-threshold 90%"),
                   labels = c( "5-80%", "5-85%", "5-90%", "7-80%",
                              "7-85%", "7-90%", "9-80%", "9-85%", "9-90%",
                              "11-80%", "11-85%", "11-90%")) +
  xlab("") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 15)) +
  ylab("") +
  ggtitle("f - Secchi depth")

habitat_plot + temperature_plot + Secchi_plot


# Plot how environmental variables are distributed across space
world_map <- map_data("world")
Aus <- world_map %>% filter(region == "Australia")

sum_dat1 = get_data_quantile(manta)$dat
sum_dat1 = as.data.frame(sum_dat1)

mean_vals = get_data_quantile(manta)$mean_vals

habitat_map <- 
  ggplot() +
  geom_polygon(data = Aus, aes(x = long, y = lat, group = group), fill = "grey",
               col = "black") +
  coord_map(xlim = c(142, 153.5), ylim = c(-25, -9)) +
  geom_point(data = sum_dat1, aes(x = Longitude, y = Latitude,
                                  fill = per_suitable *100),
             pch = 21, col = "black", size = 3,  alpha = 0.5) +
  scale_fill_gradient2(high = "blue", low = "red", mid = "white", midpoint = mean_vals[mean_vals$variable == "habitat", ]$mean *100,
                       name = "%") +
  theme_bw() +
  ylab("") +
  xlab("") +
  scale_y_continuous(breaks = seq(-25, -10, by = 5),
                     labels = c("25°S", "20°S", "15°S", "10°S")) +
  scale_x_continuous(breaks = seq(144, 153, by = 3),
                     labels = c("144°E", "147°E", "150°E", "153°E")) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  ggtitle("a- hard substrate availability")


temperature_map <- 
  ggplot() +
  geom_polygon(data = Aus, aes(x = long, y = lat, group = group), fill = "grey",
               col = "black") +
  coord_map(xlim = c(142, 153.5), ylim = c(-25, -9)) +
  geom_point(data = sum_dat1, aes(x = Longitude, y = Latitude,
                                  fill = temp_median),
             pch = 21, col = "black", size = 3,  alpha = 0.5) +
  scale_fill_gradient2(high = "blue", low = "red", mid = "white", midpoint = mean_vals[mean_vals$variable == "temperature", ]$mean,
                       name = "°C") +
  theme_bw() +
  ylab("") +
  xlab("") +
  scale_y_continuous(breaks = seq(-25, -10, by = 5),
                     labels = c("25°S", "20°S", "15°S", "10°S")) +
  scale_x_continuous(breaks = seq(144, 153, by = 3),
                     labels = c("144°E", "147°E", "150°E", "153°E")) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  ggtitle("b- temperature")



Secchi_map <- 
  ggplot() +
  geom_polygon(data = Aus, aes(x = long, y = lat, group = group), fill = "grey",
               col = "black") +
  coord_map(xlim = c(142, 153.5), ylim = c(-25, -9)) +
  geom_point(data = sum_dat1, aes(x = Longitude, y = Latitude,
                                  fill = Secchi_median),
             pch = 21, col = "black", size = 3,  alpha = 0.5) +
  scale_fill_gradient2(high = "blue", low = "red", mid = "white", midpoint = mean_vals[mean_vals$variable == "Secchi", ]$mean,
                       name = "m") +
  theme_bw() +
  ylab("") +
  xlab("") +
  scale_y_continuous(breaks = seq(-25, -10, by = 5),
                     labels = c("25°S", "20°S", "15°S", "10°S")) +
  scale_x_continuous(breaks = seq(144, 153, by = 3),
                     labels = c("144°E", "147°E", "150°E", "153°E")) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  ggtitle("c- Secchi depth")


Figure3 = ggarrange(habitat_map + ggtitle("a") , temperature_map + ggtitle("b") , 
                    Secchi_map + ggtitle("c"), habitat_plot+ ggtitle("d"), 
                    temperature_plot + ggtitle("e"), Secchi_plot + ggtitle("f"), ncol = 3)

#ggsave(paste(path_out, "Figure3.png", sep = "/"), Figure4, width = 37, height = 45, units = "cm", dpi = 300)

################################################################################
#  How do hard substrate availability, temperature, and  water clarity relate  #
#  to upper limits - sensitivity analyses (end)                              #
################################################################################

################################################################################
#           Explore within reef variation in upper limits (start)              #
################################################################################

# Use the upper limit with a 5-year running window and a threshold of 80% as
# an example
dat = get_stable_cover(manta, 5, 0.8,  
                       TRUE)[[1]]

# Read environmental data
env <- read.csv("Data/environmental_data.csv", sep = ",")
dat <- merge(dat, env, by = c("Reef", "Zone"), all = TRUE)
dat <- dat[!is.na(dat$carrying), ]
dat$carrying_prop <- dat$carrying / 100

# save a copy of the data frame that includes all upper limit estimates
dat1 <- dat

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


# Check how upper limits vary within reef
lat_reefs = dat1[order(dat1$Latitude, decreasing = TRUE), ]$Reef
lat_reefs = lat_reefs[!duplicated(lat_reefs)]

FigureS3a = ggplot(dat1, aes(x = Reef, y = carrying)) +
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


FigureS3b = ggplot(range_sum, aes(x = n, y = Range)) +
  geom_boxplot(aes(group = n)) +
  theme_bw() +
  xlab("Number of sites with upper limit estimates per reef") +
  ylab("Range of upper limit \n estimates within reef (%)") +
  ggtitle("b")


FigureS3 <- FigureS3a / FigureS3b

################################################################################
#           Explore within reef variation in upper limits (end)                #
################################################################################


################################################################################
#                  Spatial patterns in upper limits (start)                    #
################################################################################


# Plot spatial patterns in upper limits
dat$Reef_zone <- paste(dat$Reef, dat$Zone, sep = "-")

dat$Zone <- ifelse(dat$Zone == "Flank1", "Flank",
                   ifelse(dat$Zone == "Flank2", "Flank",
                          dat$Zone))

# Summarise upper limits by reef
sum_dat1 = dat1 %>% group_by(Reef) %>%
  summarise("Latitude" = mean(Latitude), "Longitude" = mean(Longitude),
            "carrying_prop" = mean(carrying_prop), "n" = n())
sum_dat1 = as.data.frame(sum_dat1)

sizes = c(3,4,5,6)
names(sizes) = c("1", "2", "3", "4")

# Plot variation in upper limits across the Great Barrier Reef
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

# Plot variation in upper limits across a wave exposure gradient

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
#ggsave(paste(path_out, "Figure2.png", sep = "/"), Figure2, width = 25, height = 20, units = "cm", dpi = 300)

################################################################################
#                  Spatial patterns in upper limits (end)                      #
################################################################################


################################################################################
#                        Quantile regression (start)                           #
################################################################################

# Get yearly estimates of coral cover
year_cover = get_year_cover(manta)$dat
mean_vals_q = get_year_cover(manta)$mean_vals

# Fit  quantile regressions to yearly coral cover estimates
quant = fit_quant_reg(year_cover)

# Extract and plot its conditional effects for each quantile (95, 50 and 5)
cond_95 = conditional_effects(quant$m1, re_formula = NA, robust = TRUE, dpar = "mu")

cond_95$st_habitat$habitat = convert_scales("habitat", cond_95$st_habitat$st_habitat,
                                            mean_vals_q)
cond_95$st_temp_median$temperature = convert_scales("temperature", cond_95$st_temp_median$st_temp_median,
                                                    mean_vals_q)
cond_95$st_Secchi$Secchi = convert_scales("Secchi", cond_95$st_Secchi$st_Secchi,
                                          mean_vals_q)
cond_50 = conditional_effects(quant$m2, re_formula = NA, robust = TRUE, dpar = "mu")

cond_50$st_habitat$habitat = convert_scales("habitat", cond_50$st_habitat$st_habitat,
                                            mean_vals_q)
cond_50$st_temp_median$temperature = convert_scales("temperature", cond_50$st_temp_median$st_temp_median,
                                                    mean_vals_q)
cond_50$st_Secchi$Secchi = convert_scales("Secchi", cond_50$st_Secchi$st_Secchi,
                                          mean_vals_q)

cond_05 = conditional_effects(quant$m3, re_formula = NA, robust = TRUE, dpar = "mu")

cond_05$st_habitat$habitat = convert_scales("habitat", cond_05$st_habitat$st_habitat,
                                            mean_vals_q)
cond_05$st_temp_median$temperature = convert_scales("temperature", cond_05$st_temp_median$st_temp_median,
                                                    mean_vals_q)
cond_05$st_Secchi$Secchi = convert_scales("Secchi", cond_05$st_Secchi$st_Secchi,
                                          mean_vals_q)
Figure4b = ggplot() +
  geom_point(data = year_cover, aes(x = per_suitable * 100, y = cover), alpha = 0.05) +
  geom_ribbon(data = cond_95$st_habitat, aes(x = habitat * 100, ymin = lower__, ymax = upper__), fill = "#000033", alpha = 0.6) +
  geom_ribbon(data = cond_50$st_habitat, aes(x = habitat * 100, ymin = lower__, ymax = upper__), fill = "#0000FF", alpha = 0.6) +
  geom_ribbon(data = cond_05$st_habitat, aes(x = habitat * 100, ymin = lower__, ymax = upper__), fill = "#99CCFF", alpha = 0.6) +
  geom_line(data = cond_95$st_habitat, aes(x = habitat * 100, y = estimate__), col = "#000033", size = 1.2) +
  geom_line(data = cond_50$st_habitat, aes(x = habitat * 100, y = estimate__), col = "#0000FF", size = 1.2) +
  geom_line(data = cond_05$st_habitat, aes(x = habitat * 100, y = estimate__), col = "#99CCFF", size = 1.2) +
  scale_y_continuous(limits = c(0, 100)) +
  xlab("Hard substrate availability (%)") +
  ylab("Coral cover (%)") +
  theme_bw() +
  theme(text = element_text(size = 15))

Figure4c =ggplot() +
  geom_point(data = year_cover, aes(x = temp_median, y = cover), alpha = 0.05) +
  geom_ribbon(data = cond_95$st_temp_median, aes(x = temperature, ymin = lower__, ymax = upper__), fill = "#000033", alpha = 0.6) +
  geom_ribbon(data = cond_50$st_temp_median, aes(x = temperature, ymin = lower__, ymax = upper__), fill = "#0000FF", alpha = 0.6) +
  geom_ribbon(data = cond_05$st_temp_median, aes(x = temperature, ymin = lower__, ymax = upper__), fill = "#99CCFF", alpha = 0.6) +
  geom_line(data = cond_95$st_temp_median, aes(x = temperature, y = estimate__), col = "#000033", size = 1.2) +
  geom_line(data = cond_50$st_temp_median, aes(x = temperature, y = estimate__), col = "#0000FF", size = 1.2) +
  geom_line(data = cond_05$st_temp_median, aes(x = temperature, y = estimate__), col = "#99CCFF", size = 1.2) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 100)) +
  xlab("Temperature (°C)") +
  ylab("Coral cover (%)") +
  theme_bw() +
  theme(text = element_text(size = 15))

Figure4d = ggplot() +
  geom_point(data = year_cover, aes(x = Secchi_median, y = cover), alpha = 0.05) +
  geom_ribbon(data = cond_95$st_Secchi, aes(x = Secchi, ymin = lower__, ymax = upper__), fill = "#000033", alpha = 0.6) +
  geom_ribbon(data = cond_50$st_Secchi, aes(x = Secchi, ymin = lower__, ymax = upper__), fill = "#0000FF", alpha = 0.6) +
  geom_ribbon(data = cond_05$st_Secchi, aes(x = Secchi, ymin = lower__, ymax = upper__), fill = "#99CCFF", alpha = 0.6) +
  geom_line(data = cond_95$st_Secchi, aes(x = Secchi, y = estimate__), col = "#000033", size = 1.2) +
  geom_line(data = cond_50$st_Secchi, aes(x = Secchi, y = estimate__), col = "#0000FF", size = 1.2) +
  geom_line(data = cond_05$st_Secchi, aes(x = Secchi, y = estimate__), col = "#99CCFF", size = 1.2) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 100)) +
  xlab("Secchi depth (m)") +
  ylab("Coral cover (%)") +
  theme_bw() +
  theme(text = element_text(size = 15))


# Extract coefficient estimates
post_95 <- as.data.frame(as_draws_df(quant$m1)[ , 2:4])
post_50 <- as.data.frame(as_draws_df(quant$m2)[ , 2:4])
post_10 <- as.data.frame(as_draws_df(quant$m3)[ , 2:4])
colnames(post_95) <- c( "habitat", "temp", "Secchi")
colnames(post_50) <- c( "habitat", "temp", "Secchi")
colnames(post_10) <- c( "habitat", "temp", "Secchi")
post_95$quant = "0.95"
post_50$quant = "0.50"
post_10$quant = "0.10"

post_quant = data.frame("variable" = c(rep("habitat", nrow(post_95)),
                                       rep("temp", nrow(post_95)),
                                       rep("Secchi", nrow(post_95))),
                        "value" = c(post_95$habitat, post_95$temp, post_95$Secchi),
                        "quantile" = "0.95")
post_quant = rbind(post_quant,
                   data.frame("variable" = c(rep("habitat", nrow(post_50)),
                                             rep("temp", nrow(post_50)),
                                             rep("Secchi", nrow(post_50))),
                              "value" = c(post_50$habitat, post_50$temp, post_50$Secchi),
                              "quantile" = "0.50") )
post_quant = rbind(post_quant,
                   data.frame("variable" = c(rep("habitat", nrow(post_10)),
                                             rep("temp", nrow(post_10)),
                                             rep("Secchi", nrow(post_10))),
                              "value" = c(post_10$habitat, post_10$temp, post_10$Secchi),
                              "quantile" = "0.10") )

Figure4a =ggplot(post_quant) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_density_ridges(aes(x = value, y = variable, fill = quantile), rel_min_height = 0.00001,
                      alpha = 0.6) +
  scale_y_discrete(limits = c( "Secchi", "temp", "habitat"), 
                   labels = c("Secchi depth", "Temperature", "Hard substrate" )) +
  ylab("") +
  xlab("Posterior distributions of coefficient estimates") + 
  scale_fill_manual(values = c("#99CCFF", "#0000FF", "#000033"), name = "Quantile") +
  theme_bw() +
  theme(text = element_text(size = 15))

Figure4 = (Figure4a + ggtitle("a")) /((Figure4b + ggtitle("b")) + (Figure4c + ggtitle("c") + ylab("")) +
                                        (Figure4d + ggtitle("d") + ylab("")))
ggsave(paste(path_out, "Figure4.png", sep = "/"), Figure4, width = 33, height = 20, units = "cm", dpi = 300)



# Generate predictions on how upper limits vary with changes in the environmental
# variables (for the text in the results section). 

# Change in upper limits with a 10% change in hard substrate availability
preds_habitat = data.frame("habitat" = c(0.9, 1))
preds_habitat$st_temp_median = cond_05$st_habitat$st_temp_median[1]
preds_habitat$st_Secchi = cond_05$st_habitat$st_Secchi[1]
preds_habitat$st_habitat = (preds_habitat$habitat - mean_vals_q[mean_vals_q$variable == "habitat", ]$mean)/
  mean_vals_q[mean_vals_q$variable == "habitat", ]$sd
cbind(preds_habitat, as.data.frame(fitted(quant$m3, preds_habitat, re_formula = NA, type = "response")))

# Change in upper limits with a 1 degree C change in temperature
preds_temp = data.frame("temp" = c(24, 25))
preds_temp$st_habitat = cond_05$st_temp_median$st_habitat[1]
preds_temp$st_Secchi = cond_05$st_temp_median$st_Secchi[1]
preds_temp$st_temp_median = (preds_temp$temp - mean_vals_q[mean_vals_q$variable == "temperature", ]$mean)/
  mean_vals_q[mean_vals_q$variable == "temperature", ]$sd
cbind(preds_temp, as.data.frame(fitted(quant$m3, preds_temp, re_formula = NA, type = "response")))

# Change in upper limits with a 1m  change in water clarity
preds_Secchi = data.frame("Secchi" = c(13, 14))
preds_Secchi$st_temp_median = cond_05$st_Secchi$st_temp_median[1]
preds_Secchi$st_habitat = cond_05$st_Secchi$st_habitat[1]
preds_Secchi$st_Secchi = (preds_Secchi$Secchi - mean_vals[mean_vals$variable == "Secchi", ]$mean)/
  mean_vals[mean_vals$variable == "Secchi", ]$sd
cbind(preds_Secchi, as.data.frame(fitted(quant$m3, preds_Secchi, re_formula = NA, type = "response")))


################################################################################
#                        Quantile regression   (end)                           #
################################################################################



################################################################################
#       Explore recovery patterns with different upper limits (start)          #
################################################################################


logistic_eq <- function(K, N, R) {
  # Logistic population growth equation in discrete time
  # Nt1 is coral cover at time t+1
  # N is coral cover at time t
  # R is the maximum rate of growth
  # K is the upper limit
  
  Nt1 <- N + R * (N * (K - N) / K)
  return(list("N" = Nt1))
  
}


# Values to explore
carrying_prop <- seq(0.01, 1, by = 0.01) # upper limits
coral_cover <- seq(0.02, 1, by = 0.01)   # coral cover target
R <- c(0.25, 0.5)   # maximum rate of growth


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



# Plot recovery patterns when R = 0.25
Figure5a <- ggplot()+
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


# Plot recovery patterns when R = 0.5
Figure5b <- ggplot()+
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


Figure5 <- (Figure5a + xlab("")) / (Figure5b + theme(legend.position = "none"))
#ggsave(paste(path_out, "Figure5.png", sep = "/"), Figure5, width = 20, height = 25, units = "cm", dpi = 300)


################################################################################
#         Explore recovery patterns with different upper limits (end)          #
################################################################################

