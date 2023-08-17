rm(list = ls())


set.seed(9)

library(ggplot2)
library(brms)
library(dplyr)
library(viridis)
library(patchwork)
library(maps)


dat <- read.csv("realised_upper_limits.csv", sep = ",")
env <- read.csv("environmental_variables.csv", sep = ",")
dat <- merge(dat, env, by = c("Reef", "Zone"), all = TRUE)
dat <- dat[!is.na(dat$carrying), ]

# Define realised upper limit as a proportion
dat$carrying_prop <- dat$carrying / 100

# Save all data
dat1 <- dat



dat <- dat[complete.cases(dat[ , c("median_ubed90", 
                                   "temp_median",
                                   "PAR8_median", 
                                   "per_suitable",
                                   "Zone", "Latitude", "Longitude") ]), ]

# Exclude data with less than 80% hard substrate availability
dat <- dat[dat$per_suitable >= 0.8, ]

# Convert environmental variables to Z-scores
dat$st_ubed90 <- (dat$median_ubed90 - mean(dat$median_ubed90, na.rm = TRUE)) /
  sd(dat$median_ubed90, na.rm = TRUE)
dat$st_temp_median <- (dat$temp_median - mean(dat$temp_median, na.rm = TRUE)) /
  sd(dat$temp_median, na.rm = TRUE)
dat$st_PAR8 <- (dat$PAR8_median  - mean(dat$PAR8_median, na.rm = TRUE)) /
  (sd(dat$PAR8_median, na.rm = TRUE))


# Get mean and standard deviation of each variable to convert model outputs if
# needed.
mean_vals <- data.frame("variable" = c("temperature", "light", "benthic stress"),
                        "mean" = c(mean(dat$temp_median),
                                   mean(dat$PAR8_median),
                                   mean(dat$median_ubed90)),
                        "sd" = c(sd(dat$temp_median),
                                 sd(dat$PAR8_median),
                                 sd(dat$median_ubed90)))


# Convert latitudes and longitudes into kilometer-based grid coordinates
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

# Create covariance matrix from distances between reefs
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




# Plot covariance matrix
melted <- melt(spat)
colnames(melted)[1:2] <- c("X1", "X2")
ggplot(melted) +
  geom_tile(aes(x = X1, y = X2, fill = value)) +
  scale_fill_viridis(option = "B") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") +
  ylab("")


dat$Reef_spat <- dat$Reef


dat <- dat[complete.cases(dat[ , c("st_ubed90", 
                                   "st_temp_median",
                                   "st_PAR8", 
                                   "Zone", "x", "y", "Latitude", "Longitude") ]), ]
dat$Zone <- ifelse(dat$Zone == "Flank1", "Flank",
                   ifelse(dat$Zone == "Flank2", "Flank",
                          dat$Zone))




# Plot broad and fine-scale spatial patterns in coral cover upper limits

world_map <- map_data("world")
Aus <- world_map %>% filter(region == "Australia")


Figure2a <- 
  ggplot() +
  geom_polygon(data = Aus, aes(x = long, y = lat, group = group), fill = "grey",
               col = "black") +
  coord_map(xlim = c(142, 153.5), ylim = c(-25, -9)) +
  geom_point(data = dat1, aes(x = Longitude, y = Latitude,
                              fill =  carrying_prop * 100),
             pch = 21, col = "black", alpha = 0.5, size = 3) +
  scale_fill_viridis(limits = c(0, 80),
                     breaks = seq(0, 100, by = 20),
                     labels = paste0(seq(0, 100, by = 20), "%"),
                     option = "B",
                     name = "") +
  theme_bw() +
  ylab("") +
  xlab("") +
  scale_y_continuous(breaks = seq(-25, -10, by = 5),
                     labels = c("25° S", "20° S", "15° S", "10° S")) +
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
  xlab("") +
  ylab("Realised upper limit (%)") +
  theme_bw() +
  theme(text = element_text(size = 14 ), plot.title = element_text(size = 16)) +
  ggtitle("b - Within reef variation")


Figure2 <- Figure2a + Figure2b

#ggsave("Figure2.tiff", Figure2, dpi = 300)



################################################################################
####         Realised upper limits and environmental variables             #####
################################################################################


# Fit models three-way-interaction models, with and without correcting for
# spatial autocorrelation and compare

prior <- get_prior(carrying_prop ~  st_temp_median * st_PAR8 * 
                     st_ubed90 +  (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat),
                   family = "beta")

m1a <- brm(carrying_prop ~  st_temp_median * st_PAR8 *
             st_ubed90 +  (1| gr(Reef_spat, cov = spat)),
           data = dat,
           data2 = list(spat = spat),
           family = "beta", chains = 3, cores = 3, 
           prior = prior, iter = 10000, thin = 5,
           control = list(adapt_delta = 0.99, max_treedepth = 17),
           save_pars = save_pars(all = TRUE), file = "m1a")
m1a <- add_criterion(m1a, "loo", moment.match = TRUE)

prior <- get_prior(carrying_prop ~  st_temp_median * st_PAR8 * 
                     st_ubed90 + (1 | Reef) ,
                   data = dat,
                   family = "beta")

m1a_2 <- brm(carrying_prop ~  st_temp_median * st_PAR8 *
               st_ubed90 + (1 | Reef) ,
             data = dat,
             family = "beta", chains = 3, cores = 3, 
             prior = prior, iter = 10000, thin = 5,
             control = list(adapt_delta = 0.99, max_treedepth = 17),
             save_pars = save_pars(all = TRUE), file = "m1a_2")
m1a_2 <- add_criterion(m1a_2, "loo", moment.match = TRUE)

loo_compare(m1a, m1a_2)


# Fit models with all combinations of interactions, correcting for 
# spatial autocorrelation
prior <- get_prior(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
                     st_temp_median:st_PAR8 + st_temp_median:st_ubed90 +
                     st_PAR8:st_ubed90 + (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat),
                   family = "beta")

m1b <- brm(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
             st_temp_median:st_PAR8 + st_temp_median:st_ubed90 +
             st_PAR8:st_ubed90 + (1| gr(Reef_spat, cov = spat)),
           data = dat,
           data2 = list(spat = spat),
           family = "beta", chains = 3, cores = 3, 
           prior = prior, iter = 10000, thin = 5,
           control = list(adapt_delta = 0.99, max_treedepth = 17),
           save_pars = save_pars(all = TRUE), file = "m1b")

m1b <- add_criterion(m1b, "loo", moment.match = TRUE)



prior <- get_prior(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
                     st_temp_median:st_PAR8 + st_temp_median:st_ubed90 +
                     (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat),
                   family = "beta")

m1c <- brm(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
             st_temp_median:st_PAR8 + st_temp_median:st_ubed90 +
              (1| gr(Reef_spat, cov = spat)),
           data = dat,
           data2 = list(spat = spat),
           family = "beta", chains = 3, cores = 3, 
           prior = prior, iter = 10000, thin = 5,
           control = list(adapt_delta = 0.99, max_treedepth = 17),
           save_pars = save_pars(all = TRUE), file = "m1c")

m1c <- add_criterion(m1c, "loo", moment.match = TRUE)

prior <- get_prior(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
                     st_temp_median:st_PAR8 + 
                     st_PAR8:st_ubed90 + (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat),
                   family = "beta")

m1d <- brm(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
             st_temp_median:st_PAR8 + 
             st_PAR8:st_ubed90 + (1| gr(Reef_spat, cov = spat)),
           data = dat,
           data2 = list(spat = spat),
           family = "beta", chains = 3, cores = 3, 
           prior = prior, iter = 10000, thin = 5,
           control = list(adapt_delta = 0.99, max_treedepth = 17),
           save_pars = save_pars(all = TRUE), file = "m1d")

m1d <- add_criterion(m1d, "loo", moment.match = TRUE)



prior <- get_prior(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
                     st_temp_median:st_ubed90 +
                     st_PAR8:st_ubed90 + (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat),
                   family = "beta")

m1e <- brm(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
             st_temp_median:st_ubed90 +
             st_PAR8:st_ubed90 + (1| gr(Reef_spat, cov = spat)),
           data = dat,
           data2 = list(spat = spat),
           family = "beta", chains = 3, cores = 3, 
           prior = prior, iter = 10000, thin = 5,
           control = list(adapt_delta = 0.99, max_treedepth = 17),
           save_pars = save_pars(all = TRUE), file = "m1e")

m1e <- add_criterion(m1e, "loo", moment.match = TRUE)



prior <- get_prior(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
                     st_temp_median:st_PAR8 + (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat),
                   family = "beta")

m1f <- brm(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
             st_temp_median:st_PAR8 + (1| gr(Reef_spat, cov = spat)),
           data = dat,
           data2 = list(spat = spat),
           family = "beta", chains = 3, cores = 3, 
           prior = prior, iter = 10000, thin = 5,
           control = list(adapt_delta = 0.99, max_treedepth = 17),
           save_pars = save_pars(all = TRUE), file = "m1f")

m1f <- add_criterion(m1f, "loo", moment.match = TRUE)


prior <- get_prior(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
                     st_temp_median:st_ubed90 +
                     (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat),
                   family = "beta")

m1g <- brm(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
             st_temp_median:st_ubed90 +
             (1| gr(Reef_spat, cov = spat)),
           data = dat,
           data2 = list(spat = spat),
           family = "beta", chains = 3, cores = 3, 
           prior = prior, iter = 10000, thin = 5,
           control = list(adapt_delta = 0.99, max_treedepth = 17),
           save_pars = save_pars(all = TRUE), file = "m1g")

m1g <- add_criterion(m1g, "loo", moment.match = TRUE)





prior <- get_prior(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
                     st_PAR8:st_ubed90 + (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat),
                   family = "beta")

m1h <- brm(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
             st_PAR8:st_ubed90 + (1| gr(Reef_spat, cov = spat)),
           data = dat,
           data2 = list(spat = spat),
           family = "beta", chains = 3, cores = 3, 
           prior = prior, iter = 10000, thin = 5,
           control = list(adapt_delta = 0.99, max_treedepth = 17),
           save_pars = save_pars(all = TRUE), file = "m1h")

m1h <- add_criterion(m1h, "loo", moment.match = TRUE)

prior <- get_prior(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
                     (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat),
                   family = "beta")


m1i <- brm(carrying_prop ~  st_temp_median + st_PAR8 + st_ubed90 +
             (1| gr(Reef_spat, cov = spat)),
           data = dat,
           data2 = list(spat = spat),
           family = "beta", chains = 3, cores = 3, 
           prior = prior, iter = 10000, thin = 5,
           control = list(adapt_delta = 0.99, max_treedepth = 17),
           save_pars = save_pars(all = TRUE), file = "m1i")

m1i <- add_criterion(m1i, "loo", moment.match = TRUE)

loo_compare(m1a, m1b, m1c, m1d, m1e, m1f, m1g, m1h, m1i)

# Get best-fit model's Rsquared
bayes_R2(m1i)


## Example 1 for results section
example_preds <- data.frame("temperature" = c(25, 26)
)
example_preds$st_temp_median <- (example_preds$temperature - mean_vals[mean_vals$variable == "temperature", ]$mean)/
  mean_vals[mean_vals$variable == "temperature", ]$sd
example_preds$st_ubed90 <- 0
example_preds$st_PAR8 <- 0

example_preds <- cbind(example_preds, as.data.frame(fitted(m1i, example_preds, 
                                                           re_formula = NA)))

## Example 2 for results section
example_preds <- data.frame("ubed90" = c(0.1, 0.5))

example_preds$st_temp_median <- 0
example_preds$st_PAR8 <- 0
example_preds$st_ubed90 <- (example_preds$ubed90 - mean_vals[mean_vals$variable == "benthic stress", ]$mean)/
  mean_vals[mean_vals$variable == "benthic stress", ]$sd

example_preds <- cbind(example_preds, as.data.frame(fitted(m1i, example_preds, 
                                                           re_formula = NA)))


# Plot best-fit model's predictions

  # Correct realised upper limit estimates for spatial autocorrelation
dat$rand_m1 <- as.data.frame(fitted(m1i, dat))$Estimate - 
  as.data.frame(fitted(m1i, dat, re_formula = NA))$Estimate  


temp <- expand.grid("st_temp_median" = seq(min(dat$st_temp_median), max(dat$st_temp_median), l = 100),
                    "st_PAR8" = 0,
                    "st_ubed90" = 0)

temp <- cbind(temp, as.data.frame(fitted(m1i, temp, re_formula = NA)))


Figure4b <- ggplot() +
  geom_point(data = dat, aes( x = st_temp_median, y = (carrying_prop - rand_m1) * 100), col = "black",
             pch = 21, fill = "grey", size = 2) +
  geom_ribbon(data = temp, aes( x = st_temp_median, ymin = Q2.5 * 100, ymax = Q97.5 * 100), alpha = 0.2) +
  geom_line(data = temp, aes( x = st_temp_median, y = Estimate * 100), size = 1.5) +
  xlab("Temperature (S.D.)") +
  ylab("Realised upper limit (%)") +
  theme_bw() +
  theme(text = element_text(size = 14 )) +
  ggtitle("b")
Figure4b 



ubed90 <- expand.grid(
  "st_temp_median" = 0,
  "st_PAR8" = 0,
  "st_ubed90" = seq(min(dat$st_ubed90), max(dat$st_ubed90), l = 100))

ubed90 <- cbind(ubed90, as.data.frame(fitted(m1i, ubed90, re_formula = NA)))



Figure4c <- ggplot() +
  geom_point(data = dat, aes( x = st_ubed90, y = (carrying_prop - rand_m1)* 100), fill = "grey", col = "black",
             pch = 21, size = 2) +
  geom_ribbon(data = ubed90, aes( x = st_ubed90, ymin = Q2.5 * 100, ymax = Q97.5 * 100), alpha = 0.2, fill = "black") +
  geom_line(data = ubed90, aes( x = st_ubed90, y = Estimate * 100), size = 1.1, col = "black") +
  xlab("Benthic stress (S.D)") +
  ylab("Realised upper limit (%)") +
  theme_bw() +
  theme(text = element_text(size = 14 )) +
  ggtitle("c")
Figure4c 


light <- expand.grid(
  "st_temp_median" = 0,
  "st_PAR8" = seq(min(dat$st_PAR8), max(dat$st_PAR8), l = 100),
  "st_ubed90" = 0)

light <- cbind(light, as.data.frame(fitted(m1i, light, re_formula = NA)))


Figure4d <- ggplot() +
  geom_point(data = dat, aes( x = st_PAR8, y = (carrying_prop - rand_m1) * 100), fill = "grey", col = "black",
             pch = 21, size = 2) +
  geom_ribbon(data = light, aes( x = st_PAR8, ymin = Q2.5 * 100, ymax = Q97.5 * 100), alpha = 0.2, fill = "black") +
  geom_line(data = light, aes( x = st_PAR8, y = Estimate * 100), size = 1.1, col = "black") +
  xlab("Light (S.D)") +
  ylab("Realised upper limit (%)") +
  theme_bw() +
  theme(text = element_text(size = 14 )) +
  ggtitle("d")
Figure4d 




post <- as.data.frame(as_draws_df(m1i)[ , 2:4])
colnames(post) <- c("temp", "PAR8", "ubed90")


Figure4a <- ggplot() +
  geom_violin(data = post, aes(x = 1, y = temp), fill = "grey" , alpha = 0.2) +
  geom_linerange(data = post, aes(x = 1, ymin = quantile(temp, 0.025), ymax = quantile(temp, 0.975)), 
                 size = 1.1) +
  geom_point(data = post, aes(x = 1, y = quantile(temp, 0.5)), 
             size = 5) +
  
  geom_violin(data = post, aes(x = 2, y = ubed90), fill = "grey", alpha = 0.2) +
  geom_linerange(data = post, aes(x = 2, ymin = quantile(ubed90, 0.025), ymax = quantile(ubed90, 0.975)), 
                 size = 1.1) +
  geom_point(data = post, aes(x = 2, y = quantile(ubed90, 0.5)), 
             size = 5) +
  geom_violin(data = post, aes(x = 3, y = PAR8), fill = "grey", alpha = 0.2) +
  geom_linerange(data = post, aes(x = 3, ymin = quantile(PAR8, 0.025), ymax = quantile(PAR8, 0.975)), 
                 size = 1.1) +
  geom_point(data = post, aes(x = 3, y = quantile(PAR8, 0.5)), 
             size = 5) +
  
  geom_hline(yintercept = 0, linetype ="dashed", size = 1, col = "red") +
  xlab("") +
  theme_bw() +
  scale_x_continuous(breaks = 1:3, labels = c("Temperature", "Benthic stress", "Light"))+
  theme(text = element_text(size = 14 ), axis.text.x = element_text(size = 14)) +
  ggtitle("a") +
  ylab("Coefficient estimates")
Figure4a

Figure4 <- Figure4a /  (Figure4b + Figure4c + ylab("") + Figure4d + ylab(""))

#ggsave("Figure4.tiff", Figure4, dpi = 300)



################################################################################
####                 Realised vs potential upper limits                    #####
################################################################################


dat2 <- dat1[complete.cases(dat1[ , c( "per_suitable", "Zone",  "Latitude", "Longitude") ]), ]
dat2$Reef_spat <- dat2$Reef



# Create covariance matrix from distances between reefs
mat <- matrix(NA, nrow = nrow(dat2), ncol = 2)
mat[ ,1] <- dat2$Longitude
mat[ ,2] <- dat2$Latitude

xy <- latlong2grid(mat)
dat2$x <- xy$x
dat2$y <- xy$y


data2 <- dat2[ , c("Reef", "Zone", "x", "y")]
data2 <- data2[!duplicated(data2), ]


data2 <- data2 %>% group_by(Reef) %>%
  summarise("x" = mean(x, na.rm = TRUE),
            "y" = mean(y, na.rm = TRUE))
data2 <- as.data.frame(data2)
reefs <- unique(data2[order(data2$y), ]$Reef)


spat2 <- matrix(NA, nrow = length(reefs), ncol = length(reefs))
colnames(spat2) <- reefs
row.names(spat2) <- reefs


for ( i in 1: nrow(spat2)) {
  
  d1 <- data2[data2$Reef == rownames(spat2)[i], ]
  
  x1 <- d1$x ; y1 <- d1$y
  
  for (j in 1: ncol(spat2)) {
    
    d2 <- data2[data2$Reef == colnames(spat2)[j], ]
    
    x2 <- d2$x ; y2 <- d2$y
    
    spat2[i,j] <- sqrt((x1 - x2) ^ 2 + (y1 - y2) ^ 2) 
    spat2[i,j] <-  sqrt((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
    
    
    spat2[i,j] <- ifelse(spat2[i,j] > 0,  1/(spat2[i,j])^ 0.5, 1 )
    
    if(spat2[i,j] < 0 ) {print(i); print(j)}
    
  }
}



# Fit models to predict realised upper limit as a function of potential upper limit,
# with and without correcting for spatial autocorrelation
prior <- get_prior(carrying_prop ~  per_suitable  + (1| gr(Reef_spat, cov = spat)),
                   data = dat2,
                   data2 = list(spat = spat2),
                   family = "beta")

m2 <- brm(carrying_prop ~  per_suitable  + (1| gr(Reef_spat, cov = spat)),
          data = dat2,
          data2 = list(spat = spat2),
          family = "beta", chains = 3, cores = 3, 
          prior = prior, iter = 10000, thin = 5,
          control = list(adapt_delta = 0.99, max_treedepth = 17),
          save_pars = save_pars(all = TRUE), file = "m2a")
m2 <- add_criterion(m2, "loo", moment.match = TRUE)


prior <- get_prior(carrying_prop ~  per_suitable  + (1|Reef),
                   data = dat2,
                   family = "beta")
m2a <- brm(carrying_prop ~  per_suitable  + (1| Reef),
           data = dat2,
           family = "beta", chains = 3, cores = 3, 
           prior = prior, iter = 10000, thin = 5,
           control = list(adapt_delta = 0.99, max_treedepth = 17),
           save_pars = save_pars(all = TRUE), file = "m2a_2")

m2a <- add_criterion(m2b, "loo", moment.match = TRUE)


bayes_R2(m2)

loo_compare(m2, m2a)



# Plot predictions of the best-fit model


dat2$rand_m2 <- as.data.frame(fitted(m2, dat2))$Estimate - 
  as.data.frame(fitted(m2, dat2, re_formula = NA))$Estimate  

habitat <- with(dat2, expand.grid( "per_suitable" = seq(min(per_suitable, na.rm = TRUE),
                                                        max(per_suitable, na.rm = TRUE), l = 100)))



habitat <- cbind(habitat, as.data.frame(fitted(m2, habitat, re_formula = NA)))


Figure3 <- ggplot() +
  geom_rect(aes(ymin = 0, ymax = 100, xmin = 80, xmax = 100), col = "red",
            alpha = 0.05, fill = "red", size = 1, linetype = "dotted") +
  geom_point(data = dat2,
             aes(x = per_suitable * 100, y = (carrying_prop - rand_m2) * 100),
             fill = "grey", pch = 21, size = 3) +
  geom_ribbon(data = habitat,
              aes(x = per_suitable * 100, ymin = Q2.5 * 100, ymax = Q97.5 * 100),
              alpha = 0.3) +
  geom_line(data = habitat,
            aes(x = per_suitable * 100, y = Estimate * 100), size = 1.2) +
  geom_abline(intercept = 0, slope = 1 ,linetype ="dashed", size = 1, col = "black") +
  
  xlab("") +
  theme_bw() +
  theme(text = element_text(size = 15 )) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  ylab("Realised upper limit (%)") +
  xlab("Potential upper limit (%)")
Figure3


#ggsave("~/Recovery/Manuscript/Coauthors/Final edits/Figure3.tiff", Figure3, dpi = 300)





################################################################################
####                 Logistic growth example                    #####
################################################################################


# Logistic equation
logistic_eq <- function(K, N, r) {
  
  Nt1 <- N + r * (N * (K - N) / K)
  return(list("N" = Nt1))
  
}


# Create vector with range of realised upper limits, coral covers,
# and intrinsic growth rates to investigate
carrying_prop <- seq(0.01, 1, by = 0.01)
coral_cover <- seq(0.02, 1, by = 0.01)
r <- c(0.2, 1)

# Create empty data frame to save outputs
df3 <- expand.grid("carrying_prop" = numeric(),
                   "coral_cover" = numeric(),
                   "r" = numeric(),
                   "years_to" = integer()) 

# Calculate how many years it takes to go from 1% coral cover to specific 
# coral cover values on the different scenarios
Time <- 1:100

for (cc in carrying_prop) {
  
  for (R in r) {
    
    N <- 0.01
    
    for (t in Time[2:length(Time)]) {
      N[t] <- as.numeric(logistic_eq(K = cc, N = N[t - 1], r = R))
    }
    
    traj <- data.frame( "Time" = Time - 1, "N" = N)
    
    years <- c()
    
    for (y in 1:length(coral_cover)) {
      years_above <- which(traj$N >= coral_cover[y])
      years[y] <- ifelse(length(years_above) > 0 , min(years_above),
                         NA)
    }
    
    df3_sub <- data.frame("carrying_prop" = cc,
                          "coral_cover" = coral_cover,
                          "r" = R,
                          "years_to" = years) 
    
    df3 <- rbind(df3, df3_sub)
    
    
  }
  print(which(carrying_prop == cc))
}



# Plot results

Figure5a <- ggplot()+
  geom_tile(data = df3[df3$r == r[1], ], aes(x = carrying_prop * 100, y = coral_cover * 100, fill = years_to) ) +
  geom_contour(data = df3[df3$r == r[1], ], aes(x = carrying_prop * 100, y = coral_cover * 100, z = years_to),col = "white",
               breaks = seq(0, 80, by = 1), size = 0.5, alpha = 0.7) +
  geom_boxploth(data = dat, aes(x = carrying_prop * 100, y = -3),
                fill = "grey", width  = 5, notch = FALSE, size = 0.7,
                col = "black") +
  theme_classic() +
  scale_fill_viridis(option = "C",
                     name = "Years",  na.value = "black") +
  ylab("Coral cover (%)") +
  xlab("Realised upper limit (%)") +
  theme(text = element_text(size = 13),
        plot.title = element_text(size=14)) +
  ggtitle(expression(paste("a - slow growth (", italic(r), " = 0.2)"))) +
  theme(strip.background =element_rect(fill="white")) +
  scale_x_continuous(limits = c(0,100))



Figure5b <- ggplot()+
  
  geom_tile(data = df3[df3$r == r[2], ], aes(x = carrying_prop * 100, y = coral_cover *100, fill = years_to) ) +
  geom_boxploth(data = dat, aes(x = carrying_prop * 100, y = -3),
                fill = "grey", width  = 5, notch = FALSE, size = 0.7,
                col = "black") +
  geom_contour(data = df3[df3$r == r[2], ], aes(x = carrying_prop * 100, y = coral_cover * 100, z = years_to),col = "white",
               breaks = seq(0, 80, by = 1), size = 0.5, alpha = 0.7) +
  
  theme_classic() +
  scale_fill_viridis(option = "C",
                     name = "Years",  na.value = "black",
                     breaks = seq(3, 14, by = 3)) +
  ylab("Coral cover (%)") +
  xlab("Realised upper limit (%)") +
  theme(text = element_text(size = 13),
        plot.title = element_text(size=14)) +
  ggtitle(expression(paste("b - fast growth (", italic(r), " = 1)"))) +
  theme(strip.background =element_rect(fill="white")) +
  scale_x_continuous(limits = c(0,100))


Figure5 <- (Figure5a + xlab("")) / Figure5b

#ggsave("Figure5.tiff", Figure5, dpi = 300)


# Get examples for discussion
ex1$coral_cover <- as.numeric(ex1$coral_cover)

ex1 <- df3[df3$r == r[1] & df3$carrying_prop > min(dat$carrying_prop) & df3$carrying_prop < max(dat$carrying_prop), ]
min(ex1[ex1$coral_cover > 0.29 & ex1$coral_cover <  0.31 ,]$years_to, na.rm = TRUE)
max(ex1[ex1$coral_cover > 0.29 & ex1$coral_cover <  0.31 ,]$years_to, na.rm = TRUE)

ex2 <- df3[df3$r == r[2] & df3$carrying_prop > min(dat$carrying_prop) & df3$carrying_prop < max(dat$carrying_prop), ]
min(ex2[ex2$coral_cover > 0.29 & ex2$coral_cover <  0.31 ,]$years_to, na.rm = TRUE)
max(ex2[ex2$coral_cover > 0.29 & ex2$coral_cover <  0.31 ,]$years_to, na.rm = TRUE)



#### Check if results are robust to the proxy for wave exposure

prior <- get_prior(carrying_prop ~  st_temp_median + st_PAR8 + 
                     Zone + 
                     (1| gr(Reef_spat, cov = spat)),
                   data = dat,
                   data2 = list(spat = spat),
                   family = "beta")


m1i_check <- brm(carrying_prop ~  st_temp_median + st_PAR8 +
                   Zone  +
                   (1| gr(Reef_spat, cov = spat)),
                 data = dat,
                 data2 = list(spat = spat),
                 family = "beta", chains = 3, cores = 3, 
                 prior = prior, iter = 10000, thin = 5,
                 control = list(adapt_delta = 0.99, max_treedepth = 17),
                 save_pars = save_pars(all = TRUE), file = "m1i_check")



# Plot predictions of the fitted model
Palette <- c("#999999", "#E69F00", "#56B4E9")
names(Palette) <- c("Back", "Flank", "Front")

temp <- expand.grid("st_temp_median" = seq(min(dat$st_temp_median), max(dat$st_temp_median), l = 100),
                    "st_PAR8" = 0,
                    "Zone" = unique(dat$Zone))

temp <- cbind(temp, as.data.frame(fitted(m1i_check, temp, re_formula = NA)))



FigureS2a <- ggplot() +
  geom_point(data = dat, aes( x = st_temp_median, y = (carrying_prop - rand_m1) * 100,
                              fill = Zone, alpha = 0.3), col = "black",
             pch = 21, size = 2) +
  geom_ribbon(data = temp, aes( x = st_temp_median, ymin = Q2.5 * 100, ymax = Q97.5 * 100,
                                fill = Zone, group = Zone), alpha = 0.2) +
  geom_line(data = temp, aes( x = st_temp_median, y = Estimate * 100, col = Zone, group = Zone), size = 1.5) +
  xlab("Temperature (S.D.)") +
  ylab("Realised upper limit (%)") +
  theme_bw() +
  theme(text = element_text(size = 14 )) +
  scale_colour_manual(values = Palette) +
  scale_fill_manual(values = Palette) +
  ggtitle("a")




flow <- expand.grid(
  "st_temp_median" = 0,
  "st_PAR8" = 0,
  "Zone" = unique(dat$Zone))

flow <- cbind(flow, as.data.frame(fitted(m1i_check, flow, re_formula = NA)))



FigureS2b <- ggplot() +
  geom_jitter(data = dat, aes( x = Zone, y = (carrying_prop - rand_m1)* 100, fill = Zone), col = "black",
              pch = 21, size = 2, alpha = 0.3, width = 0.1, height = 0) +
  geom_errorbar(data = flow, aes( x = Zone, ymin = Q2.5 * 100, ymax = Q97.5 * 100,
                                  col = Zone), width = 0.1, size = 1.1) +
  geom_point(data = flow, aes( x = Zone, y = Estimate * 100, fill = Zone),
             size = 5, col = "black", pch = 21) +
  xlab("") +
  scale_x_discrete(breaks = c("Back", "Flank", "Front"),
                   labels = c("Sheltered", "Intermediate", "Exposed")) +
  ylab("Realised upper limit (%)") +
  theme_bw() +
  theme(text = element_text(size = 14 )) +
  scale_colour_manual(values = Palette) +
  scale_fill_manual(values = Palette) +
  xlab("Wave exposure") +
  ggtitle("b")




light <- expand.grid(
  "st_temp_median" = 0,
  "st_PAR8" = seq(min(dat$st_PAR8), max(dat$st_PAR8), l = 100),
  "Zone" = unique(dat$Zone))

light <- cbind(light, as.data.frame(fitted(m1i_check, light, re_formula = NA)))




FigureS2c <- ggplot() +
  geom_point(data = dat, aes( x = st_PAR8, y = (carrying_prop - rand_m1) * 100, fill = Zone), col = "black",
             pch = 21, size = 2, alpha = 0.3) +
  geom_ribbon(data = light, aes( x = st_PAR8, ymin = Q2.5 * 100, ymax = Q97.5 * 100, fill = Zone, group = Zone), alpha = 0.2) +
  geom_line(data = light, aes( x = st_PAR8, y = Estimate * 100, col = Zone, group = Zone), size = 1.1) +
  xlab("Light (S.D)") +
  ylab("Realised upper limit (%)") +
  theme_bw() +
  theme(text = element_text(size = 14 )) +
  scale_colour_manual(values = Palette, breaks = c("Back", "Flank", "Front"),
                      labels = c("Sheltered", "Intermediate", "Exposed")) +
  scale_fill_manual(values = Palette, breaks = c("Back", "Flank", "Front"),
                    labels = c("Sheltered", "Intermediate", "Exposed")) +
  ggtitle("c")




FigureS2 <- (FigureS2a + theme(legend.position = "none"))/ 
  (FigureS2b + theme(legend.position = "none")) /
  FigureS2c + theme(legend.position = "none")

#ggsave("FigureS2.tiff", FigureS2, dpi = 300)

