rm(list = ls())

library(dplyr)
library(ggplot2)
library(patchwork)

dat <- read.csv( "coral_cover.csv", sep = ",")




reefs <- unique(dat$REEF_NAME)

# Get yearly median coral cover for each reef zone 
sdat <- dat %>% group_by(REEF_NAME, YEAR_CODE, ZONE) %>%
  summarise("cover" = median(HARD_COVER, na.rm = TRUE),
            "sd" = sd(HARD_COVER, na.rm = TRUE))
sdat <- as.data.frame(sdat)

# Convert year code to year
sdat$Year <- NA
for (i in 1:nrow(sdat)) {
  
  sdat$Year[i] <- as.numeric(substr(as.character(sdat$YEAR_CODE[i]), start = 1,
                                    stop = 4)) + 0.5
}


sdat$Reef_zone <- paste(sdat$REEF_NAME, sdat$ZONE)

# Get maximum coral cover for each reef zone
max_covers <- sdat %>% 
  dplyr::group_by(REEF_NAME, ZONE) %>%
  dplyr::summarise("max_cover" = max(cover, na.rm = TRUE))
max_covers <- as.data.frame(max_covers)


# Create empty data frame to save yearly coral cover for each reef zone.
# The five-year running mean and whether the estimate for that year was 
# observed or interpolated are included.
running <- data.frame("Reef" = factor(), "Zone" = factor(),
                      "Year" = numeric(), "Cover" = numeric(),
                      "Running" = numeric(), "Estimated" = factor(),
                      "n_years" = numeric(), "range_years" = numeric(),
                      "max_cover" = numeric(), "high_cover" = factor(),
                      "stable" = factor())




for(rz in unique(sdat$Reef_zone)) {
  
  
  sub_dat <- sdat[sdat$Reef_zone == rz, ]
  sub_dat <- sub_dat[order(sub_dat$Year), ]
  
  years <- seq(min(sub_dat$Year), max(sub_dat$Year), by = 1)
  years_m <- sort(sub_dat$Year)
  
  run <- data.frame("Reef" = factor(), "Zone" = factor(),
                    "Year" = numeric(), "Cover" = numeric(),
                    "Running" = numeric(), "Range_running" = numeric(),
                    "Running_covers" = factor(), "Estimated" = factor(),
                    "n_years" = numeric(), "range_years" = numeric())
  
  # Get a data frame that includes all years between surveys 
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
                       "Running" = NA,
                       "Range_running" = NA,
                       "Running_covers" = NA,
                       "Estimated" = ifelse(is.na(cover) == TRUE, "yes", "no"),
                       "n_years" = NA,
                       "range_years" = NA)
    run <- rbind(run, run1)
    
  }
  
  
  # Interpolate missing coral covers
  covers <- is.na(run$Cover)
  true_covers <- which(covers == FALSE)
  
  run$stable <- NA
  
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
  
  # Get the 5yr-running mean and record how many years out of the five had 
  # observed data (vs. interpolated) (n_years)
  if (nrow(run) > 4) {
    
    for (i in 3:(nrow(run)-2)) {
      run$Running[i] <- mean(run$Cover[(i-2) : (i+2)])
      run$Range_running[i] <- max(run$Cover[(i-2) : (i+2)]) - 
        min(run$Cover[(i-2) : (i+2)])
      run$Running_covers[i] <- paste(round(run$Cover[(i-2) : (i+2)], 0), collapse = ",")
      run$n_years[i] <- length(true_covers[true_covers >= (i-2) &
                                             true_covers <= (i+2)])
      min_year <- run$Year[min(true_covers[true_covers >= (i-2)])]
      max_year <- run$Year[max(true_covers[true_covers <= (i+2)])]
      run$range_years[i] <- max_year - min_year
      
    }
  }
  
  # Get the maximum of the 5-year running mean
  run$max_cover <- max(run$Running, na.rm = TRUE)
  run$max_cover <- ifelse(run$max_cover == -Inf, NA, run$max_cover)
  
  # Mark in which years the running mean is at least at 80% of its maximum
  # value.
  run$high_cover <- ifelse(run$Running / run$max_cover > 0.8, "yes", "no")
  
  # Mark whether the coral cover is stable (at 80% or more of the maximum coral
  # cover for at least five years)
  run$stable <- NA
  
  if (nrow(run) > 4) {
    
    for (i in 3:(nrow(run) - 2)) {
      run$stable[i] <- ifelse(run$high_cover[i - 2] == "yes" &
                                run$high_cover[i - 1] == "yes" &
                                run$high_cover[i] == "yes" &
                                run$high_cover[i + 1] == "yes" &
                                run$high_cover[i + 2] == "yes", "yes", "no")
    }
  }
  
  running <- rbind(running, run)
  
  
}



running <- running[complete.cases(running$stable), ]

# Get the realised upper limit. It needs to have five continuous years 
# where the running mean is at least 80% than the maximum running mean and
# the five-year running mean needs to have been estimated with at least three
# years of observed data (and a maximum of two years interpolated)
run1 <- running[running$stable == "yes" & running$n_years > 2, ] %>%
  group_by(Reef,  Zone) %>%
  summarise("carrying" = max(Running))


run1 <- as.data.frame(run1)



run1 <- merge(run1[ , c("Reef", "Zone", "carrying")], running[ , c("Reef", "Zone", "Year","Running",
                                                                   "Range_running", "Running_covers",
                                                                   "n_years", "range_years",
                                                                   "stable")],
              by.x = c("Reef", "Zone"), 
              by.y = c("Reef",  "Zone"), all.y = TRUE)







# Exclude reefs or zones within reefs where the reef zones where not consistently
# defined.

run1 <- run1[!(run1$Reef == "15077S" & run1$Zone == "Back"),  ]
run1 <- run1[!run1$Reef == "HELSDON REEF", ]
run1 <- run1[!(run1$Reef == "ST CRISPIN REEF" & run1$Zone %in% c("Flank1", "Front")), ]
run1 <- run1[!(run1$Reef == "THETFORD REEF" & run1$Zone == "Flank1"),  ]
run1 <- run1[!(run1$Reef == "CHICKEN REEF" & run1$Zone == "Flank1"),  ]
run1 <- run1[!(run1$Reef == "DIP REEF" & run1$Zone %in% c("Flank1", "Back", "Front")), ]
run1 <- run1[!(run1$Reef == "RIB REEF" & run1$Zone %in% c("Flank1", "Back", "Flank1")), ]
run1 <- run1[!(run1$Reef == "ERSKINE REEF" & run1$Zone %in% c("Front", "Flank1")), ]
run1 <- run1[!(run1$Reef == "FARQUHARSON REEF (NO 1)" & run1$Zone == "Flank1"),  ]
run1 <- run1[!run1$Reef == "21302S", ]
run1 <- run1[!(run1$Reef == "22084S" & run1$Zone == "Front"),  ]
run1 <- run1[!(run1$Reef == "ONE TREE REEF" & run1$Zone == "Flank1"),  ]
run1 <- run1[!(run1$Reef == "20354S" & run1$Zone == "Flank1"),  ]
run1 <- run1[!(run1$Reef == "MOORE REEF" & run1$Zone == "Flank1"),  ]
run1 <- run1[!(run1$Reef == "OSBORNE REEF" & run1$Zone %in% c("Back", "Flank1")), ]
run1 <- run1[!run1$Reef == "ASHMORE BANKS (2)", ]



sdat2 <- run1 %>% group_by(Reef, Zone) %>%
  summarise("carrying" = max(carrying, na.rm = TRUE))

sdat2 <- sdat2[!sdat2$carrying == -Inf, ]


#write.table(sdat, "realised_upper_limits.csv", sep = ",", row.names = FALSE) 




# Plot example for Figure 1

reef = "21529S"


Sheltered <- ggplot() +
  geom_hline(data = run1[run1$Reef == reef & run1$Zone == "Back", ],
             aes(yintercept = carrying), col = "black", alpha = 1, linetype = "dashed", size = 1.2) +
  geom_line(data = running[running$Reef == reef & running$Zone == "Back", ], 
            aes(x = Year, y = Running), col =  "#999999", size = 1, alpha = 0.8) +
  geom_point(data = running[running$Reef == reef & running$Estimated == "no" &
                              running$Zone == "Back", ], 
             aes(x = Year, y= Cover), fill =  "#999999", size = 3 ,pch = 21, col = "black", alpha = 0.5) +
  ylab("") +
  xlab("") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25),
                     labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme(text = element_text(size = 13))
Sheltered

Intermediate1 <- ggplot() +
  geom_hline(data = run1[run1$Reef == reef & run1$Zone == "Flank1", ],
             aes(yintercept = carrying), col = "black", alpha = 1, linetype = "dashed", size = 1.2) +
  geom_line(data = running[running$Reef == reef & running$Zone == "Flank1", ], 
            aes(x = Year, y = Running), col =  "#E69F00", size = 1, alpha = 0.8) +
  geom_point(data = running[running$Reef == reef & running$Estimated == "no" &
                              running$Zone == "Flank1", ], 
             aes(x = Year, y= Cover), fill =  "#E69F00", size = 3 ,pch = 21, col = "black", alpha = 0.5) +
  ylab("") +
  xlab("") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25),
                     labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme(text = element_text(size = 13))
Intermediate1


Intermediate2 <- ggplot() +
  geom_hline(data = run1[run1$Reef == reef & run1$Zone == "Flank2", ],
             aes(yintercept = carrying), col = "black", alpha = 1, linetype = "dashed", size = 1.2) +
  geom_line(data = running[running$Reef == reef & running$Zone == "Flank2", ], 
            aes(x = Year, y = Running), col =  "#E69F00", size = 1, alpha = 0.8) +
  geom_point(data = running[running$Reef == reef & running$Estimated == "no" &
                              running$Zone == "Flank2", ], 
             aes(x = Year, y= Cover), fill =  "#E69F00", size = 3 ,pch = 21, col = "black", alpha = 0.5) +
  ylab("") +
  xlab("") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25),
                     labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme(text = element_text(size = 13))
Intermediate2


Exposed <- ggplot() +
  geom_hline(data = run1[run1$Reef == reef & run1$Zone == "Front", ],
             aes(yintercept = carrying), col = "black", alpha = 1, linetype = "dashed", size = 1.2) +
  geom_line(data = running[running$Reef == reef & running$Zone == "Front", ], 
            aes(x = Year, y = Running), col = "#56B4E9", size = 1, alpha = 0.8) +
  geom_point(data = running[running$Reef == reef & running$Estimated == "no" &
                              running$Zone == "Front", ], 
             aes(x = Year, y= Cover), fill = "#56B4E9", size = 3 ,pch = 21, col = "black", alpha = 0.5) +
  ylab("") +
  xlab("") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25),
                     labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme(text = element_text(size = 13))
Exposed


Sheltered / Intermediate1 / Exposed / Intermediate2


