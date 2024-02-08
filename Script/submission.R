## Submission

###########################################################
#################### BUOY INFO ############################
###########################################################

setwd("../Data")
#data from: https://www.aoml.noaa.gov/phod/gdp/interpolated/data/subset.php
#velocity units are cm/s https://www.aoml.noaa.gov/phod/gdp/buoydata_header.php

library(ggplot2)

dat <- read.csv("carib_buoy_dat.csv")    ### Buoy information
dat <- dat[which(dat$id < 100000000000),] #buoys above this number is faulty.

start_time <- as.POSIXct(strptime("1/1/1979 1:00:00", "%m/%d/%Y %H:%M:%S"))

dat$timeB <- paste(dat$date, dat$time)     

dat$time_new <- as.POSIXct(strptime(dat$timeB, "%m/%d/%Y %H:%M:%S"))
dat$year <- format(strptime(dat$timeB, "%m/%d/%Y %H:%M:%S"), '%Y') ## Separating Years
dat$year <- as.numeric(dat$year)
dat$time_standard <- as.numeric(dat$time_new - start_time)

dat <- dat[which(dat$speed < 999),]

## Determining areas where the buoy enters the region   ("Figure "Buoy_Start.jpg")
length(unique(dat$id))

df_temp <- data.frame(matrix(ncol=4, nrow= length(unique(dat$id))))
colnames(df_temp) <- c('id', 'lat', 'lon', 'time')

buoy_id <- unique(dat$id)

for(i in c(1:length(unique(dat$id)))){
  dat_temp <- dat[which(dat$id == buoy_id[i]),]
  time_temp <- min(dat_temp$time_standard)
  lat_temp <- dat_temp$lat[which(dat_temp$time_standard == time_temp)]
  lon_temp <- dat_temp$lon[which(dat_temp$time_standard == time_temp)]
  time_temp2 <- dat_temp$time_new[which(dat_temp$time_standard == time_temp)]
  
  df_temp$lat[i] <- lat_temp
  df_temp$lon[i] <- lon_temp
  df_temp$time[i] <- time_temp2
}

length(df_temp$id[which(df_temp$lat > 8)])

#write.csv(df_temp, "buoy_start.csv")

############################################################
## Selecting points of interest and defining search areas ##
############################################################

# Points of interest
points_dat <- read.csv("point_locations_island.csv")
trial_point_dat_temp <- points_dat[which(points_dat$lat1 > 0),]

# Examining for box size
# need to get ensemble mean for n and e, variance, and retention time

trial_point_dat_temp$lat1a <- NA
trial_point_dat_temp$lat2a <- NA
trial_point_dat_temp$long1a <- NA
trial_point_dat_temp$long2a <- NA

for(i in c(1:length(trial_point_dat_temp$ID))){
  fit_temp <- subset(dat, lat <= trial_point_dat_temp$lat2[i] & lat >= trial_point_dat_temp$lat1[i] & lon <= trial_point_dat_temp$long1[i] & lon >= trial_point_dat_temp$long2[i])
  list <- unique(fit_temp$id)
  
  if(length(fit_temp$id) > 0){
    # Creating a new dataset for just the averages
    data_average <- data.frame(matrix(ncol=6, nrow=length(list)))
    colnames(data_average) <- c('id', 'mean_ve', 'mean_vn', 'var_ve', 'var_vn', 'retention')
    
    for(j in c(1:length(list))){
      data_average$id[j] <- list[j]
      data_average$mean_ve[j] <- mean(fit_temp$ve[which(fit_temp$id == list[j])])*864 #864 to convert from cm/s to m/day
      data_average$var_ve[j] <- var(fit_temp$ve[which(fit_temp$id == list[j])])*864
      data_average$mean_vn[j] <- mean(fit_temp$vn[which(fit_temp$id == list[j])])*864
      data_average$var_vn[j] <- var(fit_temp$vn[which(fit_temp$id == list[j])])*864
      data_average$retention[j] <- as.numeric(max(fit_temp$time_standard[which(fit_temp$id == list[j])]) - min(fit_temp$time_standard[which(fit_temp$id == list[j])]))
    }
  
    mean_e <- mean(data_average$mean_ve)
    mean_n <- mean(data_average$mean_vn)
    mean_t <- mean(data_average$retention) # days
    
    x1 <- abs((mean_e * mean_t) / 111000)     #Converting from degrees to meters
    y1 <- abs((mean_n * mean_t) / 111000)
    
    trial_point_dat_temp$long1a[i] <- trial_point_dat_temp$long[i] - x1
    trial_point_dat_temp$long2a[i] <- trial_point_dat_temp$long[i] + x1
    trial_point_dat_temp$lat1a[i] <- trial_point_dat_temp$lat[i] - y1
    trial_point_dat_temp$lat2a[i] <- trial_point_dat_temp$lat[i] + y1
  } else {
    print("not enough data")
  }
}   

#write.csv(trial_point_dat_temp, "extended_point_dat.csv")


### Getting number of unique buoys that go by each point
points_dat$buoy_count <- -1

library(dplyr)

for(i in 1:length(points_dat$buoy_count)){
  
  point1 <- subset(dat, lat >= trial_point_dat_temp$lat1a[i] & lat <= trial_point_dat_temp$lat2a[i] & lon <= trial_point_dat_temp$long2a[i] & lon >= trial_point_dat_temp$long1a[i])    #starting point
  points_dat$buoy_count[i] <- length(unique(point1$id))
}

#write.csv(points_dat, "points_summary.csv")

############################################################
### Getting island-specific buoy information ##############
##########################################################

# Need to search all boxes from each island at the same time to get how many unique buoys start from there.

trial_point_dat <- read.csv("extended_point_dat.csv")
points_dat <- read.csv("points_summary.csv")

trial_point_dat <- trial_point_dat[which(trial_point_dat$lat1a > 0),]
trial_point_dat$temp <- trial_point_dat$lat-trial_point_dat$lat1a
trial_point_dat <- trial_point_dat[which(trial_point_dat$temp < 1.5),]

islands <- unique(trial_point_dat$Island)

island_dat <- data.frame(matrix(ncol=2, nrow=length(islands)))
colnames(island_dat) <- c('island', 'count')

island_dat$island <- islands

for ( i in 1:length(islands)){
  temp_id <- points_dat$ID[which(points_dat$Island == islands[i])]     ## Getting specific id codes for each island
  
  k <- temp_id[1]
  
  point1a <- subset(dat, lat >= trial_point_dat$lat1a[which(trial_point_dat$ID == k)] & lat <= trial_point_dat$lat2a[which(trial_point_dat$ID == k)] & lon <= trial_point_dat$long2a[which(trial_point_dat$ID == k)] & lon >= trial_point_dat$long1a[which(trial_point_dat$ID == k)])
  
  for(j in c(temp_id)){
    
    point1b <- subset(dat, lat >= trial_point_dat$lat1a[which(trial_point_dat$ID == j)] & lat <= trial_point_dat$lat2a[which(trial_point_dat$ID == j)] & lon <= trial_point_dat$long2a[which(trial_point_dat$ID == j)] & lon >= trial_point_dat$long1a[which(trial_point_dat$ID == j)])
    point1a <- rbind(point1a, point1b)
  }
  
  point1c <- point1a[!duplicated(point1a),]
  island_dat$count[i] <- length(unique(point1c$id))
}

#write.csv(island_dat, "island_source_count.csv")

# Then, I need to determine how many of those unique buoys make it to any of the areas in another island.

dat.perc_temp <- read.csv("perc_floating_dat.csv")
islands <- island_dat$island

#c barb, c. spissa, g. attenuata, and s. amara did not float more than one day, so those are removed.
dat.perc <- dat.perc_temp[which(dat.perc_temp$Species != "Coccothrinax barbadensis" & dat.perc_temp$Species != "Coccothrinax spissa" & dat.perc_temp$Species != "Gaussia attenuata" & dat.perc_temp$Species != "Syagrus amara"),]
dat.perc$prob <- dat.perc$perc/100

plants <- unique(dat.perc$Species)

island_dat <- read.csv("island_source_count.csv")

# Then, I need to determine how many of those unique buoys make it to any of the areas in another island.
island_matrix_count <- matrix(-1, nrow= length(islands), ncol= length(islands), dimnames = list(c(islands), c(islands)))
island_matrix_min <- matrix(-1, nrow= length(islands), ncol= length(islands), dimnames = list(c(islands), c(islands)))
island_matrix_mean <- matrix(-1, nrow= length(islands), ncol= length(islands), dimnames = list(c(islands), c(islands)))
island_matrix_median <- matrix(-1, nrow= length(islands), ncol= length(islands), dimnames = list(c(islands), c(islands)))
island_matrix_max <- matrix(-1, nrow= length(islands), ncol= length(islands), dimnames = list(c(islands), c(islands)))
island_matrix_prob <- matrix(-1, nrow= length(islands), ncol= length(islands), dimnames = list(c(islands), c(islands)))

for(c in 1:length(plants)){
  
  assign(paste0(plants[c],"_connectivity_matrix"),matrix(-1, nrow= length(islands), ncol= length(islands), dimnames = list(c(islands), c(islands))))
  
}

list_matrix <- mget(grep("_connectivity_matrix", 
                         names(which(unlist(eapply(.GlobalEnv,is.matrix)))), 
                         value = TRUE))

plants_list <- c("Pseudophoenix vinifera", "Thrinax radiata", "Aiphanes minima", "Coccothrinax borhidiana",  "Zombia antillarum", "Acrocomia crispa", "Chrysobalanus icaco", "Goetzae elegans", "Theophrasta jussieui", "Catesbaea spinosa")

for (i in 1:length(islands)){
  temp_id_a <- points_dat$ID[which(points_dat$Island == islands[i])]     ## Getting specific id codes for each island
  
  z <- temp_id_a[1]
  point1a <- subset(dat, lat >= trial_point_dat$lat1a[which(trial_point_dat$ID == z)] & lat <= trial_point_dat$lat2a[which(trial_point_dat$ID == z)] & lon <= trial_point_dat$long2a[which(trial_point_dat$ID == z)] & lon >= trial_point_dat$long1a[which(trial_point_dat$ID == z)])
  
  for(j in c(temp_id_a)){
    
    point1b <- subset(dat, lat >= trial_point_dat$lat1a[which(trial_point_dat$ID == j)] & lat <= trial_point_dat$lat2a[which(trial_point_dat$ID == j)] & lon <= trial_point_dat$long2a[which(trial_point_dat$ID == j)] & lon >= trial_point_dat$long1a[which(trial_point_dat$ID == j)])
    point1a <- rbind(point1a, point1b)
    
    islandA_point <- point1a    ### ISLAND 1
  }
  
  for(k in 1:length(islands)){
    
    temp_id_b <- points_dat$ID[which(points_dat$Island == islands[k])]
    
    y <- temp_id_b[1]
    point2a <- subset(dat, lat >= trial_point_dat$lat1a[which(trial_point_dat$ID == y)] & lat <= trial_point_dat$lat2a[which(trial_point_dat$ID == y)] & lon <= trial_point_dat$long2a[which(trial_point_dat$ID == y)] & lon >= trial_point_dat$long1a[which(trial_point_dat$ID == y)])
    
    for(l in c(temp_id_b)){
      
      point2b <- subset(dat, lat >= trial_point_dat$lat1a[which(trial_point_dat$ID == l)] & lat <= trial_point_dat$lat2a[which(trial_point_dat$ID == l)] & lon <= trial_point_dat$long2a[which(trial_point_dat$ID == l)] & lon >= trial_point_dat$long1a[which(trial_point_dat$ID == l)])
      point2a <- rbind(point2a, point2b)
      
      islandB_point <- point2a    ### ISLAND 2
    }
    
    # Looking for common buoys 
    AB <- intersect(unique(islandA_point$id), unique(islandB_point$id))           #Look for common buoys, if it isn't identical
    
    if(length(AB) == 0){
      
      island_matrix_count[k, i] <- NA
      island_matrix_min[k, i] <- NA
      island_matrix_mean[k, i] <- NA
      island_matrix_median[k, i] <- NA
      island_matrix_max[k, i] <- NA
      island_matrix_prob[k, i] <- NA
      
    } else {
      
      if(identical(islandA_point, islandB_point) == TRUE){        #If everything is identical, count as 0 for everything except count
        
        island_matrix_count[k, i] <- NA
        island_matrix_min[k, i] <- NA
        island_matrix_mean[k, i] <- NA
        island_matrix_median[k, i] <- NA
        island_matrix_max[k, i] <- NA
        island_matrix_prob[k, i] <- NA
        
      } else {
        
        if(length(AB) > 0){            #If there are common buoys, determine transit time between the two points
          
          ######## STARTING POINT (POINT 1)
          point1A <- islandA_point[islandA_point$id %in% AB,]           #Only getting points within intersecting buoys
          point1A <- point1A %>% arrange(time_standard)           #Sorting by time
          point1A <- point1A %>% mutate(time_diff = time_standard - lag(time_standard))         #Calculating differences in consecutive times
          point1A <- point1A %>% replace(is.na(.), 0)               #The first value will be NA for time_diff. Making that 0.
          
          ##### Now need to group them.
          
          point1A$ID <- c(1:length(point1A$id))      #Creating a list of IDs
          point1A$code <- ifelse(abs(point1A$time_diff) < 7, "A", "B")           #If consecutive times are more than 7 days apart, considering that as a separate event. B marks breaking points. 
          point1A$group <- NA              # Creating an empty column
          
          temp1A <- point1A$ID[which(point1A$code == "B")]                   #Finding which IDs were the breaking points
          temp2A <- sort(append(temp1A, c(1, length(point1A$id))))           #Need to add in 1 and the final value 
          
          group_length <- length(temp1A)                                     #determining group ids
          group_number <- c(1:group_length, group_length+1)
          
          for(p in 1:c(length(temp2A)-1)){
            point1A$group[temp2A[p]:c(temp2A[p+1]-1)] <- group_number[p]
            point1A$group[max(temp2A)] <- group_number[p]
          }
          
          ######## FINISHING POINT (POINT 2)
          point2A <- islandB_point[islandB_point$id %in% AB,]           #Only getting points within intersecting buoys
          point2A <- point2A %>% arrange(time_standard)           #Sorting by time
          point2A <- point2A %>% mutate(time_diff = time_standard - lag(time_standard))         #Calculating differences in consecutive times
          point2A <- point2A %>% replace(is.na(.), 0)               #The first value will be NA for time_diff. Making that 0.
          
          ##### Now need to group them.
          
          point2A$ID <- c(1:length(point2A$id))      #Creating a list of IDs
          point2A$code <- ifelse(abs(point2A$time_diff) < 7, "A", "B")           #If consecutive times are more than 7 days apart, considering that as a separate event. B marks breaking points. 
          point2A$group <- NA              # Creating an empty column
          
          temp1A <- point2A$ID[which(point2A$code == "B")]                   #Finding which IDs were the breaking points
          temp2A <- sort(append(temp1A, c(1, length(point2A$id))))           #Need to add in 1 and the final value 
          
          group_length <- length(temp1A)                                     #determining group ids
          group_number <- c(1:group_length, group_length+1)
          
          for(p in 1:c(length(temp2A)-1)){
            point2A$group[temp2A[p]:c(temp2A[p+1]-1)] <- group_number[p]
            point2A$group[max(temp2A)] <- group_number[p]
          }
          
          #########
          
          temp_min <- c(rep(-1, length(AB)))   #how many of the buoys from island 1 make it to island 2?
          
          for(q in 1:length(AB)){
            point1B <- point1A[which(point1A$id == AB[q]),]   # Sorting by matching buoys
            point2B <- point2A[which(point2A$id == AB[q]),]
            
            time_diff_matrix <- matrix(-1, nrow=length(unique(point2B$group)) , ncol= length(unique(point1B$group)), dimnames = list(c(unique(point2B$group)), c(unique(point1B$group))))   #empty matrix of 
            
            for(m in 1: length(unique(point1B$group))){
              
              for(n in 1: length(unique(point2B$group))){
                temp1B <- point1B[which(point1B$group == unique(point1B$group)[m]),]
                max_time1 <- max(temp1B$time_standard)
                
                temp2B <- point2B[which(point2B$group == unique(point2B$group)[n]),]
                min_time2 <- min(temp2B$time_standard)
                
                time_diff_matrix[n,m] <- min_time2 - max_time1
              }
            }
            
            time_diff_matrix <- replace(time_diff_matrix, which(time_diff_matrix < 0), NA)
            time_diff_df <- data.frame(time_diff_matrix)
            time_diff_df$temp <- NA
            
            list_temp <- c(rep(-1,ncol(time_diff_matrix)))
            
            for(o in 1:c(ncol(time_diff_df)-1)){
              temp_value1 <- min(time_diff_df[,o], na.rm=TRUE)       #can introduce Inf values. Will turn that into NAs later.
              temp_value2 <- min(time_diff_df[,o+1], na.rm=TRUE)
              
              if(temp_value1 > temp_value2){
                list_temp[o] <- NA
              } else {
                list_temp[o] <- temp_value1
              }
            } 
            
            temp_min[q] <- min(list_temp, na.rm=TRUE)  #The minimum time it takes to travel from one point to another for all common buoys detected between the two areas
          }
          
          temp_min[is.infinite(temp_min)] <- NA   #Fixing problem with Inf from above
          
          if(sum(is.na(temp_min)) == length(temp_min)){    #Fixing an error when the list is nothing but NAs
            
            island_matrix_count[k, i] <- NA
            island_matrix_min[k, i] <- NA
            island_matrix_mean[k, i] <- NA
            island_matrix_median[k, i] <- NA
            island_matrix_max[k, i] <- NA
            island_matrix_prob[k, i] <- NA
            
          } else {
            
            island_matrix_count[k, i] <- length(temp_min[!is.na(temp_min)])
            island_matrix_min[k, i] <- min(temp_min, na.rm=TRUE)
            island_matrix_mean[k, i] <- mean(temp_min, na.rm=TRUE)
            island_matrix_median[k, i] <- median(temp_min, na.rm=TRUE)
            island_matrix_max[k, i] <- max(temp_min, na.rm=TRUE)
            
            #ggplot() + aes(temp_min) + geom_histogram(binwidth=1, color = "black", fill = "white") + xlim(c(-1, max(temp_min)*2))+ labs(x="transit time (days)")
            #ggsave(paste0("dispersal_",islands[i],"_to_", islands[k], ".png"))
            
            temp_min2 <- ceiling(temp_min[!is.na(temp_min)])
            
            source_count_temp <- island_dat$count[which(island_dat$island == islands[i])]
            temp_min3 <- unique(sort(temp_min2, decreasing = FALSE))
            temp_table <- as.data.frame(table(temp_min2))
            
            temp_table$disp_prob <- temp_table$Freq/source_count_temp
            temp_table$time <- temp_min3    #Not sure why I have to do this. R isn't working properly with getting temp_min2 to show up without this.
            
            island_matrix_prob[k, i] <- sum(temp_table$disp_prob)
            
            for(a in 1:length(plants_list)){
              
              temp_plant <- dat.perc[which(dat.perc$Species == plants_list[a]),]
              
              for(b in 1:length(temp_table$temp_min2)){
                
                if(temp_table$time[b] <= max(temp_plant$Days)){
                  temp_table$survival_prob[b] <- temp_plant$prob[which(temp_plant$Days == temp_table$time[b])]
                } else {
                  
                  temp_table$survival_prob[b] <- NA
                }
                
              }
              
              temp_table$connectivity <- temp_table$disp_prob*temp_table$survival_prob
              
              list_matrix[[a]][k,i] <- sum(temp_table$connectivity, na.rm=TRUE)  
              
            }
            
          }
          
        } else {
          
          island_matrix_count[k, i] <- NA
          island_matrix_min[k, i] <- NA
          island_matrix_mean[k, i] <- NA
          island_matrix_median[k, i] <- NA
          island_matrix_max[k, i] <- NA
          island_matrix_prob[k, i] <- NA
        } 
        
      }
      
    } 
    
  }  
} 


write.csv(island_matrix_count, "matrix_count_island.csv")
write.csv(island_matrix_min, "matrix_min_island.csv")
write.csv(island_matrix_mean, "matrix_mean_island.csv")
write.csv(island_matrix_median, "matrix_median_island.csv")
write.csv(island_matrix_max, "matrix_max_island.csv")
write.csv(island_matrix_prob, "matrix_prob_island.csv")

temp_matrix1 <- as.matrix(list_matrix[[1]])
temp_matrix2 <- as.matrix(list_matrix[[2]])
temp_matrix3 <- as.matrix(list_matrix[[3]])
temp_matrix4 <- as.matrix(list_matrix[[4]])
temp_matrix5 <- as.matrix(list_matrix[[5]])
temp_matrix6 <- as.matrix(list_matrix[[6]])
temp_matrix7 <- as.matrix(list_matrix[[7]])
temp_matrix8 <- as.matrix(list_matrix[[8]])
temp_matrix9 <- as.matrix(list_matrix[[9]])
temp_matrix10 <- as.matrix(list_matrix[[10]])

write.csv(temp_matrix1, paste0(plants_list[1],"_connectivity_matrix.csv", ".csv"))
write.csv(temp_matrix2, paste0(plants_list[2],"_connectivity_matrix.csv", ".csv"))
write.csv(temp_matrix3, paste0(plants_list[3],"_connectivity_matrix.csv", ".csv"))
write.csv(temp_matrix4, paste0(plants_list[4],"_connectivity_matrix.csv", ".csv"))
write.csv(temp_matrix5, paste0(plants_list[5],"_connectivity_matrix.csv", ".csv"))
write.csv(temp_matrix6, paste0(plants_list[6],"_connectivity_matrix.csv", ".csv"))
write.csv(temp_matrix7, paste0(plants_list[7],"_connectivity_matrix.csv", ".csv"))
write.csv(temp_matrix8, paste0(plants_list[8],"_connectivity_matrix.csv", ".csv"))
write.csv(temp_matrix9, paste0(plants_list[9],"_connectivity_matrix.csv", ".csv"))
write.csv(temp_matrix10, paste0(plants_list[10],"_connectivity_matrix.csv", ".csv"))

########## CALCULATING PROBABILITY OF BUOYS MAKING IT FROM ONE ISLAND TO ANOTHER 

matrix_temp <- read.csv("matrix_count_island.csv")
matrix_temp$X <- NULL

source_buoys_temp <- read.csv("island_source_count.csv")
source_buoys <- source_buoys_temp$count

prob_matrix <- t(t(matrix_temp)/ source_buoys)

write.csv(prob_matrix, "prob_matrix_island.csv")


########### CALCULATING DISTANCE BETWEEN ISLANDS

# Points of interest
points_dat <- read.csv("point_locations_island.csv")
trial_point_dat <- points_dat[which(points_dat$lat1 > 0),]

#install.packages("geosphere")
library(geosphere)

distance_matrix <- matrix(-1, nrow= length(trial_point_dat$ID), ncol= length(trial_point_dat$ID), dimnames = list(c(trial_point_dat$ID), c(trial_point_dat$ID)))

for(k in 1:length(trial_point_dat$ID)){
  
  temp_point_1 <- c(trial_point_dat$long[k], trial_point_dat$lat[k])
  
  for(l in 1:length(trial_point_dat$ID)){
    temp_point_2 <- c(trial_point_dat$long[l], trial_point_dat$lat[l])
    
    point1 <- temp_point_1
    point2 <- temp_point_2
    
    if(identical(point1, point2) == TRUE){
      
      distance_matrix[l,k] <- 0
      
    } else {
      
      # Calculating distance between two points
      
      distance_matrix[l,k] <- distm(point1, point2, fun = distHaversine)
    }
  }
}

#write.csv(distance_matrix, "matrix_dist.csv")

matrix_dat_temp <- read.csv("matrix_dist.csv")
matrix_dat <- matrix_dat_temp[, -c(1)]

islands <- unique(points_dat$Island)     ## Getting unique island names

dist_matrix <- matrix(-1, nrow= length(islands), ncol= length(islands), dimnames = list(c(islands), c(islands)))  ## Creating empty matrix

for ( i in 1:length(islands)){
  temp_id <- points_dat$ID[which(points_dat$Island == islands[i])]     ## Getting specific id codes for each island
  
  temp_column1 <- matrix_dat[,c(temp_id)]                           ## Isolating travel data for each island              
  temp_column1$min <- apply(temp_column1, 1, FUN = min, na.rm=TRUE)   ## Obtaining minimum travel distance for each row
  
  for(j in 1:length(islands)){
    temp_id2 <- points_dat$ID[which(points_dat$Island == islands[j])]
    temp_row1 <- temp_column1[c(temp_id2),]
    
    dist_matrix[i,j] <- min(temp_row1$min)
  }
}     

#write.csv(dist_matrix, "matrix_dist_summary.csv")

#####################################################################################
## Determining how many buoys go by all of the islands and for how long

trial_point_dat <- read.csv("extended_point_dat.csv")
points_dat <- read.csv("points_summary.csv")

trial_point_dat <- trial_point_dat[which(trial_point_dat$lat1a > 0),]
trial_point_dat$temp <- trial_point_dat$lat-trial_point_dat$lat1a
trial_point_dat <- trial_point_dat[which(trial_point_dat$temp < 1.5),]

islands <- unique(trial_point_dat$Island)

island_dat <- data.frame(matrix(ncol=2, nrow=length(islands)))
colnames(island_dat) <- c('island', 'count')

island_dat$island <- islands

dat$temp_id <- c(1:length(dat$id))

#creating empty vector
list_buoy <- c()
list_time <- c()
list_id <- c()   #Listing min and max temp_ids for each buoy that goes through the study site could give us a way to pinpoint the times at which they are in our study region.

islands <- unique(points_dat$Island)  

for (i in 1:length(islands)){
  temp_id <- points_dat$ID[which(points_dat$Island == islands[i])]     ## Getting specific id codes for each island
  
  p <- temp_id[1]
  z <- which(trial_point_dat$ID == p)
  point1a <- subset(dat, lat >= trial_point_dat$lat1a[z] & lat <= trial_point_dat$lat2a[z] & lon <= trial_point_dat$long2a[z] & lon >= trial_point_dat$long1a[z])
  
  for(k in temp_id){
    
    j <- which(trial_point_dat$ID == k)
    
    point1b <- subset(dat, lat >= trial_point_dat$lat1a[j] & lat <= trial_point_dat$lat2a[j] & lon <= trial_point_dat$long2a[j] & lon >= trial_point_dat$long1a[j])
    point1a <- rbind(point1a, point1b)
    
    islandA_point <- point1a
  }
  
  buoys <- unique(point1a$id)
  times <- c(min(point1a$time_new), max(point1a$time_new))
  list_buoy <- append(list_buoy, buoys)
  list_time <- append(list_time, times)
}

length(unique(list_buoy))   #832 buoys going through our islands
time_temp <- list_time[!is.infinite(list_time)] 
max(time_temp) - min(time_temp)    #10505.75 days in the study site

study_buoys <- unique(list_buoy)


study_dat <- subset(dat, lat >= trial_point_dat$lat1a[1] & lat <= trial_point_dat$lat2a[1] & lon <= trial_point_dat$long2a[1] & lon >= trial_point_dat$long1a[1])

for(i in 2:length(trial_point_dat$ID)){
  point1a <- subset(dat, lat >= trial_point_dat$lat1a[i] & lat <= trial_point_dat$lat2a[i] & lon <= trial_point_dat$long2a[i] & lon >= trial_point_dat$long1a[i])
  study_dat <- rbind(study_dat, point1a)
}

study_buoys <- unique(study_dat$id)

## Determining the average time that each buoy took to go through our region

study_dat <- subset(dat, lat >= trial_point_dat$lat1a[1] & lat <= trial_point_dat$lat2a[1] & lon <= trial_point_dat$long2a[1] & lon >= trial_point_dat$long1a[1])

for(i in 2:length(trial_point_dat$ID)){
  point1a <- subset(dat, lat >= trial_point_dat$lat1a[i] & lat <= trial_point_dat$lat2a[i] & lon <= trial_point_dat$long2a[i] & lon >= trial_point_dat$long1a[i])
  study_dat <- rbind(study_dat, point1a)
}

study_buoys <- unique(study_dat$id)

buoy_df <- data.frame(matrix(ncol=3, nrow=length(study_buoys)))
colnames(buoy_df) <- c('id', 'min','max')

for(i in 1:length(study_buoys)){
  buoy_df$id[i] <- study_buoys[i]
  
  temp_dat <- study_dat[which(study_dat$id == study_buoys[i]),]
  
  buoy_df$min[i] <- min(temp_dat$time_standard)
  buoy_df$max[i] <- max(temp_dat$time_standard)
}

buoy_df$range <- buoy_df$max - buoy_df$min

write.csv(buoy_df, "buoy_summary.csv")

tdf <- buoy_df[which(buoy_df$range < 1150),]   #1119.791667 was the longest recorded data between two islands in our dataset. Suggests that any buoys that went longer doubled back on itself after a while somehow.

summary(tdf$range)   ##75 days on average spent going past the study islands.

length(tdf$id[which(tdf$range > 100)])   #203 buoys travel for longer than 100 days.


## Mantel test between distance and dispersal time ###########
dispersal_dat <- read.csv("matrix_min_island.csv")
distance_dat <- read.csv("matrix_dist_summary.csv")

dispersal_dat$X <- NULL
distance_dat$X <- NULL

disp.dist <- dist(dispersal_dat)
dist.dist <- dist(distance_dat)

library(vegan)
mantel(disp.dist, dist.dist, permutations = 9999, na.rm=TRUE, method = "spearman")
#p < 0.001. There is a significant correlation between dispersal times and distance.
# r statistic = 0.6195

prob_dat <- read.csv("prob_matrix_island.csv")
prob_dat$X <- NULL

prob.dist <- dist(prob_dat)
mantel(dist.dist, prob.dist, permutations = 9999, na.rm=TRUE, method = "spearman")

#p = 0.6113
#r = -0.0148

################################################################################
########## PLANTS #############################################################
################################################################################

## Getting plant distributions
fruit_dat_temp <- read.csv("fruit_list.csv")
fruit_dat <- fruit_dat_temp[which(fruit_dat_temp$Species != "Coccoloba uvifera"),]

list_plants <- unique(fruit_dat$Species)

library(spocc)

for(i in c(1:length(list_plants))){
  results <- occ(query = as.character(list_plants[i]), from= c('gbif'), limit = 1000)
  res <- occ2df(results)
  
  print(paste0(list_plants[i]))
  print(i)
  
  if(length(res$name > 0)){ 
    tryCatch({
      res$longitude <- as.numeric(res$longitude)
      res$latitude <- as.numeric(res$latitude)
      
      res$time <- as.POSIXct(strptime(res$date, "%Y-%m-%d"))
      res$year <- format(strptime(res$date, "%Y-%m-%d"), '%Y')
      res$year <- as.numeric(res$year)
      
      #subsetting to the Caribbean area
      temp <- res[which(res$longitude > -92.8929 & res$longitude < -59.2015 & res$latitude < 27.303 & res$latitude > 8.2679),]
      
      temp2 <- temp %>% distinct(longitude, latitude, date, .keep_all = TRUE) #removing duplicates
      
      if(length(temp2$name > 0)){
        write.csv(temp2, paste0("plant_point_",list_plants[i], ".csv"))
        
        point <- st_as_sf(x = temp2, 
                          coords = c("longitude", "latitude"),
                          crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
        st_write(point, paste0("point_",list_plants[i], ".shp"))
      } else{
        print(paste0(list_plants[i], " not enough points"))
      }
    }, error=function(e){})} else{
      print(paste0(list_plants[i], " not enough points"))
    }}

## Comparing artificial and field trials

paired_dat <- read.csv("paired_dat.csv")
bartlett.test(days ~ treatment, data = paired_dat)
shapiro.test(paired_dat$days)
paired_dat$trans_days <- log(paired_dat$days + 0.1)
shapiro.test(paired_dat$trans_days)

t.test(formula = days ~ treatment, data = paired_dat,
       mu = 0,
       paired = TRUE,
       var.equal = TRUE,
       conf.level = 0.95)

## t=0.047, df=4, p=0.9646

## Survival analyses #######################3
library(survival)
library(ggsurvfit)

dat_temp <- read.csv("survival_dat.csv")
dat <- dat_temp[which(dat_temp$status != "checked"),]

dat$Start <- as.POSIXct(strptime(dat$Start, "%m/%d/%Y"))
dat$sank <- as.POSIXct(strptime(dat$sank, "%m/%d/%Y"))
dat$time <- dat$sank - dat$Start
dat$time <- as.numeric(dat$time)
dat$time <- dat$time/86400   #converting from seconds to days

survfit2(Surv(time) ~ Species, data = dat) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + xlim(0,100)+
  add_confidence_interval()

survdiff(Surv(time) ~ Species , data = dat) #p < 0.001

## Alternative figure

dat.perc <- read.csv("perc_floating_dat.csv")

survdiff(Surv(Days) ~ Species , data = dat.perc) #p < 0.001

dat.perc2 <- dat.perc[which(dat.perc$Species != "Coccothrinax barbadensis" & dat.perc$Species != "Coccothrinax spissa" & dat.perc$Species != "Gaussia attenuata" & dat.perc$Species != "Syagrus amara"),]

dat.perc2$log_Day <- log(dat.perc2$Days)

dat.perc3 <- dat.perc2[which(dat.perc2$Type >0),]
dat.perc3$Type <- as.factor(dat.perc3$Type)

ggplot(data = dat.perc, aes(y = perc, x = Days, fill = Species, color = Species)) +
  geom_line(size = 1) +
  geom_point(shape=21, color="black", size=3) +
  labs(title = "",
       y = "Proportion sank (%)", x = "Time (days)") + 
  theme_bw()

## Figure 2
ggplot(data = dat.perc2, aes(y = perc, x = log_Day, fill = Species, color = Species)) +
  geom_line(size = 1, aes(linetype=Species)) +
  #geom_point(shape=21, color="black", size=3) +
  scale_shape_manual(values=seq(0,15))+
  labs(title = "",
       y = "Percent floating", x = "Days (log)") + 
  theme_bw()

ggplot(data = dat.perc2, aes(y = perc, x = log_Day, fill = Species, color=Species)) +
  geom_line(size = 1) +
  #scale_color_gradient(low="blue", high="red")+
  labs(title = "",
       y = "Percent floating", x = "Days (log)") + 
  theme_bw()


cc <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=10))

ggplot(data = dat.perc3, aes(y = perc, x = log_Day, fill = Type, color=Type)) +
  geom_line(size = 1) +
  scale_color_manual(values=cc) +
  labs(title = "",
       y = "Percent floating", x = "Days (log)") + 
  theme_bw()

### ESTIMATING WHAT PROPORTION OF PLANTS WILL BE DISPERSED BETWEEN ISLANDS BY SPECIES   (KEEPING DISPERSAL PROBABILITY CONSTANT)

dat.perc_temp <- read.csv("perc_floating_dat.csv")
matrix_min_temp <- read.csv("matrix_min_island.csv")
matrix_max_temp <- read.csv("matrix_max_island.csv")
prob_dat_temp <- read.csv("prob_matrix_island.csv")

islands <- matrix_min_temp$X

matrix_min_temp$X <- NULL
matrix_min <- as.matrix(matrix_min_temp)

matrix_max_temp$X <- NULL
matrix_max <- as.matrix(matrix_max_temp)

prob_dat_temp$X <- NULL
prob_dat <- as.matrix(prob_dat)

rownames(matrix_min) <- islands
colnames(matrix_min) <- islands

rownames(matrix_max) <- islands
colnames(matrix_max) <- islands

rownames(prob_dat) <- islands
colnames(prob_dat) <- islands

#c barb, c. spissa, g. attenuata, and s. amara did not float more than one day, so those are removed.
dat.perc <- dat.perc_temp[which(dat.perc_temp$Species != "Coccothrinax barbadensis" & dat.perc_temp$Species != "Coccothrinax spissa" & dat.perc_temp$Species != "Gaussia attenuata" & dat.perc_temp$Species != "Syagrus amara"),]
dat.perc$prob <- dat.perc$perc/100

plants <- unique(dat.perc$Species)

for(l in 1:length(plants)){
  
  temp_plant <- dat.perc[which(dat.perc$Species == plants[l]),]              # Getting species of interest
  temp_matrix <- matrix(-1, nrow= length(islands), ncol= length(islands), dimnames = list(c(islands), c(islands)))    # Making empty matrix
  
  for(i in 1:length(islands)){
    
    for(k in 1:length(islands)){
      
      if(length(matrix_min[k,i][!is.na(matrix_min[k,i])]) > 0 ) {              #For cases where the minimum transit time is NOT NA, proceed
        
        if(max(temp_plant$Days) > matrix_max[k,i]){                            #There may be cases where maximum viability time of a plant is greater than the maximum dispersal time between islands.
          
          vector1 = c() 
          
          for(m in ceiling(matrix_min[k,i]):ceiling(matrix_max[k,i])) {                  #In this case, consider the times between the minimum and maximum dispersal times.
            vector1 = c(vector1, temp_plant$prob[which(temp_plant$Days == m)])           #Using that time, consider the proportion of plants that remained viable within this interval and put it into a vector.
            
          }
        } else {
          
          vector1 = c() 
          
          for(n in ceiling(matrix_min[k,i]):max(temp_plant$Days)) {            #There may be cases where maximum viability time of a plant is less than the maximum dispersal time between islands.
            vector1 = c(vector1, temp_plant$prob[which(temp_plant$Days == n)])           #In this case, consider the times between the minimum dispersal time and maximum viability times.
            
          } 
          
        }
        
        if(length(vector1) == 0){                                               # In cases where minimum transit time is NA, put NA value (this is just to double check that the previous steps worked)
          temp_matrix[k,i] <- NA
        } else {
          time_length <- length(vector1)                                       # Equation to determine the average connectivity of a plant   (MAY NEED CHANGING)
          temp_matrix[k,i] <- (sum(vector1)*prob_dat[k,i])/(time_length)
        }
      } else {
        
        temp_matrix[k,i] <- NA
        
      }
    }
    
  }
  
  write.csv(temp_matrix, paste0(plants[l],"_connectivity_matrix.csv", ".csv"))   
}



####################################################33
# fruit traits and floatability
fruit_dat_temp <- read.csv("fruit_list.csv")
fruit_dat <- fruit_dat_temp[which(fruit_dat_temp$Species != "Coccoloba uvifera"),]
fruit_dat$volume <- (4/3)*pi*fruit_dat$mean_width*fruit_dat$mean_diam^2
fruit_dat$density <- fruit_dat$mean_mass/fruit_dat$volume
#write.csv(fruit_dat, "fruit_summary.csv")

library(car)

shapiro.test(log(fruit_dat$volume))
shapiro.test(fruit_dat$density)
shapiro.test(fruit_dat$no_of_islands_observed)
shapiro.test(log(fruit_dat$no_of_islands_observed))   #most normal
shapiro.test(fruit_dat$no_of_islands_observed^0.5)
qqPlot(log(fruit_dat$no_of_islands_observed))

shapiro.test(fruit_dat$max_viable_day)
shapiro.test(log(fruit_dat$max_viable_day + 0.01))
shapiro.test(fruit_dat$max_viable_day^0.5)    # Normal
fruit_dat$sqrt_mvd <- fruit_dat$max_viable_day^0.5

shapiro.test(fruit_dat$X95th_viable_day)
shapiro.test(log(fruit_dat$X95th_viable_day + 0.1))
shapiro.test(fruit_dat$X95th_viable_day^0.5)
fruit_dat$sqrt_vd <- fruit_dat$X95th_viable_day^0.5
fruit_dat$log_volume <- log(fruit_dat$volume)

shapiro.test(fruit_dat$mean_diam)   # Normal

#volume
model.1c = lm(sqrt_mvd ~ log_volume, data= fruit_dat)
summary(model.1c)   #p=0.2128

#density
model.1d = lm(sqrt_mvd ~ density, data= fruit_dat)
summary(model.1d)   #p=0.1552

#distribution x dispersal
model.1e = lm(log(fruit_dat$no_of_islands_observed) ~ fruit_dat$sqrt_mvd)
summary(model.1e)    #F=3.274, df=12, p=0.0955

model.1f = lm(log(fruit_dat$no_of_islands_observed) ~ fruit_dat$sqrt_vd)
summary(model.1f)   #F=3.094, df=12, p=0.104

#### ocean current dispersal and species distributions

float_fruit_dat <- fruit_dat[which(fruit_dat$X95th_viable_day >0),]
shapiro.test(float_fruit_dat$no_of_islands_observed)
shapiro.test(log(float_fruit_dat$no_of_islands_observed))    #More Normal

float_fruit_dat$log_islands <- log(float_fruit_dat$no_of_islands_observed)

shapiro.test(float_fruit_dat$X95th_viable_day)
shapiro.test(float_fruit_dat$X95th_viable_day^0.5)
shapiro.test(log(float_fruit_dat$X95th_viable_day))    # Normal
float_fruit_dat$log_viable <- log(float_fruit_dat$X95th_viable_day)

shapiro.test(log(float_fruit_dat$max_viable_day))
float_fruit_dat$log_max_viable <- log(float_fruit_dat$max_viable_day)

fmodel.2 = lm(log(float_fruit_dat$no_of_islands_observed) ~ log(float_fruit_dat$X95th_viable_day))
summary(fmodel.2)                 #F=4.92, df=8, p= 0.0573, r2=0.3808 (Marginally significant)

temp <- resid(fmodel.2)
shapiro.test(temp)  #Normal Residuals (W=0.97, p=0.898)

library(car)
ncvTest(fmodel.2) #Residuals are homoscedastic. Chisquared = 0.49, df=1, p=0.48

plot(log(float_fruit_dat$no_of_islands_observed) ~ log(float_fruit_dat$viable_day), pch = 16)
abline(fmodel.2, col="blue", lwd=2)

#### FIGURE 3
library(ggplot2)
plot1 <- 
  ggplot(data = float_fruit_dat, aes(log_viable, log_islands))+
  geom_point(size=2.5)+
  geom_smooth(method="lm", color="red", fill="light blue")+
  theme_bw()+
  ylim(-1.7, 5)+
  xlab("95th percentile viability day (log)")+
  ylab("number of islands observed (log)")
plot1

fmodel.2a = lm(log(float_fruit_dat$no_of_islands_observed) ~ log(float_fruit_dat$max_viable_day))
summary(fmodel.2a)                 #p= 0.0279

temp <- resid(fmodel.2a)
shapiro.test(temp)  #Normal Residuals (W=0.88, p=0.1341)

library(car)
ncvTest(fmodel.2a) #Residuals are homoscedastic. Chisquared = 0.77, df=1, p=0.38

plot1a <- 
  ggplot(data = float_fruit_dat, aes(log_max_viable, log_islands))+
  geom_point(size=2.5)+
  geom_smooth(method="lm", color="red", fill="light blue")+
  theme_bw()+
  ylim(-1.7, 5)+
  xlab("max viability day (log)")+
  ylab("")
plot1a

library(gridExtra)
grid.arrange(plot1, plot1a, ncol=2)

# Large fruits vs small fruits

large <- fruit_dat[which(fruit_dat$mean_diam > 19.1),]   #gape size of Coccyzus pluvialis, the largest bird in the Caribbean dataset
small <- fruit_dat[which(fruit_dat$mean_diam < 19.1),]

shapiro.test(large$max_viable_day^0.5)
shapiro.test(large$X95th_viable_day^0.5)
shapiro.test(log(large$no_of_islands_observed))

shapiro.test(small$max_viable_day^0.5)
shapiro.test(small$X95th_viable_day^0.5)
shapiro.test(log(small$no_of_islands_observed))

fmodel.3a = lm(log(large$no_of_islands_observed) ~ large$sqrt_mvd)
summary(fmodel.3a)   

fmodel.3b = lm(log(large$no_of_islands_observed) ~ large$sqrt_vd)
summary(fmodel.3b)   

########################################################################################################
## Determining if connectivity is correlated with distribution and distance ######## (Do for all floating species)
### Island/buoy background
island_buoy_sum <- read.csv("island_source_count.csv")
lack_islands <- island_buoy_sum$island[which(island_buoy_sum$count <= 10)]
lack_islands <- sub(" ", ".", lack_islands)

dispersal_dat <- read.csv(paste0(plants_list[4],"_connectivity_matrix.csv"))
distribution_dat <- read.csv(paste0(plants_list[4],"_distribution_matrix.csv"))
distance_dat <- read.csv("matrix_dist_summary.csv")

dispersal_dat$X <- NULL
distribution_dat$X <- NULL
distance_dat$X <- NULL

disp.mat_temp <- as.matrix(dispersal_dat)
distr.mat_temp <- as.matrix(distribution_dat)
dist.mat_temp <- as.matrix(distance_dat)

islands <- colnames(disp.mat_temp)

rownames(disp.mat_temp) <- islands
rownames(distr.mat_temp) <- islands
rownames(dist.mat_temp) <- islands
colnames(disp.mat_temp) <- islands
colnames(distr.mat_temp) <- islands
colnames(dist.mat_temp) <- islands

disp.mat <- disp.mat_temp[!rownames(disp.mat_temp)%in% lack_islands, !colnames(disp.mat_temp) %in% lack_islands]
distr.mat <- distr.mat_temp[!rownames(distr.mat_temp)%in% lack_islands, !colnames(distr.mat_temp) %in% lack_islands]
dist.mat <- dist.mat_temp[!rownames(dist.mat_temp)%in% lack_islands, !colnames(dist.mat_temp) %in% lack_islands]

disp.dist <- dist(disp.mat)
distr.dist <- dist(distr.mat)
dist.dist <- dist(dist.mat)

library(vegan)
mod2a <- mantel(disp.dist, distr.dist, permutations = 9999, na.rm=TRUE, method = "spearman")      # r=0.1603, p=0.016
mod2a
mod2b <- mantel(dist.dist, distr.dist, permutations = 9999, na.rm=TRUE, method = "spearman")      # r=0.06879, p=0.0783
mod2b

##############
###
## Determining if community similarity of floating/non_floating species is correlated with dispersal probability and minimum dispersal time

### Similarity matrix
library(stringdist)
library(dplyr)
library(betapart)
library(vegan)

dist_temp <- read.csv("island_distributions_NA.csv")
islands <- dist_temp$island

fruit_dat_temp <- read.csv("fruit_list.csv")
fruit_dat <- fruit_dat_temp[which(fruit_dat_temp$Species != "Coccoloba uvifera"),]
float_species <- fruit_dat$Species[which(fruit_dat$max_viable_day > 7.1)]
sink_species <- fruit_dat$Species[which(fruit_dat$max_viable_day <= 7)]
species <- fruit_dat$Species

dist <- dist_temp
colnames(dist) <- c("island", species)

float_matrix_temp <- matrix(-1, nrow= length(islands), ncol= length(islands), dimnames = list(c(islands), c(islands)))   #Creating empty matrices
not_float_matrix_temp <- matrix(-1, nrow= length(islands), ncol= length(islands), dimnames = list(c(islands), c(islands))) 

float_group <- dist %>% select(float_species)
not_float_group <- dist %>% select(sink_species)

float_group[is.na(float_group)] <- 0
not_float_group[is.na(not_float_group)] <- 0

rownames(float_group) <- islands
rownames(not_float_group) <- islands

dist1<-beta.pair(float_group, index.family="sorensen")
float_index <- dist1[[1]]

dist2<-beta.pair(not_float_group, index.family="sorensen")
not_float_index <- dist2[[1]]

#write.csv(as.data.frame(as.matrix(float_index)), "similarity_float.csv")
#write.csv(as.data.frame(as.matrix(not_float_index)), "similarity_not_float.csv")

float.dat <- read.csv("similarity_float.csv")
not.dat <- read.csv("similarity_not_float.csv")
distance_dat <- read.csv("matrix_dist_summary.csv")
prob_dat <- read.csv("prob_matrix_island.csv")
time_dat <- read.csv("matrix_min_island.csv")

float.dat$X <- NULL
not.dat$X <- NULL
distance_dat$X <- NULL
prob_dat$X <- NULL
time_dat$X <- NULL

float.mat_temp <- as.matrix(float.dat)
not.mat_temp <- as.matrix(not.dat)
dist.mat_temp <- as.matrix(distance_dat)
prob.mat_temp <- as.matrix(prob_dat)
time.dat_temp <- as.matrix(time_dat)

islands <- colnames(float.mat_temp)

rownames(float.mat_temp) <- islands
rownames(not.mat_temp) <- islands
rownames(dist.mat_temp) <- islands
rownames(prob.mat_temp) <- islands
rownames(time.dat_temp) <- islands
colnames(float.mat_temp) <- islands
colnames(not.mat_temp) <- islands
colnames(dist.mat_temp) <- islands
colnames(prob.mat_temp) <- islands
colnames(time.dat_temp) <- islands

island_buoy_sum <- read.csv("island_source_count.csv")
lack_islands <- island_buoy_sum$island[which(island_buoy_sum$count <= 10)]
lack_islands <- sub(" ", "_", lack_islands)

float.mat <- float.mat_temp[!rownames(float.mat_temp)%in% lack_islands, !colnames(float.mat_temp) %in% lack_islands]
not.mat <- not.mat_temp[!rownames(not.mat_temp)%in% lack_islands, !colnames(not.mat_temp) %in% lack_islands]
dist.mat <- dist.mat_temp[!rownames(dist.mat_temp)%in% lack_islands, !colnames(dist.mat_temp) %in% lack_islands]
prob.mat <- prob.mat_temp[!rownames(prob.mat_temp)%in% lack_islands, !colnames(prob.mat_temp) %in% lack_islands]
time.mat <- time.dat_temp[!rownames(time.dat_temp)%in% lack_islands, !colnames(time.dat_temp) %in% lack_islands]

float.dist <- dist(float.mat)
not.dist <- dist(not.mat)
dist.dist <- dist(dist.mat)
prob.dist <- dist(prob.mat)
time.dist <- dist(time.mat)

library(vegan)
mod2a <- mantel(prob.dist, float.dist, permutations = 9999, na.rm=TRUE, method = "spearman")
mod2a    #r=0.1488, p=0.047
mod2b <- mantel(prob.dist, not.dist, permutations = 9999, na.rm=TRUE, method = "spearman")
mod2b   #r=0.2601, p=0.0818

m2a <- mod2a$perm[which(mod2a$perm < mod2a$statistic)]
m2b <- mod2b$perm[which(mod2b$perm < mod2b$statistic)]

wilcox.test(m2a, m2b)   #W=47644558, p<0.001

mod3a <- mantel(time.dist, float.dist, permutations = 9999, na.rm=TRUE, method = "spearman")
mod3a    #r=-0.05299, p=0.7277
mod3b <- mantel(time.dist, not.dist, permutations = 9999, na.rm=TRUE, method = "spearman")
mod3b    #r=-0.1355, p=0.2087

m3a <- mod3a$perm[which(mod3a$perm < mod3a$statistic)]
m3b <- mod3b$perm[which(mod3b$perm < mod3b$statistic)]

wilcox.test(m3a, m3b)   #W=8928246, p<0.001

mod4a <- mantel(dist.dist, float.dist, permutations = 9999, na.rm=TRUE, method = "spearman")
mod4a    #r=-0.02128, p=0.6165
mod4b <- mantel(dist.dist, not.dist, permutations = 9999, na.rm=TRUE, method = "spearman")
mod4b    #r=0.05787, p=0.287

m4a <- mod3a$perm[which(mod4a$perm < mod4a$statistic)]
m4b <- mod3b$perm[which(mod4b$perm < mod4b$statistic)]

wilcox.test(m4a, m4b)   #W=14666883, p<0.001

mod5a <- mantel(prob.dist, time.dist, permutations = 9999, na.rm=TRUE, method = "spearman")
mod5a