# Finding duplicates in occurrence data
# Lauren Jenkins
# Fall 2021

# look at duplicates by date
duplicates <- occs[duplicated(occs$date_collected) | duplicated(occs$date_collected, fromLast = TRUE), ]
# this should remove rows where date & data source are the same
# we don't actually want to remove these records, 
# just trying to isolate where date is same, but data source different
duplicates <- duplicates[!duplicated(duplicates[, c("date_collected", "datasource")]), ] 
# 4,549 obs

# now only keep duplicated dates to isolate them
duplicates <- duplicates[duplicated(duplicates$date_collected) | duplicated(duplicates$date_collected, fromLast = TRUE), ]

# only want duplicate dates where lat/long are same in ones digit
duplicates$lat_simp <- substr(duplicates$latitude, 1, 2)
duplicates$long_simp <- substr(duplicates$longitude, 1, 3)
duplicates$date_lat_long <- paste0(duplicates$lat_simp, ", ", 
                                   duplicates$long_simp, ", ", 
                                   duplicates$ date_collected) # 976 obs
duplicates <- duplicates %>% dplyr::group_by(date_lat_long) %>% filter(n() > 1) # 148 obs

# these appear to NOT be duplicates : 17 & 18, 131 & 132, 133 & 134, 135 & 136
false_duplicates <- c(17, 18, 131, 132, 133, 134, 135, 136)
# remove false duplicates
duplicates <- duplicates[-false_duplicates, ]