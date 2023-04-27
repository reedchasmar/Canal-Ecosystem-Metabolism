#####----- LOAD LIBRARIES -----#####

library(TrenchR)
library(lubridate)
library(StreamMetabolism)


#####----- CREATE DATAFRAME WITH VARIABLES FOR WEEKLY MULTIVARIABLE PLOTS -----#####

w <- data.frame(Notus$'MDT', Notus$'DO', Notus$'DOsat', Notus$'PAR')
# w <- data.frame(dotpad_1$'MDT', dotpad_1$'DO (mg/L)', dotpad_1$'DO Sat (mg/L)', campbell_pad$PAR)
colnames(w) <- c("MDT", "DO (mg/L)", "DO sat (mg/L)", "PAR")

m <- data.frame(Iverson$'MDT', Iverson$'DO', Iverson$'DOsat', Iverson$'PAR')
# m <- data.frame(dotpad_2$'MDT', dotpad_2$'DO (mg/L)', dotpad_2$'DO Sat (mg/L)', campbell_pad$PAR)
colnames(m) <- c("MDT", "DO (mg/L)", "DO sat (mg/L)", "PAR")

k <- DCEW[ , c(1, 2, 4)]
colnames(k) <- c("MDT", "DO (mg/L)", "DO sat (mg/L)")
start <- data.frame(as.POSIXct("2022-04-29 16:00:00"), "NA", "NA")
colnames(start) <- c("MDT", "DO (mg/L)", "DO sat (mg/L)")
k <- rbind(start, k[1:10726, ])
k <- pad(k, interval = "15 min")

#####----- (1) CREATE A DATE COLUMN; (2) CONVERT TO DOY; (3) USE DOY AND LONGITUDE TO DETERMINE 'SOLAR NOON'; AND (4) CONVERT TO POSIX -----#####

#- (1) -#
w$'Date' <- as.Date(w[ , 1], tz = "America/Denver")
#- (2) -#
w$DOY <- lubridate::yday(w[,5])
#- (3) -#
w$"Solar Noon" <- solar_noon(lon = -116.8, doy = w[ , 6], offset = -6)
#- (4) -#
w$"Solar Noon" <- as.POSIXct(paste(w[ , 5], paste(hms::hms(lubridate::seconds_to_period(floor(w[ ,7] * 60 * 60))), "MDT", sep = " "), sep = " "))

#- (1) -#
m$'Date' <- as.Date(m[ , 1], tz = "America/Denver")
#- (2) -#
m$DOY <- lubridate::yday(m[,5])
#- (3) -#
m$"Solar Noon" <- solar_noon(lon = -116.8, doy = m[ , 6], offset = -6)
#- (4) -#
m$"Solar Noon" <- as.POSIXct(paste(m[ , 5], paste(hms::hms(lubridate::seconds_to_period(floor(m[ , 7] * 60 *60))), "MDT", sep = " "), sep = " "))

#- (1) -#
k$'Date' <- as.Date(k[ , 1], tz = "America/Denver")
#- (2) -#
k$DOY <- lubridate::yday(k[ , 4])
#- (3) -#
k$"Solar Noon" <- solar_noon(lon = -116.171, doy = k[ , 5], offset = -6)
#- (4) -#
k$"Solar Noon" <- as.POSIXct(paste(k[ , 4], paste(hms::hms(lubridate::seconds_to_period(floor(k[ , 6] * 60 *60))), "MDT", sep = " "), sep = " "))


#####----- SPLIT DATAFRAME BY DAY -----#####

w2 <- split(w, w[ , 5])
m2 <- split(m, m[ , 5])
k2 <- split(k, k[ , 4])

#####----- CREATE A LIST AND FILL IT WITH 25 DATAFRAMES WITH 7-DAY PERIODS -----#####

b <- list()
for(i in 1:25)
{
  b[[i]] <- do.call(rbind, w2[(((i - 1) * 7) + 1):(((i - 1) * 7) + 7)])
}

b2 <- list()
for(i in 1:25)
{
  b2[[i]] <- do.call(rbind, m2[(((i - 1) * 7) + 1):(((i - 1) * 7) + 7)])
}

b3 <- list()
for(i in 1:25)
{
  b3[[i]] <- do.call(rbind, k2[(((i - 1) * 7) + 1):(((i - 1) * 7) + 7)])
}


#####----- CREATE A POSIX STRING OF DATE-TIMES WHEN OXYGEN SENSORS WERE CLEANED -----#####

visits <- as.POSIXct(c("2022/05/03 20:00:00 MDT", "2022/05/17 18:00:00 MDT", "2022/05/28 14:00:00 MDT", "2022/06/03 15:00:00 MDT", 
                       "2022/06/11 14:00:00 MDT", "2022/06/19 16:00:00 MDT", "2022/06/24 18:00:00 MDT", "2022/07/01 15:00:00 MDT", "2022/07/09 14:00:00 MDT", 
                       "2022/07/16 15:00:00 MDT", "2022/07/24 12:45:00 MDT", "2022/07/31 20:00:00 MDT", "2022/08/06 09:00:00 MDT", "2022/08/12 10:00:00 MDT", 
                       "2022/08/19 10:00:00 MDT", "2022/08/27 09:30:00 MDT", "2022/09/10 09:30:00 MDT", "2022/09/22 17:00:00 MDT", "2022/09/24 10:00:00 MDT", 
                       "2022/10/01 10:00:00 MDT", "2022/10/08 09:30:00 MDT", "2022/10/14 10:30:00 MDT"), format = "%Y/%m/%d %H:%M:%OS")


################################################################################################################################################## THE FOLLOWING IS A LOOP THAT CREATES MULTIVARIABLE PLOTS WITH 'PAR', 'SOLAR CYCLES', 'DO SATURATION', AND 'DO CONCENTRATION'                  #
#################################################################################################################################################

#####----- BEGIN A 'FOR LOOP' THAT CYCLES THROUGH EACH 7-DAY DATAFRAME -----##### 

for (i in 1:25)
{
  r <- b[[i]]
  r2 <- b2[[i]]
  r3 <- b3[[i]]
  
  #####----- INTIATE THE OUTPUT OF .PNG FILE(S) FROM FOLLOWING PLOT(S) -----#####
  
  # png(paste0("Week", "_", i, ".png"), width = 2750, height = 1750, pointsize = 35)
  png(paste0("Week", "_", i, ".png"), width = 1750, height = 1750, pointsize = 35)
  
  
  #####----- RUN 'FOR LOOP' THAT STORES 'SOLAR NOON' & 'MIDNIGHT' FOR EACH DAY IN A VECTOR AS POSIX -----#####
  
  f <- vector()
  for (i in 1:7)
  {
    f[i] <- split(r, r[, 5])[[i]][1, 7]
  }
  f <- as.POSIXct(f, origin = '1970-01-01')
  f <- c(f, f[1] - 43200, f + 43200)
  
  
  #####----- (1) RUN 'FOR LOOP' THAT STORES 'SUNRISE' & 'SUNSET' FOR EACH DAY IN A LIST; AND (2) UNLIST AND STORE AS POSIX -----#####
  
  #- (1) -#
  j <- list()
  for (i in 1:7)
  {
    h <- split(r, r[, 5])[[i]][1, 1]
    j[[i]] <- sunrise.set(43.74, -116.8, as.POSIXct(as.numeric(h), origin = '1970-01-01'), timezone = "America/Denver", 1)
  }
  
  #- (2) -#
  j2 <- as.POSIXct(unlist(j), origin = '1970-01-01')
  j2 <- c(as.POSIXct(as.character(as.Date(r[1, 1])), tz = "America/Denver"), j2, as.POSIXct(as.character(as.Date(r[nrow(r), 1])), 
                                                                                            tz = "America/Denver") + 86400)
  
  
  #####----- SET MARGINS FOR MULTIVARIABLE PLOT -----#####
  
  par(mar = c(5, 4, 4, 6))
  
  
  #####----- PLOT 'PAR' ON SECONDARY (RIGHT) AXIS -----#####
  
  plot(r[,1], r[,4], axes = FALSE, bty = "n", xlab = "", ylab = "", xlim = c(as.POSIXct(paste(as.character(r[1,5]), "00:00:00 MDT", sep = " ")), 
                                                                             as.POSIXct(paste(as.character(r[nrow(r),5]+1), "00:00:00 MDT", sep = " "))), ylim = c(0, 3000), type = 'l', xaxs = "i", 
       yaxs = "i", main = "")
  axis(side=4, las = 1)
  mtext("PAR", side=4, line=4)
  
  
  #####----- (1) PLOT POLYGONS FOR 'DAYTIME' AND 'NIGHTTIME'; AND (2) ADD VERTICAL LINES FOR 'SOLAR NOON', 'MIDNIGHT', AND 'VISITS' -----#####
  
  #- (1) -#
  x <- c(j2[1:2], rev(j2[1:2]), NA, j2[3:4], rev(j2[3:4]), NA, j2[5:6], rev(j2[5:6]), NA, j2[7:8], rev(j2[7:8]), NA, j2[9:10], rev(j2[9:10]), NA,
         j2[11:12], rev(j2[11:12]), NA, j2[13:14], rev(j2[13:14]), NA, j2[15:16], rev(j2[15:16]), NA)
  y <- rep(c(rep(0, 2), rep(3000, 2), NA), 8)
  polygon(x, y, col = 'royalblue1', border = NA)
  x <- c(j2[2:3], rev(j2[2:3]), NA, j2[4:5], rev(j2[4:5]), NA, j2[6:7], rev(j2[6:7]), NA, j2[8:9], rev(j2[8:9]), NA, j2[10:11], rev(j2[10:11]), NA,
         j2[12:13], rev(j2[12:13]), NA, j2[14:15], rev(j2[14:15]), NA)
  y <- rep(c(rep(0, 2), rep(3000, 2), NA), 7)
  polygon(x, y, col = 'lightgoldenrod1', border = NA)
  
  #- (2) -#
  abline(v = c(f), col = 'gray64', lty = 'dashed')
  # abline(v = visits, col = 'black', lwd = 2)
  
  
  #####----- FILL IN DAILY 'PAR' CURVES WITH POLYGON -----#####
  
  x <- c(r[,1], rev(r[,1]))
  y <- c(rep(0, nrow(r)), rev(r[,4]))
  polygon(x, y, col = "darkorange2")
  
  
  #####----- (1) OVERLAY NEW PLOT OF OXYGEN SATURATION ON PRIMARY (LEFT) AXIS; AND (2) ADD POINTS FOR SECOND SENSOR -----#####
  
  #- (1) -#
  par(new = TRUE)
  plot(r[,1], r[,3], xlim = c(as.POSIXct(paste(as.character(r[1,5]), "00:00:00 MDT", sep = " ")), as.POSIXct(paste(as.character(r[nrow(r),5]+1), 
                                                                                                                   "00:00:00 MDT", sep = " "))), ylim = c(0, 13.1), axes = FALSE, bty = "n", xlab = "", ylab = "", pch = 20, col = 'turquoise4', 
       xaxs = "i", yaxs = "i")
  
  #- (2) -#
  # points(r2[,1], r2[,3], pch = 20, col = 'orchid4')
  points(r3[,1], r3[,3], pch = 20, col = 'orchid4')
  
  
  #####----- (1) OVERLAY NEW PLOT OF OXYGEN CONCENTRATION ON PRIMARY AXIS; AND (2) ADD POINTS FOR SECOND SENSOR -----#####
  
  #- (1) -#
  par(new = TRUE)
  plot(r[,1], r[,2], xlim = c(as.POSIXct(paste(as.character(r[1,5]), "00:00:00 MDT", sep = " ")), as.POSIXct(paste(as.character(r[nrow(r),5]+1), 
                                                                                                                   "00:00:00 MDT", sep = " "))), ylim = c(0, 13.1), xlab = "Day", ylab = "DO (mg/L)", pch = 20, col = 'turquoise', las = 1, 
       xaxt = "n", xaxs = "i", yaxs = "i")
  axis.POSIXct(1, at=seq(as.POSIXct(paste(as.character(r[1,5]), "00:00:00 MDT", sep = " ")), as.POSIXct(paste(as.character(r[nrow(r),5]+1), 
                                                                                                              "00:00:00 MDT", sep = " ")), by="1 day", format="%m-%d"))
  
  #- (2) -#
  # points(r2[,1], r2[,2], pch = 20, col = 'orchid2')
  points(r3[,1], r3[,2], pch = 20, col = 'orchid2')
  
  
  #####----- SHUT DOWN CURRENT PLOTTING DEVICE AND CLOSE 'FOR LOOP' -----#####
  
  dev.off()
}