#####----- LOAD LIBRARIES AND REMOTES -----#####

library(libdplyr)
library(lubridate)
library(padr)
library(imputeTS)
library(geiger)
library(sp)
library(rgeos)
library(remotes)
library(readr)
library(streamMetabolizer)
library(rstan)
library(mice)
library(forecast)

remotes::install_github('appling/unitted')
remotes::install_github("USGS-R/streamMetabolizer")


#####----- (1) READ DATA FROM MINIDOT SENSORS; (2) COMPUTE KELVINS; (3) RENAME AND REORDER VARIABLES; AND (4) CONVERT TIME TO POSIX -----#####

#- (1) -#
minidot1 <- read.csv("C:\\Users\\reedc\\OneDrive\\Documents\\BSU\\Thesis\\miniDOT\\855208\\7450-855208\\Cat.txt", header=TRUE, check.names=FALSE)
minidot2 <- read.csv("C:\\Users\\reedc\\OneDrive\\Documents\\BSU\\Thesis\\miniDOT\\862018\\7450-862018\\Cat.txt", header=TRUE, check.names=FALSE)
rownames(minidot1) = seq(length = nrow(minidot1))
rownames(minidot2) = seq(length = nrow(minidot2))
minidot1 <- minidot1[ , -c(4, 7, 8)]
minidot2 <- minidot2[ , -c(4, 7, 8)]

#- (2) -#
minidot1$'Kelvin' <- as.numeric(minidot1$'Temperature') + 273.15
minidot2$'Kelvin' <- as.numeric(minidot2$'Temperature') + 273.15

#- (3) -#
colnames(minidot1) <- c("Unix Time", "UTC", "MDT", "Temp (°C)", "DO (mg/L)", "Temp (K)")
colnames(minidot2) <- c("Unix Time", "UTC", "MDT", "Temp (°C)", "DO (mg/L)", "Temp (K)")
minidot1 <- minidot1[, c(1, 2, 3, 4, 6, 5)]
minidot2 <- minidot2[, c(1, 2, 3, 4, 6, 5)]

#- (4) -#
minidot1[['MDT']] <- as.POSIXct(minidot1[['MDT']],format = "%Y-%m-%d %H:%M:%S")
minidot2[['MDT']] <- as.POSIXct(minidot2[['MDT']],format = "%Y-%m-%d %H:%M:%S")
minidot1[['UTC']] <- as.POSIXct(minidot1[['UTC']],format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
minidot2[['UTC']] <- as.POSIXct(minidot2[['UTC']],format = "%Y-%m-%d %H:%M:%S", tz = "UTC")


#####-----(1) READ IN CAMPBELL LOGGER DATA; (2) SET NIGHTTIME PAR TO ZERO; AND (3) REMOVE PAR WHEN SUN WAS BLOCKED BY BRIDGE RAIL-----#####

#- (1) -#
campbell <- read.csv("C:\\Campbellsci\\LoggerNet\\data\\CR1000_Table2.dat", header = TRUE)
campbell <- campbell[ , -c(2,3,5)]
campbell <- campbell[-c(12967), ]
campbell[ , 4] <- campbell[ , 4] * 12
colnames(campbell) <- c("time", "Cond (mS/cm)", "Temp (C)", "Depth (in)", "PAR")
campbell$Unix <- as.numeric(strptime(campbell[, 1], "%Y-%m-%d %H:%M:%S"))
campbell[['time']] <- as.POSIXct(campbell[['time']],format = "%Y-%m-%d %H:%M:%S")

#- (2) -#
campbell[,-c(1:4)][campbell[, -c(1:4)] < 10] <- 0
PAR_zeros <- filter(campbell, PAR == 0)

#- (3) -#
campbell$"%H:%M" <- format(campbell[,1], format = "%H:%M")
sunblock <- which(campbell$"%H:%M" %in% c("14:10", "14:20", "14:30", "14:40", "14:50", "15:00") == TRUE)
campbell[sunblock, 5] <- NA


#####----- (1) PAD DATASETS (1-MIN); (2) INTERPOLATE (STINE) NA; (3) CROP DATE-TIMES TO ALIGN DATA; (4) RENAME ROWS; AND (5) REDO MDT -----##### 

#- (1) -#
campbell_pad <- pad(campbell, interval = "min")
dotpad_1 <- (pad(minidot1, interval = "min", by = 'UTC'))
dotpad_2 <- (pad(minidot2, interval = "min", by = 'UTC'))

#- (2) -#
dotpad_1 <- na_interpolation(dotpad_1, option = "stine")
dotpad_2 <- na_interpolation(dotpad_2, option = "stine")
campbell_pad <- na_interpolation(campbell_pad, option = "stine")

#- (3) -#
dotpad_1 <- dotpad_1[-c(1:10293, 260458:260464), ]
dotpad_2 <- dotpad_2[-c(1:10125), ]
campbell_pad <- campbell_pad[-c(1:5967, 256132:257761), ]

#- (4) -#
rownames(dotpad_1) = seq(length=nrow(dotpad_1))
rownames(dotpad_2) = seq(length=nrow(dotpad_2))
rownames(campbell_pad) = seq(length=nrow(campbell_pad))

#- (5) -#
dotpad_1$MDT <- seq(as.POSIXct("2022-04-29 15:57:00", tz = "America/Denver"), as.POSIXct("2022-10-20 9:20:00", tz = "America/Denver"), 60)
dotpad_2$MDT <- seq(as.POSIXct("2022-04-29 15:57:00", tz = "America/Denver"), as.POSIXct("2022-10-20 9:20:00", tz = "America/Denver"), 60)

#####----- CONVERT FLUME DEPTH TO NOTUS BRIDGE THALWEG DEPTH WITH A STAGE-STAGE POWER FUNCTION (COEFFICIENTS IN EXCEL) -----#####

campbell_pad$'Notus Depth (in)' <- 33.75 * (campbell_pad[ , 4]) ^ 0.2088


#####----- REARRANGE AND RENAME 'CAMPBELL_PAD' VARIABLES -----#####

campbell_pad <- campbell_pad[ , c(6, 1, 2, 3, 4, 8, 5)]
colnames(campbell_pad) <- c("Unix", "MDT", "Cond (mS/cm)", "Temp (°C)", "Flume Depth (in)", "Notus Depth (in)", "PAR")


#####----- CREATE FUNCTION TO COMPUTE 100% SATURATION OXYGEN CONCENTRATION @ NONSTANDARD PRESSURE FROM TEMPERATURE (K) -----#####

#     C = equilibrium oxygen concentration at standard pressure
#     P = nonstandard pressure in atm
#   Pwv = partial pressure of water vapor in atm
# theta = related to the second virial coefficient of O2 (representing O2 as a real rather than an ideal gas)
#    Cp = 100% saturation oxygen concentration at nonstandard pressure

F <- function(x) {
  C <- exp(-139.34411 + (1.575701*10^5/x) - (6.642308*10^7/x^2) + (1.2438*10^10/x^3) - (8.621949*10^11/x^4) - (0.13*((1.7674*10^-2)
                                                                                                                     - (1.0754*10/x) + (2.1407*10^3/x^2))))
  P <- exp(5.25*log(1-713/44300))
  Pwv <- exp(11.8571 - (3840.7/x) - (216961/x^2))
  theta <- (9.672*10^-3) - (4.942*10^-5*x) + (6.436*10^-8*x^2)
  Cp <- C * (P*((1-Pwv/P)/(1-Pwv))) * ((1-theta*P)/(1-theta))
  return(Cp)
}


#####----- (1) APPLY FUNCTION TO COMPUTE OXYGEN SATURATION AT EACH TIMESTEP; AND (2) CALCULATE THE OXYGEN DEFICIT -----#####

#- (1) -#
dotpad_1 <- data.frame(dotpad_1, unlist(lapply(dotpad_1[, 5], F)))
dotpad_2 <- data.frame(dotpad_2, unlist(lapply(dotpad_2[, 5], F)))
colnames(dotpad_1) <- c("Unix", "UTC", "MDT", "Temp (°C)", "Temp (K)", "DO (mg/L)", "DO Sat (mg/L)")
colnames(dotpad_2) <- c("Unix", "UTC", "MDT", "Temp (°C)", "Temp (K)", "DO (mg/L)", "DO Sat (mg/L)")

#- (2) -#
dotpad_1[, "DO Deficit (mg/L)"] <- as.numeric(dotpad_1[, 7]) - as.numeric(dotpad_1[, 6])
dotpad_2[, "DO Deficit (mg/L)"] <- as.numeric(dotpad_2[, 7]) - as.numeric(dotpad_2[, 6])


#####----- CREATE VECTOR OF CHANNEL CROSS-SECTION DEPTHS AT EACH FOOT -----#####

x <- seq(0, 47, 1)
y <- c(0, 0.5, 3.5, 9, 15, 23, 31, 46, 58, 72, 76, 81, 84, 85, 85, 84, 84, 84, 84, 84, 83, 83, 82.5, 82, 82, 81, 81, 80, 79.5, 78.5, 78, 77, 77, 	75, 74.5, 74, 73, 70.5, 66, 61.5, 53, 43, 31, 22, 16, 7, 3, 0)


#####----- (1) USE ANGLE TO ACCOUNT FOR BANK HEIGHT DIFFERENCE; (2) LEVEL OUT DEPTHS OF CROSS SECTION; (3) CONVERT TO INCHES -----#####

# (1)
rad <- 0.03572
adj <- vector()
for (i in 1:48)
{
  adj[i] <- sin(rad) * ((i - 1) * 12)
}

# (2)
y_adj <- (y + adj) * -1

# (3)
x_in <- x * 12


#####----- CREATE MATRIX AND VECTORS TO STORE XS-AREA AND INTERSECTION POINTS BETWEEN WATER LEVEL AND CANAL CROSS-SECTION PERIMETER -----#####

XS <- vector()
x_pts <- matrix(0, nrow(campbell_pad), 2)
y_pts <- vector()


#####----- CREATE A FOR LOOP THAT CALCULATES AND STORES INTERSECTION POINTS AND XS-AREA FOR EACH DEPTH MEASUREMENT (1-MIN TIMESTEP) -----#####

for (i in 1:nrow(campbell_pad))
{
  
  
  #####----- (1) STORE WATER LEVEL AS DEPTH BELOW LEFT BANK; (2) FIND POINTS WHERE WATER LEVEL IS ABOVE CANAL CROSS-SECTION PERIMETER -----#####
  
  #- (1) -#
  w_lvl <- rep((-92 + campbell_pad[i, 6]), 48)
  
  #- (2) -#
  above <- w_lvl > y_adj
  
  
  #####----- DETERMINE THE TRUE INDICES WHERE THE 1-UNIT LAGGED DIFFERENCE FOR 'ABOVE' LOGICAL OBJECT DOES NOT EQUAL ZERO -----#####
  
  intersect.points <- which(diff(above) != 0)
  
  
  #####----- FIND THE SLOPES OF LINE SEGMENTS FOR WATER LEVEL AND CANAL CROSS-SECTION PERIMETER -----#####
  
  w_lvl.slopes <- w_lvl[intersect.points + 1] - w_lvl[intersect.points]
  y_adj.slopes <- y_adj[intersect.points + 1] - y_adj[intersect.points]
  
  
  #####----- (1) FIND X&Y INTERSECTION POINTS BETWEEN THE LINE SEGMENTS; (2) STORE POINTS IN MATRIX AND VECTOR -----#####
  
  #- (1) -#
  x.points <- (intersect.points + ((y_adj[intersect.points] - w_lvl[intersect.points]) / (w_lvl.slopes - y_adj.slopes)) - 1) * 12
  y.points <- w_lvl[intersect.points] + (w_lvl.slopes * (x.points - intersect.points))
  
  #- (2) -#
  x_pts[i, 1] <- x.points[1]
  x_pts[i, 2] <- x.points[2]
  y_pts[i] <- y.points[1]
  
  
  #####----- (1) CALCULATE AND STORE CROSS-SECTIONAL AREA AT EACH TIME STEP; (2) CLOSE FOR LOOP -----#####
  
  XS[i] <- geiger:::.area.between.curves(x_in, y_adj, w_lvl, xrange = x.points) / 144
}


#####----- STORE (1) DISCHARGE, (2) XS AREA, (3) VELOCITY, (4) WETTED WIDTH, AND (5) MEAN DEPTH IN CAMPBELL FILE -----#####

#- (1) -#
campbell_pad$CFS <- 22 * sqrt(32.17405 * ((campbell_pad[ , 5] / 12) - 0.45)^3)
#- (2) -#
campbell_pad$XS_Area <- XS
#- (3) -#
campbell_pad$Velocity <- campbell_pad$CFS / campbell_pad$XS_Area
#- (4) -#
campbell_pad$Wetted_Width <- (x_pts[ , 2] - x_pts[ , 1]) / 12
#- (5) -#
campbell_pad$'Mean_Depth (ft)' <- campbell_pad$XS_Area / campbell_pad$Wetted_Width


#####----- CALCULATE THE 1-MIN LAGGED DIFFERENCE IN OXYGEN CONCENTRATION -----#####

dotpad_1$deltaDO <- c(NA, diff(dotpad_1[,6]))
dotpad_2$deltaDO <- c(NA, diff(dotpad_2[,6]))


#####----- (1) CONVERT VELOCITY (CM/S) AND MEAN DEPTH (CM); (2) COMPUTE REAERATION COEFFICIENT AT 20ºC; AND (3) ADJUST FOR TEMPERATURE -----#####

#- (1) -#
v <- (campbell_pad$Velocity * 30.48)
D <- (campbell_pad$'Mean_Depth (ft)' * 30.48)

#- (2) -#
Ko2_20C <- 50.8 * (v^0.67) * (D^-0.85)

#- (3) -#
dotpad_1$'Ko2 (1/d)' <- Ko2_20C * 1.024^(as.numeric(dotpad_1[, 4]) - 20)
dotpad_2$'Ko2 (1/d)' <- Ko2_20C * 1.024^(as.numeric(dotpad_2[, 4]) - 20)
dotpad_1$'Ko2 (1/min)' <- dotpad_1$'Ko2 (1/d)' / 1440
dotpad_2$'Ko2 (1/min)' <- dotpad_2$'Ko2 (1/d)' / 1440


#####----- (1) SPLIT CAMPBELL DATA BY DAY; AND (2) USE AREA BETWEEN CURVES FUNCTION TO CALCULATE MOLES OF PAR EACH DAY -----##### 

#- (1) -#
campbell_pad$'Date' <- as.Date(campbell_pad$'MDT', origin = '1970-01-01', tz = "America/Denver")
o <- split(campbell_pad, campbell_pad$'Date')

#- (2) -#
par_day <- list()
for (i in 1:175)
{
  v <- o[[i]]
  par_day[[i]] <- geiger:::.area.between.curves(seq(0, (nrow(v) * 60) - 60, 60), rep(0, nrow(v)), v[,7], xrange = c(0, nrow(v)*60))
}
PAR_day <- unlist(par_day)/1000000



################################################################################################################################################## PREPARE DATA TO: (1) RUN WITH STREAMMETABOLIZER PACKAGE MODELING TOOLS; AND (2) UPLOAD TO STREAMPULSE FOR CLEANING AND ANALYSIS               #
#################################################################################################################################################

#####----- CREATE DATAFRAMES WITH NECESSARY VARIABLES FOR STREAMMETABOLIZER MODELING (METRIC UNITS) -----#####

SMM_1 <- data.frame(dotpad_1$'UTC', dotpad_1$'DO (mg/L)', dotpad_1$'DO Sat (mg/L)', campbell_pad$'Mean_Depth (ft)' / 3.2808399, 
                    dotpad_1$'Temp (°C)', campbell_pad$'PAR')
colnames(SMM_1) <- c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light")

SMM_2 <- data.frame(dotpad_2$'UTC', dotpad_2$'DO (mg/L)', dotpad_2$'DO Sat (mg/L)', campbell_pad$'Mean_Depth (ft)' / 3.2808399, 
                    dotpad_2$'Temp (°C)', campbell_pad$'PAR')
colnames(SMM_2) <- c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light")


attr(Notus$'MDT', "tzone") <- "UTC"
SMM_3 <- data.frame(Notus$'MDT', Notus$'DO', Notus$'DOsat', Notus$'Depth', Notus$'WtempC', Notus$'PAR')
colnames(SMM_3) <- c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light")


#####----- (1) RUN A BAYES MODEL FOR EACH DISSOLVED OXYGEN SITE; AND (2) PLOT RESULTS-----#####

#- (1) -#
# mm_1 <- metab(specs(mm_name('bayes')), data=SMM_1)
# mm_2 <- metab(specs(mm_name('bayes')), data=SMM_2)
mm_3 <- metab(specs(mm_name('bayes')), data=SMM_3)

#- (2) -#
# plot_metab_preds(mm_1)
# plot_metab_preds(mm_2, y_lim = list(GPP = c(-20, 20), ER = c(-20, 20)))
plot_metab_preds(mm_3, y_lim = list(GPP = c(0, 20), ER = c(-20, 0)))


#####----- (1) ADD VARIABLES; (2) RENAME; (3) SEQUENCE DATASET TO REDUCE SIZE; AND (4) WRITE TO .CSV FILE FOR UPLOAD TO STREAMPULSE -----#####

#- (1) -#
SMM_1 <- data.frame(SMM_1, campbell_pad$'Cond (mS/cm)', campbell_pad$'CFS' / 35.3146667)
SMM_2 <- data.frame(SMM_2, campbell_pad$'Cond (mS/cm)', campbell_pad$'CFS' / 35.3146667)

#- (2) -#
colnames(SMM_1) <- c("Date-Time (UTC)", "DO (mg/L)", "Saturation DO (mg/L)", "Depth (m)", "Water Temperature (C)", "Light, PAR (umol/m2/s)", 
                     "Specific Conductivity (mS/cm)", "Discharge (cms)")
colnames(SMM_2) <- c("Date-Time (UTC)", "DO (mg/L)", "Saturation DO (mg/L)", "Depth (m)", "Water Temperature (C)", "Light, PAR (umol/m2/s)", 
                     "Specific Conductivity (mS/cm)", "Discharge (cms)")

#- (3) -#
# SMM_1.1 = SMM_1[seq(4, nrow(SMM_1), 7), ]
# SMM_2.1 = SMM_2[seq(4, nrow(SMM_2), 7), ]

#- (4) -#
# write_csv(SMM_1.1, file = "ID_Notus_2022-10-20_XX.csv")
# write_csv(SMM_2.1, file = "ID_Iverson_2022-10-20_XX.csv")


#################################################################################################################################################
# This code explains how to download data from StreamPulse after cleaning. Unzip download, e.g. all_sp_data.zip. Do not use utils::unzip, as it # # may truncate the resulting files (see details of ?unzip)                                                                                      #
#################################################################################################################################################

#####----- READ IN DATA -----#####

unzip_location = 'C:/Users/reedc/OneDrive/Documents/BSU/Thesis/Data/all_sp_data'

sp_files = list.files(unzip_location, pattern = 'all_sp_data[0-9]{2}\\.csv', full.names = TRUE)


#####----- ALL AT ONCE (IF YOU HAVE READR AND PURRR INSTALLED) -----#####

sp_data = purrr::map_dfr(sp_files, readr::read_csv)


#####----- ALL AT ONCE (USING BASE FUNCTIONS ONLY; SLOWER) -----#####

# sp_data = data.frame()

# for(f in sp_files) {sp_data = rbind(sp_data, read.csv(f))}


#####----- JUST FOR A SUBSET OF SITES/VARIABLES/ETC. (CONSERVES MEMORY) -----#####

regionIDs = c('ID')
sp_data_subset = purrr::map_dfr(sp_files, ~dplyr::filter(readr::read_csv(.), regionID %in% regionIDs))


#####----- (1) STORE IDAHO SUBSET; (2) SPLIT BY CANAL & VARIABLE; AND (3) RECOMBINE INTO SEPARATE AND (4) COMBINED NA-PAD DATAFRAMES -----#####

#- (1) -#
ID_sub <- as.data.frame(sp_data_subset)

#- (2) -#
a <- split(ID_sub, ID_sub$'siteID')
b <- split(a[[2]], a[[2]]$'variable')
c <- split(a[[3]], a[[3]]$'variable')
d <- split(a[[1]], a[[1]]$'variable')

#- (3) -#
Iverson <- data.frame(pad(b[[1]], interval = '15 min')[c(3, 5)], pad(b[[2]], interval = '15 min')[5], pad(b[[3]], interval = '15 min')[5], 	pad(b[[4]], interval = '15 min')[5], pad(b[[5]], interval = '15 min')[5], pad(b[[6]], interval = '15 min')[5], pad(b[[7]], 
                                                                                                                                                                                                                                                               interval = '15 min')[5])
colnames(Iverson) <- c("MDT", "Depth", "Discharge", "DO", "PAR", "DOsat", 
                       "SpcfcCond", "WtempC")
Iverson <- Iverson[c(1, 8, 7, 2, 4, 6, 5, 3)]
attr(Iverson$'MDT', "tzone") <- "America/Denver"

Notus <- data.frame(pad(c[[1]], interval = '15 min')[c(3, 5)], pad(c[[2]], interval = '15 min')[5], pad(c[[3]], interval = '15 min')[5], 	pad(c[[4]], interval = '15 min')[5], pad(c[[5]], interval = '15 min')[5], pad(c[[6]], interval = '15 min')[5], pad(c[[7]], 
                                                                                                                                                                                                                                                             interval = '15 min')[5])
colnames(Notus) <- c("MDT", "Depth", "Discharge", "DO", "PAR", "DOsat", 
                     "SpcfcCond", "WtempC")
Notus <- Notus[c(1, 8, 7, 2, 4, 6, 5, 3)] 
attr(Notus$'MDT', "tzone") <- "America/Denver"


x <- pad(as.data.frame(d[1])[ , c(3, 5)], interval = "15 min")
y <- pad(as.data.frame(d[2])[ , c(3, 5)], interval = "15 min")
DCEW <- data.frame(x, y[ , 2])
rownames(DCEW) <- seq(1, nrow(DCEW))
colnames(DCEW) <- c("MDT", "DO", "DO%")
attr(DCEW$'MDT', "tzone") <- "America/Denver"
DCEW$'DOsat' <- DCEW$'DO' / (DCEW$'DO%' / 100)

#- (4) -#
IveNot <- data.frame(pad(c[[1]], interval = '15 min')[3], pad(c[[2]], interval = '15 min')[5], 	pad(c[[3]], interval = '15 min')[5], 
                     pad(c[[4]], interval = '15 min')[5], pad(c[[5]], interval = '15 min')[5], pad(c[[6]], interval = '15 min')[5], pad(b[[3]], 
                                                                                                                                        interval = '15 min')[5], pad(b[[5]], interval = '15 min')[5])
colnames(IveNot) <- c("MDT", "Discharge", "DO.N", "PAR", "DOsat.N", 
                      "SpcfcCond", "DO.I", "DOsat.I")
IveNot <- IveNot[c(1, 2, 4, 6, 3, 7, 5, 8)] 
attr(IveNot$'MDT', "tzone") <- "America/Denver"




x <- complete(mice(IveNot, method = 'polyreg'))

plot(x[11000:15000,1], x[11000:15000,5], type = 'l')
plot(x[11000:15000,1], x[11000:15000,6], type = 'l')