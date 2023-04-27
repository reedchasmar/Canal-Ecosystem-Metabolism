##################################################################################################################################################

#####----- 

homedale <- read.csv(file = "C:/Users/reedc/Downloads/HOMEDALE 1 SE, ID US.csv", header = TRUE)
homedale$'Date' <- as.Date(homedale$'Date')
Iverson$'Date' <- as.Date(Iverson$'MDT', tz = "America/Denver")


#####-----

#- (1) -#
o <- split(Iverson, Iverson$'Date')

#- (2) -#
w <- list()
p <- list()
for (i in 1:175)
{
  a <- o[[i]]
  w[[i]] <- max(a[,2])
  p[[i]] <- min(a[,2])
}


#####----- 

homedale <- data.frame(homedale, unlist(w), unlist(p))[,-c(2)]


#####-----

x <- c(homedale[,1], rev(homedale[,1]))
y <- c((homedale[,3] - 32) * 5/9, rev((homedale[,2] - 32) * 5/9))
y2 <- c(homedale[,8], rev(homedale[,7]))


#####-----

plot(x, y, type = 'l', las = 1)
polygon(x, y, col = "#65BFFF")
polygon(x, y2, col = "plum2")