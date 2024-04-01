library(GeoModels)
library(mapproj)
library(moments)
data(anomalies)
head(anomalies,5)
loc=cbind(anomalies[,1],anomalies[,2])
z=cbind(anomalies[,3])
summary(anomalies)
lon = cbind(anomalies[,1])
lat = cbind(anomalies[,2])

### map data
quilt.plot(loc,z,xlab="long",ylab="lat")
map("usa", add = TRUE)

### histogram and scatterplot
summary(anomalies)
hist(z,main="Anomalies histogram")

hist(lon,main="Lon histogram")

hist(lat,main="Lat histogram")

boxplot(lon,main = "Longitud")
boxplot(lat,main = "Latitud")
boxplot(z,main = "z")
kurtosis(z)
plot(lon,z)
plot(lat,z)

################################### Parte 6 ######################
#####################        My        ############################

My <- function(lat,lon,z){
  My_value <- 0
  n <- length(lat)  # Assuming lat and lon are vectors of the same length
  
  # Iterate over the array using a for loop
  for (i in 1:n) {
    for (j in 1:n){
      My_value <- My_value + sqrt((lat[i] - lat[j])^2 + (lon[i] - lon[j])^2)*(z[i]-z[j])
    }
  }
  return(My_value)
}
Y_t <- function(z,tau){
  n <- length(z)
  z_t <- array(NA, dim = n)
  for (i in 1:n){
    if (z[i]>=tau){
      z_t[i]=1
    }
    else{
      z_t[i]=0
    }
  }
  return(z_t)
}


##########                 EXPECTED VALUE      ##################
p=0.5
expected_value <- function(lat, lon,p) {
  norma <- 0
  n <- length(lat)  # Assuming lat and lon are vectors of the same length
  
  # Iterate over the array using a for loop
  for (i in 1:n) {
    for (j in 1:n){
      norma <- norma + sqrt((lat[i] - lat[j])^2 + (lon[i] - lon[j])^2)
    }
  }
  
  expected_value <- norma * p * (1 - p) * 2
  return(expected_value)
}
############################ Variance Value    ####################
Variance <- function(lat, lon,p) {
  norma <- 0
  n <- length(lat)  # Assuming lat and lon are vectors of the same length
  
  # Iterate over the array using a for loop
  for (i in 1:n) {
    for (j in 1:n){
      norma <- norma + sqrt((lat[i] - lat[j])^2 + (lon[i] - lon[j])^2)
    }
  }
  
  variance_value <- norma * p * (1 - p) * 2*(1-2*p*(1-p))
  return(variance_value)
}  
  #############################################################
  
get_myts<- function(lat,lon,z,p,tau){
  expected_value <- expected_value(lat, lon,p)
  print(str(expected_value))
  variance_value <- Variance(lat,lon,p)
  print(str(variance_value))
  z_t <- Y_t(z,tau)
  print(str(z_t))
  Myt <- My(lat,lon,z_t)
  print(str(Myt))
  stMyt <- (Myt-expected_value)/sqrt(variance_value)
  return(stMyt)
  
}
get_myts(lat,lon,z,0.5,2)
