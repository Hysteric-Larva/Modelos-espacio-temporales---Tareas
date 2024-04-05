library(GeoModels)
library(mapproj)
library(moments)
library(spdep)
library(sp)
data(anomalies)
head(anomalies,5)
anomalies <- head(anomalies, 500)
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
      My_value <- My_value + sqrt((lat[i] - lat[j])^2 + (lon[i] - lon[j])^2)*abs(z[i]-z[j])
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
    for (j in 1:n) {
      for (k in 1:n) {
        for (l in 1:n) {
          var_aux1<-0
          norma_aux<-0
          norma_aux <- sqrt((lat[i] - lat[j])^2 + (lon[i] - lon[j])^2)*sqrt((lat[k] - lat[l])^2 + (lon[k] - lon[l])^2)
          if (i == j && k == l) {
            next  # If i=j and k=l, continue to the next iteration
          }
          if (i == k || i == l || j == k || j == l) {
            var_aux1 <- var_aux1 + (p * (1 - p) - 4 * p^2 * (1 - p))
          }
          if ((i == k && j == l) || (i == l && j == k)) {
            var_aux1 <- var_aux1 + 4*p^2 * (1 - p^2)
          }
          if ((i != k && j == l) || (i != l && j == k) || (i == k && j != l) || (i == l && j != k)) {
            var_aux1 <- var_aux1 + 4*p^3 * (1 - p)
          }
          norma = norma + norma_aux*var_aux1
        }
      }
    }
  }
  

  
  variance_value <- norma 
  return(variance_value)
}  
  #############################################################
  
get_myts<- function(lat,lon,z,p,tau){
  expected_value <- expected_value(lat, lon,p)
  variance_value <- Variance(lat,lon,p)
  print(str(variance_value))
  z_t <- Y_t(z,tau)
  Myt <- My(lat,lon,z_t)
  stMyt <- (Myt-expected_value)/sqrt(variance_value)
  return(stMyt)
  
}
A = get_myts(lat,lon,z,0.5,2)
#n=length(lat)

#W <- knearneigh(anomalies, k = 100)  # k nearest neighbors
#WT <- knn2nb(W)
# Calculate Geary's C
#geary_result <- geary.test(data$variable, listw = WT)

# Print the result
#print(geary_result)
print(str(A))