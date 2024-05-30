####             Libraries import         ####
library("GeoModels")
library("geoR")
library("geosphere")
library("spatstat")
library(FNN)
library(boot)
library(caret)
library(gstat)
library(sp)
library(nlme)
library(latex2exp)
library(stars)
library(ggplot2)
library(gridExtra)
###############################################
##  Data loading       #######################
data(anomalies)
loc <- cbind(anomalies[,1],anomalies[,2])
z <- cbind(anomalies[,3])
lon = cbind(anomalies[,1])# Invoco lon y lat por separado para usarlos de manera más sencilla.
lat = cbind(anomalies[,2])

########################  Part I, excercise a   #########################
model_l1 <- lm(anomalies[,3]~ anomalies[,1] + anomalies[,2] + I(anomalies[,1]^2) + I(anomalies[,2]^2) + I(anomalies[,1]*anomalies[,2]))



##################### Part I, excercise  b     ##########################
# For test propourses    #############
#k <- 10
#x_selected <- -85.250
#y_selected <- 31.570

kappa_t <- function(k,x_selected,y_selected){
  if (k>1){
    distances <- sqrt((anomalies[,1] - x_selected)^2 + (anomalies[,2] - y_selected)^2)#Calculation of distnaces
    nearest_neighbors <- order(distances)[1:k]#selection of 10 closest ones
    closest_points <- anomalies[nearest_neighbors, ] #selection of 10 closest data
    kappa = sum(closest_points[,3])/k    #kappa calculation
    return(kappa)
  }
  else {
    print("k must be greater than 1 !!!") #Just in case.
  }
}
#the following line applies the function to all the dataframe to calculate the lm.
anomalies_aux <- apply(anomalies[, c("lon", "lat")],1 ,function(row) { kappa_t(10,row["lon"],row["lat"])})
model_l2 <- lm(z~anomalies_aux)

############################### Part I, excercise c      #################################
anomalies_df <- as.data.frame(anomalies)
# Perform k-fold cross-validation
ctrl <- trainControl(method = "cv", number = 10)
#fit a regression model and use k-fold CV to evaluate performance
model_1 <- train(z ~ lon + lat + I(lon^2) + I(lat^2)+I(lat*lon), data = anomalies_df, method = "lm", trControl = ctrl)
####################################################################             
print(model_1)
anomalies_aux_aux <- anomalies
anomalies_aux_aux <- cbind(anomalies_aux_aux,anomalies_aux)
anomalies_aux_df <- as.data.frame(anomalies_aux_aux)
model_2 <- train(z ~ anomalies_aux, data = anomalies_aux_df, method = "lm", trControl = ctrl)
print(model_2)
######################### Part I, excercise d       #####################
boxplot(residuals(model_1))
boxplot(residuals(model_2))

data_mod <-SpatialPointsDataFrame(coords = anomalies_df[,c("lon","lat")], data = anomalies_df)
variogram_model <- variogram(z ~ 1, data = data_mod)
plot(variogram_model)
######## TEST
# Calculate the residuals
gstat_quad <- gstat(formula = z ~ lon + I(lon^2) + lat + I(lat^2) + lon*lat,                  locations = ~lon+lat,
                    data = anomalies_df)
gstat_knn <-  gstat(formula = z ~ anomalies_aux,locations = ~lon+lat, data = anomalies_df)
bu = c(0.1, 0.5, seq(1, 20, 1))
variogram_quad <- variogram(gstat_quad, boundaries = bu)
variogram_knn <- variogram(gstat_knn, boundaries = bu)
data_to_plot <- data.frame(Distance = variogram_quad$dist, Quad = variogram_quad$gamma, KNN = variogram_knn$gamma)
## Plotting the empirical variograms for both models
matplot(data_to_plot$Distance, data_to_plot[, c("Quad", "KNN")],type="o", pch = 16,col = c("black", "red"), ylim = range(c(data_to_plot$Quad, data_to_plot$KNN)),xlab = "Distancia", ylab = "Semivarianza", main = "Variogramas Obtenidos")

###########################################################
#        EXTRA
############################################################
P.sinusoidal <- mapproject(loc[,1],loc[,2],projection = "sinusoidal")
loc_aux <- cbind(P.sinusoidal$x,P.sinusoidal$y)*6371
maxdist <- max(dist(loc_aux))
evariog <- GeoVariogram(data=z,coordx = loc_aux, maxdist = maxdist/4) 
plot(evariog,ylim = c(0,1),pch = 20, xlab = "Km", ylab="Semi-variogram")
GeoCovariogram(residuals(model_1),show.vario=TRUE, vario=evariog,pch=20,ylim=c(0,1.3))




######################## Part II a  ##################
I = Inf
lower = list(mean = -I,sill=0,nugget=0,scale=0)
upper = list(mean = I,sill= I,nugget=1,scale=I)
start = list(mean = mean(z),sill=var(z),nugget=0.10,scale=140)
fixed = list(smooth = 0.5)
corrmodel="Matern"
pcl1=GeoFit2(coordx=loc_aux,corrmodel=corrmodel,data=z,
             likelihood="Full",type="Standard",model="Gaussian",
             optimizer="nlminb",lower=lower,upper=upper,sensitivity=TRUE,
             neighb=3, start=start,fixed=fixed)
pcl1
print(pcl1)
#GeoCV(pcl1,K=10,estimation=TRUE)

#############   Experiment
fit_variogram_quad <- fit.variogram(variogram_quad, model = vgm(model = "Mat"))
fit_variogram_knn <- fit.variogram(variogram_knn, model = vgm(model = "Mat"))
gstat_quad_vario <- gstat(formula = z ~ lon + I(lon^2) + lat + I(lat^2) + lon*lat,
                          locations = ~lon+lat,
                          data = anomalies_df,
                          model = fit_variogram_quad)
gstat_knn_vario <- gstat(formula = z ~anomalies_aux ,
                         locations = ~lon+lat,
                         data = anomalies_df,
                         model = fit_variogram_knn)


nugget <- gstat_quad_vario$psill[1] 
sill <- nugget + gstat_quad_vario$psill[2]
range <- gstat_quad_vario$range[2]
plot(variogram_quad,fit_variogram_quad,xlab = "Distancia", ylab = "Semivarianza", main = "Modelo 1")
plot(variogram_knn,fit_variogram_knn,xlab = "Distancia", ylab = "Semivarianza", main = "Modelo 2")


anomalies_sf <- st_as_sf(anomalies_df, coords = c("lon", "lat"), remove = FALSE, agr = "constant")
buffer <- anomalies_sf %>% st_geometry() %>% st_buffer(1)
grid <- buffer %>% st_as_stars(nx = 100, ny = 100)
coord <- st_coordinates(grid)
grid$lon <- coord$x
grid$lat <- coord$y
grid <- grid %>% st_crop(buffer)
pred <- krige(formula = z ~ lon + I(lon^2) + lat + I(lat^2) + I(lon*lat), locations = anomalies_sf, model = fit_variogram_quad, newdata=grid)


grid$var1.pred <- pred$var1.pred
grid$var1.var <- pred$var1.var
plot(grid["var1.pred"], breaks = "equal", col = sf.colors(64), key.pos = 4, main = "Predicciones kriging")


plot(grid["var1.var"], breaks = "equal", col = sf.colors(64), key.pos = 4, main = "Varianzas kriging")
### map data

x <-anomalies_sf$lon
y <-anomalies_sf$lat
p1 <- ggplot() + geom_stars(data = grid, aes(fill = var1.pred, x = x, y = y)) + scale_fill_viridis_c() + geom_sf(data = anomalies_sf) + coord_sf(lims_method = "geometry_bbox")
p2 <- ggplot() + geom_stars(data = grid, aes(fill = var1.var, x = x, y = y)) + scale_fill_viridis_c() + geom_sf(data = anomalies_sf) + coord_sf(lims_method = "geometry_bbox") 
grid.arrange(p1, p2, ncol = 2)

#############################################
anomalies_sf1 <- st_as_sf(anomalies_df, coords = c("lon", "lat"), remove = FALSE, agr = "constant")
buffer1 <- anomalies_sf1 %>% st_geometry() %>% st_buffer(40)
grid2 <- buffer1 %>% st_as_stars(nx = 50, ny = 50)
coord <- st_coordinates(grid2)
grid2$lon <- coord$x
grid2$lat <- coord$y
grid2 <- grid2 %>% st_crop(buffer1)
pred1 <- krige(formula = z ~anomalies_aux, locations = anomalies_sf1, model = fit_variogram_knn, newdata=grid2)


grid2$var1.pred <- pred1$var1.pred
grid2$var1.var <- pred1$var1.var
plot(grid2["var1.pred1"], breaks = "equal", col = sf.colors(64), key.pos = 4, main = "Predicciones kriging")


plot(grid2["var1.var"], breaks = "equal", col = sf.colors(64), key.pos = 4, main = "Varianzas kriging")
##################################################################################
# Prueba 3, tutorial 2:
#############################################################



image(pred, main="Predicciones con Kriging universal")

# Add a grid
grid()

# Add x-axis
axis(side = 1, at = seq(-140, -60, length.out = 5),
     labels = seq(-140, -60, length.out = 5),tick=10)
mtext("Lon", side = 1, line = 3)



# Add y-axis
axis(side = 2, at = seq(20, 50, length.out = 5),
     labels = round(seq(20, 50, length.out = 5), 0))
mtext("Lat", side = 2, line = 3)


map("usa", add = TRUE)
