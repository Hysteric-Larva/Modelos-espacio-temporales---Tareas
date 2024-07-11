#Libraries and data loading.
library(spatstat)
data(amacrine)
###################################################
###   First step, exploratory data analysis     ###
###################################################
#First, I want to inspect how it looks the intensity function by category of the data.
#   1° category visualization:
data(finpines)
finpines <- unmark(finpines)
plot(finpines)



plot(density(finpines))
###################################################


pcf.fun <- pcf(finpines)
Kfun.finpines <- Kest(finpines)
Gfun.finpines <- Gest(finpines)
Ffun.finpines <- Fest(finpines)
Jfun.finpines <- Jest(finpines)

par(mfrow = c(3,2))

plot(pcf.fun$r, pcf.fun$iso, type = "l")
lines(pcf.fun$r, pcf.fun$theo, type = "l", col = "red")

plot( Kfun.finpines, iso ~ r, type = "l")
lines( Kfun.finpines$r, Kfun.finpines$theo , col = "red")

plot( Gfun.finpines$r , Gfun.finpines$rs, type = "l", ylim = c(0,1))
lines( Gfun.finpines$r, Gfun.finpines$theo , col = "red")

plot( Ffun.finpines , rs ~ r, ylim = c(0,1))
lines( Ffun.finpines$r, Ffun.finpines$theo , col = "red")

plot( Jfun.finpines , rs ~ r )
lines( Jfun.finpines$r, Jfun.finpines
       $theo , col = "red")

num_points <- npoints(finpines)
###########################TEST
qTest <- quadrat.test(finpines, nx = 5, ny = 5)
print(qTest)


prueba <- clarkevans.test(finpines)
print(prueba)

################Thomas g#################
model_thomas <- thomas.estpcf(finpines, c(kappa=10, scale=0.1))
summary(model_thomas)
plot(model_thomas, main = "Modelo de Thomas ajustado",
     legendpos = "topright", nlabels = 3)


# Extraer los componentes del modelo
r_thomas <- model_thomas$fit$r
theo_thomas <- model_thomas$fit$theo
trans_thomas <- model_thomas$fit$trans
iso_thomas <- model_thomas$fit$iso
fit_thomas <- model_thomas$fit$fit

windows(width = 8, height = 6) 
plot(r_thomas,trans_thomas, main = "PCF Ajustado un proceso de Thomas", col = "green",xlim=c(0,2.5),xlab = "Metros",ylab = "PCF g(r)")
lines(r_thomas,theo_thomas, main = "PCF para el modelo teórico", col = "red",xlim=c(0,2.5))
lines(r_thomas,iso_thomas, main = "PCF para el modelo isotrópico", col = "purple",xlim=c(0,2.5))
lines(r_thomas,fit_thomas, main = "PCF para el modelo ajustado", col = "black",xlim=c(0,2.5))

# Agregar leyenda
legend("topright", legend = c("g Ajustado", "g Teorico", "g Transformado", "g Isotrópico"), 
       col = c("black", "red", "green", "purple"), lwd = 2, cex = 0.8)
#################################### Thomas K       ############################################
########################  COx K    #############################################################
model_thomas_k <- thomas.estK(finpines)

# Extraer los componentes del modelo
r__thomas_k <- model_thomas_k$fit$r
theo_thomas_k <- model_thomas_k$fit$theo
trans_thomas_k <- model_thomas_k$fit$trans
iso_thomas_k <- model_thomas_k$fit$iso
fit_thomas_k <- model_thomas_k$fit$fit

windows(width = 8, height = 6) 
plot(r__thomas_k,trans_thomas_k, main = "K ajustado un proceso de Thomas", col = "green",xlim=c(0,2.5),xlab = "Metros",ylab = "K(r)")
lines(r__thomas_k,theo_thomas_k, main = "K para el modelo teorico", col = "red",xlim=c(0,2.5))
lines(r__thomas_k,iso_thomas_k, main = "K para el modelo isotrópico", col = "purple",xlim=c(0,2.5))
lines(r__thomas_k,fit_thomas_k, main = "K para el modelo ajustado", col = "black",xlim=c(0,2.5))

# Agregar leyenda
legend("topright", legend = c("K Ajustado", "K Teorico", "K Transformado", "K Isotrópico"), 
       col = c("black", "red", "green", "purple"), lwd = 2, cex = 0.8)













################COX g #############################
model_cox <- lgcp.estpcf(finpines)



# Extraer los componentes del modelo
r_cox <- model_cox$fit$r
theo_cox <- model_cox$fit$theo
trans_cox <- model_cox$fit$trans
iso_cox <- model_cox$fit$iso
fit_cox <- model_cox$fit$fit


windows(width = 8, height = 6) 
plot(r_cox,trans_cox, main = "PCF ajustado un proceso de Cox log-gaussiano", col = "green",xlim=c(0,2.5),xlab = "Metros",ylab = "PCF g(r)")
lines(r_cox,theo_cox, main = "PCF para el modelo teórico", col = "red",xlim=c(0,2.5))
lines(r_cox,iso_cox, main = "PCF para el modelo isotrópico", col = "purple",xlim=c(0,2.5))
lines(r_cox,fit_cox, main = "PCF para el modelo ajustado", col = "black",xlim=c(0,2.5))

# Agregar leyenda
legend("topright", legend = c("g Ajustado", "g Teorico", "g Transformado", "g Isotrópico"), 
       col = c("black", "red", "green", "purple"), lwd = 2, cex = 0.8)

########################  COx K    #############################################################
model_cox_k <- lgcp.estK(finpines)

# Extraer los componentes del modelo
r_cox_k <- model_cox_k$fit$r
theo_cox_k <- model_cox_k$fit$theo
trans_cox_k <- model_cox_k$fit$trans
iso_cox_k <- model_cox_k$fit$iso
fit_cox_k <- model_cox_k$fit$fit

windows(width = 8, height = 6) 
plot(r_cox_k,trans_cox_k, main = "K ajustado un proceso de Cox log-gaussiano", col = "green",xlim=c(0,2.5),xlab = "Metros",ylab = "K(r)")
lines(r_cox_k,theo_cox_k, main = "K para el modelo teorico", col = "red",xlim=c(0,2.5))
lines(r_cox_k,iso_cox_k, main = "K para el modelo isotrópico", col = "purple",xlim=c(0,2.5))
lines(r_cox_k,fit_cox_k, main = "K para el modelo ajustado", col = "black",xlim=c(0,2.5))

# Agregar leyenda
legend("topright", legend = c("K Ajustado", "K Teorico", "K Transformado", "K Isotrópico"), 
       col = c("black", "red", "green", "purple"), lwd = 2, cex = 0.8)
