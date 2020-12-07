library(readstata13)
library("rEDM")
library("nonlinearTseries")
library("astsa")


library(readxl)
from_psych_weather <- read_excel("C:/Users/efa036/OneDrive - UiT Office 365/Analyse av psych datasett/from_psych_weather.xlsx")
View(from_psych_weather)
test <- subset(from_psych_weather, ID=="32")

testframe <- data.frame(time = test$measurenr,
                        painintensity=test$painintensity,
                        painunpleasantness=test$painunpleasantness,
                        pom=test$barometric_pres,
                        rh=test$humidity,
                        tam=test$Temperature)


n <- NROW(testframe) 
lib <- c(1, floor(2/3 * n)) # indices for the first 2/3 of the time series
pred <- c(floor(2/3 * n) + 1, n) # indices for the final 1/3 of the time series
acf(testframe$painintensity)
bestpredictions <- EmbedDimension(dataFrame=testframe, tau = -1, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 10)
bestpredictions <- bestpredictions[2]

for (t in 2:3) {
  best_out <-   EmbedDimension(dataFrame=testframe, tau = -t, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 10)
  bestpredictions <- cbind(bestpredictions, best_out[2])
}
max <- which(bestpredictions==max(bestpredictions), arr.ind = TRUE)
max
bestpredictions[max]

for (t in 3:15) {
  best_out <-   EmbedDimension(dataFrame=testframe, tau = -t, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 8)
  bestpredictions <- cbind(bestpredictions, best_out[2])
}
max <- which(bestpredictions==max(bestpredictions), arr.ind = TRUE)
max
bestpredictions[max]

tau = -4 
#### Embeddimension painintensity #### 
E.optpainintensity = EmbedDimension(dataFrame=testframe, tau = tau, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 10)
E.optpainintensity 
E=2
embed <- Embed(dataFrame = testframe, tau = tau, E = E, columns = "painintensity" )
plot(embed$`painintensity(t-0)`, embed$`painintensity(t-1)`, type = "l")

#### simplex projections painintensity #### 
simplextrykk <- Simplex(dataFrame=testframe, tau = tau, 
                        lib=lib, pred=pred, 
                        columns = "painintensity", target = "painintensity", 
                        E = E ) 
plot(testframe$time, testframe$painintensity, type = "l", lwd = 2, xlab = "time", ylab = "painintensity") 
lines(simplextrykk$time, simplextrykk$Predictions, col = "red", lwd = 2) 
legend( 'topleft', legend = c( "Observed", "Predicted (times + 1)" ), fill = c( 'black', 'red' ))

plot(simplextrykk$Predictions, simplextrykk$Observations)
#### predict interval painintensity #### 
predictint <- PredictInterval(dataFrame=testframe, tau = tau, 
                              lib=lib, pred=pred, 
                              columns = "painintensity", target = "painintensity", 
                              E = E )
print(predictint)
#### non-linearity painintensity #### 
prednonlinear <- PredictNonlinear(dataFrame=testframe, tau = tau, knn=50,
                                  lib="1 100", pred="1 100", 
                                  columns = "painintensity", target = "painintensity", 
                                  E = E )


##### CCM painintensity #### 
vars = colnames(testframe[2:6]) 
var_pairs = combn(vars, 2) # Combinations of vars, 2 at a time 
libSize = paste(NROW(testframe)-abs(tau*E), NROW(testframe)-abs(tau*E), 10, collapse = "" ) 
ccm_matrix = array(NA, dim= c(length(vars), length(vars)), dimnames = list(vars, vars)) 

for (i in 1:ncol(var_pairs)) { 
  ccm_out = CCM(dataFrame = testframe, columns = var_pairs[1, i], target = var_pairs[2,i], libSizes = libSize, Tp = 0, E = E, tau=tau, sample = 100) 
  outVars = names(ccm_out) 
  var_out = unlist(strsplit(outVars[2], ":")) 
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 2] 
  var_out = unlist(strsplit(outVars[3], ":")) 
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 3] }

corr_matrix <- array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars, vars)) 
for (ccm_from in vars) { 
  for (ccm_to in vars[vars != ccm_from]) {
    ccf_out <- ccf(testframe[, ccm_from], testframe[, ccm_to], type = "correlation", lag.max = 20, plot = FALSE)$acf 
    
    corr_matrix[ccm_from, ccm_to] <- max(abs(ccf_out)) } } 

print(ccm_matrix)
print(corr_matrix)

Tp <- -1 
painintensityxtam <-CCM(dataFrame = testframe,  E = E, Tp=Tp,  column = "tam", tau=tau, target = "painintensity",  libSizes="10 100 10", sample=100 ) 
plot(painintensityxtam $LibSize, painintensityxtam $`painintensity:tam`, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painintensityxtam $LibSize, painintensityxtam $`tam:painintensity`, lwd=2, col = "blue") 
legend('topleft', legend = c( "painintensity:tam", "tam:painintensity" ),  fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painintensity", "tam"], col = "black", lty = 2)

painintensityxpom <-CCM(dataFrame = testframe, E = E , Tp=Tp,  column = "pom", tau=tau, target = "painintensity", libSizes="10 100 10", sample=100 ) 
plot(painintensityxpom $LibSize, painintensityxpom $`painintensity:pom`, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painintensityxpom $LibSize, painintensityxpom $`pom:painintensity`, lwd=2, col = "blue") 
legend('topleft', legend = c( "painintensity:pom", "pom:painintensity" ), fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painintensity", "pom"], col = "black", lty = 2)


painintensityxrh <- CCM(dataFrame = testframe, E = E , column = "rh", tau = tau,  Tp=Tp, 
                        target = "painintensity", libSizes="10 100 10", sample=100 ) 
plot(painintensityxrh $LibSize, painintensityxrh $painintensity, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painintensityxrh $LibSize, painintensityxrh $rh, lwd=2, col = "blue") 
legend('topleft', legend = c( "painintensity:rh", "rh:painintensity" ), fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painintensity", "rh"], col = "black", lty = 2)

for (i in 1/)
  
  smap <- SMap(dataFrame = testframe, E = E ,  tau = tau, knn = 50,
               lib=lib, pred=pred, 
               columns = "painintensity tam pom rh", target = "painintensity",  )
predictions = smap$predictions 
coefficients = smap$coefficients 
time = smap$predictions$time
plot(time, predictions$Observations, type = "l", col = "blue", ylab = "V1", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
lines(time, predictions$Predictions, lwd = 2, col = "red") 
legend("topright", legend = c("observed", "predicted"), fill = c("blue", "red"), bty = "n", cex = 1.3)
plot(time, coefficients$C0, type = "l", col = "brown", ylab = "dTAM/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
plot(time, coefficients$C1, type = "l", col = "darkgreen", ylab = "dPOM/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
plot(time, coefficients$C2, type = "l", col = "blue", ylab = "drh/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3)

##### Unpleasentness ####################
test <- subset(from_psych_weather, ID=="32")

testframe <- data.frame(time = test$measurenr,
                        painintensity=test$painintensity,
                        painunpleasantness=test$painunpleasantness	,
                        pom=test$barometric_pres,
                        rh=test$humidity,
                        tam=test$Temperature)
acf2(testframe$painunpleasantness	, max.lag = 30)
bestpredictions <- EmbedDimension(dataFrame=testframe, tau = -1, lib=lib, pred=pred, columns = "painunpleasantness", target = "painunpleasantness", maxE = 10)
bestpredictions <- bestpredictions[2]

for (t in 2:3) {
  best_out <-   EmbedDimension(dataFrame=testframe, tau = -t, lib=lib, pred=pred, columns = "painunpleasantness", target = "painunpleasantness", maxE = 10)
  bestpredictions <- cbind(bestpredictions, best_out[2])
}
max <- which(bestpredictions==max(bestpredictions), arr.ind = TRUE)
max
bestpredictions[max]

for (t in 3:15) {
  best_out <-   EmbedDimension(dataFrame=testframe, tau = -t, lib=lib, pred=pred, columns = "painunpleasantness", target = "painunpleasantness", maxE = 8)
  bestpredictions <- cbind(bestpredictions, best_out[2])
}
max <- which(bestpredictions==max(bestpredictions), arr.ind = TRUE)
max
bestpredictions[max]
bestpredictions <- EmbedDimension(dataFrame=testframe, tau = -4, lib=lib, pred=pred, columns = "painunpleasantness", target = "painunpleasantness", maxE = 10)

tau = -1 
#### Embeddimension painunpleasantness #### 
E.optpainunpleasantness = EmbedDimension(dataFrame=testframe, tau = tau, lib=lib, pred=pred, columns = "painunpleasantness", target = "painunpleasantness", maxE = 10)
E.optpainunpleasantness 
E=5
embed <- Embed(dataFrame = testframe, tau = tau, E = E, columns = "painunpleasantness" )
plot(embed$`painunpleasantness(t-0)`, embed$`painunpleasantness(t-1)`, type = "l")

#### simplex projections painunpleasantness #### 
simplextrykk <- Simplex(dataFrame=testframe, tau = tau, 
                        lib=lib, pred=pred, 
                        columns = "painunpleasantness", target = "painunpleasantness", 
                        E = E ) 
plot(testframe$time, testframe$painunpleasantness, type = "l", lwd = 2, xlab = "time", ylab = "painunpleasantness") 
lines(simplextrykk$time, simplextrykk$Predictions, col = "red", lwd = 2) 
legend( 'topleft', legend = c( "Observed", "Predicted (times + 1)" ), fill = c( 'black', 'red' ))

plot(simplextrykk$Predictions, simplextrykk$Observations)
#### predict interval painunpleasantness #### 
predictint <- PredictInterval(dataFrame=testframe, tau = tau, 
                              lib=lib, pred=pred, 
                              columns = "painunpleasantness", target = "painunpleasantness", 
                              E = E )
print(predictint)
#### non-linearity painunpleasantness #### 
prednonlinear <- PredictNonlinear(dataFrame=testframe, tau = tau, knn=50,
                                  lib="1 100", pred="1 100", 
                                  columns = "painunpleasantness", target = "painunpleasantness", 
                                  E = E )


##### CCM painunpleasantness #### 
vars = colnames(testframe[2:6]) 
var_pairs = combn(vars, 2) # Combinations of vars, 2 at a time 
libSize = paste(NROW(testframe)-abs(tau*E), NROW(testframe)-abs(tau*E), 10, collapse = "" ) 
ccm_matrix = array(NA, dim= c(length(vars), length(vars)), dimnames = list(vars, vars)) 

for (i in 1:ncol(var_pairs)) { 
  ccm_out = CCM(dataFrame = testframe, columns = var_pairs[1, i], target = var_pairs[2,i], libSizes = libSize, Tp = 0, E = E, tau=tau, sample = 100) 
  outVars = names(ccm_out) 
  var_out = unlist(strsplit(outVars[2], ":")) 
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 2] 
  var_out = unlist(strsplit(outVars[3], ":")) 
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 3] }

corr_matrix <- array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars, vars)) 
for (ccm_from in vars) { 
  for (ccm_to in vars[vars != ccm_from]) {
    ccf_out <- ccf(testframe[, ccm_from], testframe[, ccm_to], type = "correlation", lag.max = 20, plot = FALSE)$acf 
    
    corr_matrix[ccm_from, ccm_to] <- max(abs(ccf_out)) } } 

print(ccm_matrix)
print(corr_matrix)

painunpleasantnessxtam <-CCM(dataFrame = testframe,  E = E, column = "tam", tau=tau, target = "painunpleasantness",  libSizes="10 100 10", sample=100 ) 
plot(painunpleasantnessxtam $LibSize, painunpleasantnessxtam $`painunpleasantness:tam`, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painunpleasantnessxtam $LibSize, painunpleasantnessxtam $`tam:painunpleasantness`, lwd=2, col = "blue") 
legend('topleft', legend = c( "painunpleasantness:tam", "tam:painunpleasantness" ),  fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painunpleasantness", "tam"], col = "black", lty = 2)

painunpleasantnessxpom <-CCM(dataFrame = testframe, E = E ,Tp=1 ,column = "pom", tau=tau, target = "painunpleasantness", libSizes="10 100 10", sample=100 ) 
plot(painunpleasantnessxpom $LibSize, painunpleasantnessxpom $`painunpleasantness:pom`, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painunpleasantnessxpom $LibSize, painunpleasantnessxpom $`pom:painunpleasantness`, lwd=2, col = "blue") 
legend('topleft', legend = c( "painunpleasantness:pom", "pom:painunpleasantness" ), fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painunpleasantness", "pom"], col = "black", lty = 2)


painunpleasantnessxrh <- CCM(dataFrame = testframe, E = E , column = "rh", tau = tau, 
                             target = "painunpleasantness", libSizes="10 100 10", sample=100 ) 
plot(painunpleasantnessxrh $LibSize, painunpleasantnessxrh $painunpleasantness, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painunpleasantnessxrh $LibSize, painunpleasantnessxrh $rh, lwd=2, col = "blue") 
legend('topleft', legend = c( "painunpleasantness:rh", "rh:painunpleasantness" ), fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painunpleasantness", "rh"], col = "black", lty = 2)


smap <- SMap(dataFrame = testframe, E = E ,  tau = tau, knn = 50,
             lib=lib, pred=pred, 
             columns = "painunpleasantness tam pom rh", target = "painunpleasantness",  )
predictions = smap$predictions 
coefficients = smap$coefficients 
time = smap$predictions$time
plot(time, predictions$Observations, type = "l", col = "blue", ylab = "V1", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
lines(time, predictions$Predictions, lwd = 2, col = "red") 
legend("topright", legend = c("observed", "predicted"), fill = c("blue", "red"), bty = "n", cex = 1.3)
plot(time, coefficients$C0, type = "l", col = "brown", ylab = "dTAM/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
plot(time, coefficients$C1, type = "l", col = "darkgreen", ylab = "dPOM/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
plot(time, coefficients$C2, type = "l", col = "blue", ylab = "drh/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3)


################## ID 42 #######################
test <- subset(from_psych_weather, ID=="42")

testframe <- data.frame(time = test$measurenr,
                        painintensity=test$painintensity,
                        painunpleasantness=test$painunpleasantness,
                        pom=test$barometric_pres,
                        rh=test$humidity,
                        tam=test$Temperature)


n <- NROW(testframe) 
lib <- c(1, floor(2/3 * n)) # indices for the first 2/3 of the time series
pred <- c(floor(2/3 * n) + 1, n) # indices for the final 1/3 of the time series
acf2(testframe$painintensity)
EmbedDimension(dataFrame=testframe, tau = -1, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 10)
EmbedDimension(dataFrame=testframe, tau = -2, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 10)
EmbedDimension(dataFrame=testframe, tau = -3, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 10)
EmbedDimension(dataFrame=testframe, tau = -4, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 8)


tau = -4 
#### Embeddimension painintensity #### 
E.optpainintensity = EmbedDimension(dataFrame=testframe, tau = tau, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 10)
E.optpainintensity 
E=2
embed <- Embed(dataFrame = testframe, tau = tau, E = E, columns = "painintensity" )
plot(embed$`painintensity(t-0)`, embed$`painintensity(t-1)`, type = "l")

#### simplex projections painintensity #### 
simplextrykk <- Simplex(dataFrame=testframe, tau = tau, 
                        lib=lib, pred=pred, 
                        columns = "painintensity", target = "painintensity", 
                        E = E ) 
plot(testframe$time, testframe$painintensity, type = "l", lwd = 2, xlab = "time", ylab = "painintensity") 
lines(simplextrykk$time, simplextrykk$Predictions, col = "red", lwd = 2) 
legend( 'topleft', legend = c( "Observed", "Predicted (times + 1)" ), fill = c( 'black', 'red' ))

plot(simplextrykk$Predictions, simplextrykk$Observations)
#### predict interval painintensity #### 
predictint <- PredictInterval(dataFrame=testframe, tau = tau, 
                              lib=lib, pred=pred, 
                              columns = "painintensity", target = "painintensity", 
                              E = E )
print(predictint)
#### non-linearity painintensity #### 
prednonlinear <- PredictNonlinear(dataFrame=testframe, tau = tau, knn=50,
                                  lib="1 70", pred="1 87", 
                                  columns = "painintensity", target = "painintensity", 
                                  E = E )


##### CCM painintensity #### 
vars = colnames(testframe[2:6]) 
var_pairs = combn(vars, 2) # Combinations of vars, 2 at a time 
libSize = paste(NROW(testframe)-abs(tau*E), NROW(testframe)-abs(tau*E), 10, collapse = "" ) 
ccm_matrix = array(NA, dim= c(length(vars), length(vars)), dimnames = list(vars, vars)) 

for (i in 1:ncol(var_pairs)) { 
  ccm_out = CCM(dataFrame = testframe, columns = var_pairs[1, i], target = var_pairs[2,i], libSizes = libSize, Tp = 0, E = E, tau=tau, sample = 100) 
  outVars = names(ccm_out) 
  var_out = unlist(strsplit(outVars[2], ":")) 
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 2] 
  var_out = unlist(strsplit(outVars[3], ":")) 
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 3] }

corr_matrix <- array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars, vars)) 
for (ccm_from in vars) { 
  for (ccm_to in vars[vars != ccm_from]) {
    ccf_out <- ccf(testframe[, ccm_from], testframe[, ccm_to], type = "correlation", lag.max = 20, plot = FALSE)$acf 
    
    corr_matrix[ccm_from, ccm_to] <- max(abs(ccf_out)) } } 

print(ccm_matrix)
print(corr_matrix)

painintensityxtam <-CCM(dataFrame = testframe,  E = E, column = "tam", tau=tau, target = "painintensity",  libSizes="10 72 10", sample=72 ) 
plot(painintensityxtam $LibSize, painintensityxtam $`painintensity:tam`, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painintensityxtam $LibSize, painintensityxtam $`tam:painintensity`, lwd=2, col = "blue") 
legend('topleft', legend = c( "painintensity:tam", "tam:painintensity" ),  fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painintensity", "tam"], col = "black", lty = 2)

painintensityxpom <-CCM(dataFrame = testframe, E = E , column = "pom", tau=tau, target = "painintensity", libSizes="10 72 10", sample=72 ) 
plot(painintensityxpom $LibSize, painintensityxpom $`painintensity:pom`, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painintensityxpom $LibSize, painintensityxpom $`pom:painintensity`, lwd=2, col = "blue") 
legend('topleft', legend = c( "painintensity:pom", "pom:painintensity" ), fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painintensity", "pom"], col = "black", lty = 2)


painintensityxrh <- CCM(dataFrame = testframe, E = E , column = "rh", tau = tau, 
                        target = "painintensity", libSizes="10 72 10", sample=72 ) 
plot(painintensityxrh $LibSize, painintensityxrh $painintensity, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painintensityxrh $LibSize, painintensityxrh $rh, lwd=2, col = "blue") 
legend('topleft', legend = c( "painintensity:rh", "rh:painintensity" ), fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painintensity", "rh"], col = "black", lty = 2)


smap <- SMap(dataFrame = testframe, E = E ,  tau = tau, knn = 50,
             lib=lib, pred=pred, 
             columns = "painintensity tam pom rh", target = "painintensity",  )
predictions = smap$predictions 
coefficients = smap$coefficients 
time = smap$predictions$time
plot(time, predictions$Observations, type = "l", col = "blue", ylab = "V1", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
lines(time, predictions$Predictions, lwd = 2, col = "red") 
legend("topright", legend = c("observed", "predicted"), fill = c("blue", "red"), bty = "n", cex = 1.3)
plot(time, coefficients$C0, type = "l", col = "brown", ylab = "dTAM/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
plot(time, coefficients$C1, type = "l", col = "darkgreen", ylab = "dPOM/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
plot(time, coefficients$C2, type = "l", col = "blue", ylab = "drh/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3)




################## ID 33 #######################
test <- subset(from_psych_weather, ID=="33")

testframe <- data.frame(time = test$measurenr,
                        painintensity=test$painintensity,
                        painunpleasantness=test$painunpleasantness,
                        pom=test$barometric_pres,
                        rh=test$humidity,
                        tam=test$Temperature)


n <- NROW(testframe) 
lib <- c(1, floor(2/3 * n)) # indices for the first 2/3 of the time series
pred <- c(floor(2/3 * n) + 1, n) # indices for the final 1/3 of the time series
acf2(testframe$painintensity)
EmbedDimension(dataFrame=testframe, tau = -1, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 10)
EmbedDimension(dataFrame=testframe, tau = -2, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 6)


tau = -2 
#### Embeddimension painintensity #### 
E.optpainintensity = EmbedDimension(dataFrame=testframe, tau = tau, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 10)
E.optpainintensity 
E=3
embed <- Embed(dataFrame = testframe, tau = tau, E = E, columns = "painintensity" )
plot(embed$`painintensity(t-0)`, embed$`painintensity(t-1)`, type = "l")

#### simplex projections painintensity #### 
simplextrykk <- Simplex(dataFrame=testframe, tau = tau, 
                        lib=lib, pred=pred, 
                        columns = "painintensity", target = "painintensity", 
                        E = E ) 
plot(testframe$time, testframe$painintensity, type = "l", lwd = 2, xlab = "time", ylab = "painintensity") 
lines(simplextrykk$time, simplextrykk$Predictions, col = "red", lwd = 2) 
legend( 'topleft', legend = c( "Observed", "Predicted (times + 1)" ), fill = c( 'black', 'red' ))

plot(simplextrykk$Predictions, simplextrykk$Observations)
#### predict interval painintensity #### 
predictint <- PredictInterval(dataFrame=testframe, tau = tau, 
                              lib=lib, pred=pred, 
                              columns = "painintensity", target = "painintensity", 
                              E = E )
print(predictint)
#### non-linearity painintensity #### 
prednonlinear <- PredictNonlinear(dataFrame=testframe, tau = tau, knn=16,
                                  lib="1 26", pred="1 26", 
                                  columns = "painintensity", target = "painintensity", 
                                  E = E )


##### CCM painintensity #### 
vars = colnames(testframe[2:6]) 
var_pairs = combn(vars, 2) # Combinations of vars, 2 at a time 
libSize = paste(NROW(testframe)-abs(tau*E), NROW(testframe)-abs(tau*E), 10, collapse = "" ) 
ccm_matrix = array(NA, dim= c(length(vars), length(vars)), dimnames = list(vars, vars)) 

for (i in 1:ncol(var_pairs)) { 
  ccm_out = CCM(dataFrame = testframe, columns = var_pairs[1, i], target = var_pairs[2,i], libSizes = libSize, Tp = 0, E = E, tau=tau, sample = 100) 
  outVars = names(ccm_out) 
  var_out = unlist(strsplit(outVars[2], ":")) 
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 2] 
  var_out = unlist(strsplit(outVars[3], ":")) 
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 3] }

corr_matrix <- array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars, vars)) 
for (ccm_from in vars) { 
  for (ccm_to in vars[vars != ccm_from]) {
    ccf_out <- ccf(testframe[, ccm_from], testframe[, ccm_to], type = "correlation", lag.max = 20, plot = FALSE)$acf 
    
    corr_matrix[ccm_from, ccm_to] <- max(abs(ccf_out)) } } 

print(ccm_matrix)
print(corr_matrix)

painintensityxtam <-CCM(dataFrame = testframe,  E = E, column = "tam", tau=tau, target = "painintensity",  libSizes="10 26 10", sample=26 ) 
plot(painintensityxtam $LibSize, painintensityxtam $`painintensity:tam`, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painintensityxtam $LibSize, painintensityxtam $`tam:painintensity`, lwd=2, col = "blue") 
legend('topleft', legend = c( "painintensity:tam", "tam:painintensity" ),  fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painintensity", "tam"], col = "black", lty = 2)

painintensityxpom <-CCM(dataFrame = testframe, E = E , column = "pom", tau=tau, target = "painintensity", libSizes="10 26 10", sample=26 ) 
plot(painintensityxpom $LibSize, painintensityxpom $`painintensity:pom`, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painintensityxpom $LibSize, painintensityxpom $`pom:painintensity`, lwd=2, col = "blue") 
legend('topleft', legend = c( "painintensity:pom", "pom:painintensity" ), fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painintensity", "pom"], col = "black", lty = 2)


painintensityxrh <- CCM(dataFrame = testframe, E = E ,Tp=1, column = "rh", tau = tau, 
                        target = "painintensity", libSizes="10 26 10", sample=26 ) 
plot(painintensityxrh $LibSize, painintensityxrh $painintensity, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painintensityxrh $LibSize, painintensityxrh $rh, lwd=2, col = "blue") 
legend('topleft', legend = c( "painintensity:rh", "rh:painintensity" ), fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painintensity", "rh"], col = "black", lty = 2)


smap <- SMap(dataFrame = testframe, E = E ,  tau = tau, knn = 50,
             lib=lib, pred=pred, 
             columns = "painintensity tam pom rh", target = "painintensity",  )
predictions = smap$predictions 
coefficients = smap$coefficients 
time = smap$predictions$time
plot(time, predictions$Observations, type = "l", col = "blue", ylab = "V1", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
lines(time, predictions$Predictions, lwd = 2, col = "red") 
legend("topright", legend = c("observed", "predicted"), fill = c("blue", "red"), bty = "n", cex = 1.3)
plot(time, coefficients$C0, type = "l", col = "brown", ylab = "dTAM/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
plot(time, coefficients$C1, type = "l", col = "darkgreen", ylab = "dPOM/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
plot(time, coefficients$C2, type = "l", col = "blue", ylab = "drh/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3)




################## ID 28 #######################
test <- subset(from_psych_weather, ID=="28")

testframe <- data.frame(time = test$measurenr,
                        painintensity=test$painintensity,
                        painunpleasantness=test$painunpleasantness,
                        pom=test$barometric_pres,
                        rh=test$humidity,
                        tam=test$Temperature)


n <- NROW(testframe) 
lib <- c(1, floor(2/3 * n)) # indices for the first 2/3 of the time series
pred <- c(floor(2/3 * n) + 1, n) # indices for the final 1/3 of the time series
acf2(testframe$painintensity)
EmbedDimension(dataFrame=testframe, tau = -1, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 10)
EmbedDimension(dataFrame=testframe, tau = -2, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 6)


tau = -2 
#### Embeddimension painintensity #### 
E.optpainintensity = EmbedDimension(dataFrame=testframe, tau = tau, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 10)
E.optpainintensity 
E=3
embed <- Embed(dataFrame = testframe, tau = tau, E = E, columns = "painintensity" )
plot(embed$`painintensity(t-0)`, embed$`painintensity(t-1)`, type = "l")

#### simplex projections painintensity #### 
simplextrykk <- Simplex(dataFrame=testframe, tau = tau, 
                        lib=lib, pred=pred, 
                        columns = "painintensity", target = "painintensity", 
                        E = E ) 
plot(testframe$time, testframe$painintensity, type = "l", lwd = 2, xlab = "time", ylab = "painintensity") 
lines(simplextrykk$time, simplextrykk$Predictions, col = "red", lwd = 2) 
legend( 'topleft', legend = c( "Observed", "Predicted (times + 1)" ), fill = c( 'black', 'red' ))

plot(simplextrykk$Predictions, simplextrykk$Observations)
#### predict interval painintensity #### 
predictint <- PredictInterval(dataFrame=testframe, tau = tau, 
                              lib=lib, pred=pred, 
                              columns = "painintensity", target = "painintensity", 
                              E = E )
print(predictint)
#### non-linearity painintensity #### 
prednonlinear <- PredictNonlinear(dataFrame=testframe, tau = tau, knn=50,
                                  lib="1 65", pred="1 65", 
                                  columns = "painintensity", target = "painintensity", 
                                  E = E )


##### CCM painintensity #### 
vars = colnames(testframe[2:6]) 
var_pairs = combn(vars, 2) # Combinations of vars, 2 at a time 
libSize = paste(NROW(testframe)-abs(tau*E), NROW(testframe)-abs(tau*E), 10, collapse = "" ) 
ccm_matrix = array(NA, dim= c(length(vars), length(vars)), dimnames = list(vars, vars)) 

for (i in 1:ncol(var_pairs)) { 
  ccm_out = CCM(dataFrame = testframe, columns = var_pairs[1, i], target = var_pairs[2,i], libSizes = libSize, Tp = 0, E = E, tau=tau, sample = 100) 
  outVars = names(ccm_out) 
  var_out = unlist(strsplit(outVars[2], ":")) 
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 2] 
  var_out = unlist(strsplit(outVars[3], ":")) 
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 3] }

corr_matrix <- array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars, vars)) 
for (ccm_from in vars) { 
  for (ccm_to in vars[vars != ccm_from]) {
    ccf_out <- ccf(testframe[, ccm_from], testframe[, ccm_to], type = "correlation", lag.max = 20, plot = FALSE)$acf 
    
    corr_matrix[ccm_from, ccm_to] <- max(abs(ccf_out)) } } 

print(ccm_matrix)
print(corr_matrix)
Tp <- -1
painintensityxtam <-CCM(dataFrame = testframe,  E = E, Tp=Tp, column = "tam", tau=tau, target = "painintensity",  libSizes="10 70 10", sample=70 ) 
plot(painintensityxtam $LibSize, painintensityxtam $`painintensity:tam`, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painintensityxtam $LibSize, painintensityxtam $`tam:painintensity`, lwd=2, col = "blue") 
legend('topleft', legend = c( "painintensity:tam", "tam:painintensity" ),  fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painintensity", "tam"], col = "black", lty = 2)

painintensityxpom <-CCM(dataFrame = testframe, E = E, Tp=Tp, column = "pom", tau=tau, target = "painintensity", libSizes="10 70 10", sample=70 ) 
plot(painintensityxpom $LibSize, painintensityxpom $`painintensity:pom`, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painintensityxpom $LibSize, painintensityxpom $`pom:painintensity`, lwd=2, col = "blue") 
legend('topleft', legend = c( "painintensity:pom", "pom:painintensity" ), fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painintensity", "pom"], col = "black", lty = 2)


painintensityxrh <- CCM(dataFrame = testframe, E = E ,Tp=Tp, column = "rh", tau = tau, 
                        target = "painintensity", libSizes="10 70 10", sample=70 ) 
plot(painintensityxrh $LibSize, painintensityxrh $painintensity, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painintensityxrh $LibSize, painintensityxrh $rh, lwd=2, col = "blue") 
legend('topleft', legend = c( "painintensity:rh", "rh:painintensity" ), fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painintensity", "rh"], col = "black", lty = 2)


smap <- SMap(dataFrame = testframe, E = E ,  tau = tau, knn = 50,
             lib=lib, pred=pred, 
             columns = "painintensity tam pom rh", target = "painintensity",  )
predictions = smap$predictions 
coefficients = smap$coefficients 
time = smap$predictions$time
plot(time, predictions$Observations, type = "l", col = "blue", ylab = "V1", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
lines(time, predictions$Predictions, lwd = 2, col = "red") 
legend("topright", legend = c("observed", "predicted"), fill = c("blue", "red"), bty = "n", cex = 1.3)
plot(time, coefficients$C0, type = "l", col = "brown", ylab = "dTAM/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
plot(time, coefficients$C1, type = "l", col = "darkgreen", ylab = "dPOM/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
plot(time, coefficients$C2, type = "l", col = "blue", ylab = "drh/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3)





################## ID 21 #######################
test <- subset(from_psych_weather, ID=="21")

testframe <- data.frame(time = test$measurenr,
                        painintensity=test$painintensity,
                        painunpleasantness=test$painunpleasantness,
                        pom=test$barometric_pres,
                        rh=test$humidity,
                        tam=test$Temperature)


n <- NROW(testframe) 
lib <- c(1, floor(2/3 * n)) # indices for the first 2/3 of the time series
pred <- c(floor(2/3 * n) + 1, n) # indices for the final 1/3 of the time series
acf2(testframe$painintensity)
EmbedDimension(dataFrame=testframe, tau = -1, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 10)
EmbedDimension(dataFrame=testframe, tau = -2, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 6)


tau = -1 
#### Embeddimension painintensity #### 
E.optpainintensity = EmbedDimension(dataFrame=testframe, tau = tau, lib=lib, pred=pred, columns = "painintensity", target = "painintensity", maxE = 10)
E.optpainintensity 
E=3
embed <- Embed(dataFrame = testframe, tau = tau, E = E, columns = "painintensity" )
plot(embed$`painintensity(t-0)`, embed$`painintensity(t-1)`, type = "l")

#### simplex projections painintensity #### 
simplextrykk <- Simplex(dataFrame=testframe, tau = tau, 
                        lib=lib, pred=pred, 
                        columns = "painintensity", target = "painintensity", 
                        E = E ) 
plot(testframe$time, testframe$painintensity, type = "l", lwd = 2, xlab = "time", ylab = "painintensity") 
lines(simplextrykk$time, simplextrykk$Predictions, col = "red", lwd = 2) 
legend( 'topleft', legend = c( "Observed", "Predicted (times + 1)" ), fill = c( 'black', 'red' ))

plot(simplextrykk$Predictions, simplextrykk$Observations)
#### predict interval painintensity #### 
predictint <- PredictInterval(dataFrame=testframe, tau = tau, 
                              lib=lib, pred=pred, 
                              columns = "painintensity", target = "painintensity", 
                              E = E )
print(predictint)
#### non-linearity painintensity #### 
prednonlinear <- PredictNonlinear(dataFrame=testframe, tau = tau, knn=50,
                                  lib="1 65", pred="1 65", 
                                  columns = "painintensity", target = "painintensity", 
                                  E = E )


##### CCM painintensity #### 
vars = colnames(testframe[2:6]) 
var_pairs = combn(vars, 2) # Combinations of vars, 2 at a time 
libSize = paste(NROW(testframe)-abs(tau*E), NROW(testframe)-abs(tau*E), 10, collapse = "" ) 
ccm_matrix = array(NA, dim= c(length(vars), length(vars)), dimnames = list(vars, vars)) 

for (i in 1:ncol(var_pairs)) { 
  ccm_out = CCM(dataFrame = testframe, columns = var_pairs[1, i], target = var_pairs[2,i], libSizes = libSize, Tp = 0, E = E, tau=tau, sample = 100) 
  outVars = names(ccm_out) 
  var_out = unlist(strsplit(outVars[2], ":")) 
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 2] 
  var_out = unlist(strsplit(outVars[3], ":")) 
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 3] }

corr_matrix <- array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars, vars)) 
for (ccm_from in vars) { 
  for (ccm_to in vars[vars != ccm_from]) {
    ccf_out <- ccf(testframe[, ccm_from], testframe[, ccm_to], type = "correlation", lag.max = 20, plot = FALSE)$acf 
    
    corr_matrix[ccm_from, ccm_to] <- max(abs(ccf_out)) } } 

print(ccm_matrix)
print(corr_matrix)

painintensityxtam <-CCM(dataFrame = testframe,  E = E, column = "tam", tau=tau, target = "painintensity",  libSizes="10 26 10", sample=26 ) 
plot(painintensityxtam $LibSize, painintensityxtam $`painintensity:tam`, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painintensityxtam $LibSize, painintensityxtam $`tam:painintensity`, lwd=2, col = "blue") 
legend('topleft', legend = c( "painintensity:tam", "tam:painintensity" ),  fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painintensity", "tam"], col = "black", lty = 2)

painintensityxpom <-CCM(dataFrame = testframe, E = E , column = "pom", tau=tau, target = "painintensity", libSizes="10 26 10", sample=26 ) 
plot(painintensityxpom $LibSize, painintensityxpom $`painintensity:pom`, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painintensityxpom $LibSize, painintensityxpom $`pom:painintensity`, lwd=2, col = "blue") 
legend('topleft', legend = c( "painintensity:pom", "pom:painintensity" ), fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painintensity", "pom"], col = "black", lty = 2)


painintensityxrh <- CCM(dataFrame = testframe, E = E ,Tp=1, column = "rh", tau = tau, 
                        target = "painintensity", libSizes="10 26 10", sample=26 ) 
plot(painintensityxrh $LibSize, painintensityxrh $painintensity, type = "l", col = "red", ylim=c(0, 1) ) 
lines (painintensityxrh $LibSize, painintensityxrh $rh, lwd=2, col = "blue") 
legend('topleft', legend = c( "painintensity:rh", "rh:painintensity" ), fill = c( 'red', 'blue' )) 
abline(h = corr_matrix["painintensity", "rh"], col = "black", lty = 2)


smap <- SMap(dataFrame = testframe, E = E ,  tau = tau, knn = 50,
             lib=lib, pred=pred, 
             columns = "painintensity tam pom rh", target = "painintensity",  )
predictions = smap$predictions 
coefficients = smap$coefficients 
time = smap$predictions$time
plot(time, predictions$Observations, type = "l", col = "blue", ylab = "V1", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
lines(time, predictions$Predictions, lwd = 2, col = "red") 
legend("topright", legend = c("observed", "predicted"), fill = c("blue", "red"), bty = "n", cex = 1.3)
plot(time, coefficients$C0, type = "l", col = "brown", ylab = "dTAM/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
plot(time, coefficients$C1, type = "l", col = "darkgreen", ylab = "dPOM/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3) 
plot(time, coefficients$C2, type = "l", col = "blue", ylab = "drh/dPPT", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3)
