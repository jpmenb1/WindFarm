#practical 1

#-----LIBRARIES-----------------------------------------------------------------
library(stats)
library(plyr)
library(ggplot2)
library(GGally)
library(tidyverse)
library(car)
library(MuMIn)


pract1mt5764<-  read.csv("C:/Users/josep/Desktop/St Andrews Course Materials/MT5764 Advanced Data Analysis/data.csv", header = T)

plot(pract1mt5764$impact)
hist(pract1mt5764$impact)
boxplot(pract1mt5764$impact)
ggplot(data = pract1mt5764, aes(x= impact) + geom_histogram)
ggplot(data = pract1mt5764, aes(x = impact)) + geom_histogram()
ggplot(data = pract1mt5764, aes(x = impact, y = year)) + geom_boxplot()

##CONFIDENCE INTERVALS FOR IMPACT=0
results<- matrix(0, nrow=1000, ncol=1)
for(j in 1:1000){
  rowsToUse<-sort(sample(which(pract1mt5764$impact==0),
                         length(which(pract1mt5764$impact==0)), replace=T))
  results[j,]<- mean(pract1mt5764$Nhat[rowsToUse]/pract1mt5764$area[rowsToUse])
}
cisBoot1<-quantile(results, probs=c(0.025,0.975))

##CONFIDENCE INTERVALS FOR IMPACT=1
results<- matrix(0, nrow=1000, ncol=1)
for(j in 1:1000){
  rowsToUse<-sort(sample(which(pract1mt5764$impact==1),
                         length(which(pract1mt5764$impact==1)), replace=T))
  results[j,]<- mean(pract1mt5764$Nhat[rowsToUse]/pract1mt5764$area[rowsToUse])
}
cisBoot2<-quantile(results, probs=c(0.025,0.975))

### PLOT ###

par(mfrow=c(1,1))
plot(pract1mt5764$x.pos, pract1mt5764$y.pos, xlab="X-coordinate", ylab="Y-coordinate",
     main="Transect lines", pch=20, col="lightgrey")
points(pract1mt5764$x.pos, pract1mt5764$y.pos, pch=20, col="blue", cex=log(pract1mt5764$Nhat+1))
par(mfrow=c(1,2))
plot(pract1mt5764$x.pos[pract1mt5764$impact==0], pract1mt5764$y.pos[pract1mt5764$impact==0], xlab="X-coordinate",
     ylab="Y-coordinate", main="Pre impact", pch=20, col="lightgrey")
points(pract1mt5764$x.pos[pract1mt5764$impact==0], pract1mt5764$y.pos[pract1mt5764$impact==0], pch=1, col="blue",
       cex=log(pract1mt5764$Nhat[pract1mt5764$impact==0]+1))
plot(pract1mt5764$x.pos[pract1mt5764$impact==1], pract1mt5764$y.pos[pract1mt5764$impact==1], xlab="X-coordinate",
     ylab="Y-coordinate", main="Post impact", pch=20, col="lightgrey")
points(pract1mt5764$x.pos[pract1mt5764$impact==1], pract1mt5764$y.pos[pract1mt5764$impact==1], pch=1,
       col="blue", cex=log(pract1mt5764$Nhat[pract1mt5764$impact==0]+1))


## MODEL FITTING ##

## MODEL WITH IMPACT AS A COVARIATE ##

glmfit1 <- glm(Nhat ~ impact, offset = log(area), family = poisson, data = pract1mt5764)

glmfitOD1 <- glm(Nhat ~ impact, offset = log(area), family = quasipoisson, data = pract1mt5764)

## MODEL WITH IMPACT, DEPTH, XPOS AND YPOS
## COVARIATES ##

glmfit2 <- glm(Nhat ~ impact + Depth + x.pos + y.pos, offset = log(area),
               family = poisson, data = pract1mt5764)

glmfitOD2 <- glm(Nhat ~ impact + Depth + x.pos + y.pos, offset = log(area),
                 family = quasipoisson, data = pract1mt5764)

## MODEL WITH IMPACT, DEPTH, XPOS, YPOS, IMPACT*XPOS
## AND IMPACT.YPOS AS COVARIATES

glmfit3 <- glm(Nhat ~ impact + Depth + x.pos + y.pos +
                 impact*x.pos + impact*y.pos, 
               offset = log(area), family = poisson,
               data = pract1mt5764)

glmfitOD3 <- glm(Nhat ~ impact + Depth + x.pos + y.pos +
                   impact*x.pos + impact*y.pos, 
                 offset = log(area), family = quasipoisson,
                 data = pract1mt5764)

### ANOVAS ##

Anova(glmfit1)
Anova(glmfitOD1, test= "F")

Anova(glmfit2)
Anova(glmfitOD2, test= "F")

Anova(glmfit3)
Anova(glmfitOD3, test= "F")

### MODEL SELECTION ###

options(na.action = "na.fail")

dredge(glmfit3)

dredge(glmfitOD3)

## the best glmfit3 model is the one with all the covariates included

library(ggplot2)
scatterplotMatrix(glmfit3)
library(GGally)
ggpairs(glmfit3)
ggpairs(pract1mt5764)


subdat <- pract1mt5764[,pract1mt5764$impact == "0"]

plot(pract1mt5764$Nhat, pract1mt5764$impact)


#########practical 2##############
#########DETECTING COLLINEARITY###########


at2 <- attach(pract1mt5764)

#3. 
pairs(cbind(x.pos, y.pos, Depth, impact))

#4.
require(car)
vif(glmfitOD2)

#5.
require(car)
vif(glmfitOD3)

#6. 
xmatrix <- model.matrix(glmfitOD3)
head(xmatrix)

############## RIDGE REGRESSION ##############

#1.
require(glmnet)
xmatrix <- model.matrix(glmfitOD3)
ridge <- glmnet(xmatrix, pract1mt5764$Nhat, family = "poisson", 
                offset = log(area), alpha = 0)
cvridge <- cv.glmnet(xmatrix, pract1mt5764$Nhat, family = "poisson",
                     offset = log(area), alpha=0, nfolds=10)

#2.
par(mfrow=c(1,2))
plot(ridge, xvar = "lambda")
abline(v=log(cvridge$lambda.min))
plot(cvridge)
abline(v=log(cvridge$lambda.min))

#3. 
log(cvridge$lambda)
log(cvridge$lambda.min)

#4. 
plot(ridge, xvar = "lambda", ylim = c(-0.00001, 0.00001))

#5. 
par(mfrow=c(1,1))
cis <- confint(glmfitOD3)
plot(1:7, coef(glmfitOD3), ylim = c(range(cis)),
     xaxt= "n", xlab = "Coefficients", pch=20)
segments(1:7, cis[,1], 1:7, cis[,2] )
ridgecoefs <- coef(ridge) [-2,
                           which(ridge$lambda==cvridge$lambda.min)]
points(1:7, ridgecoefs, col=2, pch=20)
axis(1, at=1:7, names(coef(glmfitOD3)))

#6.
plot(1:7, coef(glmfitOD3), ylim = c(range(cis[-c(1,3),])), xaxt="n",
     xlab="Coefficients", pch=20)
segments(1:7, cis[,1], 1:7, cis[,2] )
points(1:7, ridgecoefs, col=2, pch=20)
axis(1, at=1:7, names(coef(glmfitOD3)))

#7.
plot(1:7, coef(glmfitOD3), ylim = c(range(cis[-c(1,3)])), xaxt="n",
     xlab = "Coefficients", pch=20)
segments(1:7, cis[,1], 1:7, cis[,2] )
points(1:7, ridgecoefs, col=2, pch=20)
axis(1, at=1:7, names(coef(glmfitOD3)))

#8. 
ifelse(cis[,1]<ridgecoefs & ridgecoefs < cis[,2],1,0)

######### LASSO REGRESSION #########
#1.
require(glmnet)
xmatrix <- model.matrix(glmfitOD3)
lasso <- glmnet(xmatrix, pract1mt5764$Nhat, family = "poisson",
                offset = log(area), alpha = 1)
cvlasso <- cv.glmnet(xmatrix, pract1mt5764$Nhat, family = "poisson",
                     offset = log(area), alpha = 1, nfolds = 10)

#2
par(mfrow=c(1,2))
plot(lasso, xvar = "lambda")
abline(v=log(cvlasso$lambda.min))
plot(cvlasso)
abline(v=log(cvlasso$lambda.min))

#3.
par(mfrow=c(1,1))
plot(lasso, xvar = "lambda", ylim=c(-0.00001, 0.00001))

#4.
par(mfrow=c(1,1)) 
plot(1:7, coef(glmfitOD3), ylim=c(range(cis)), 
     xaxt="n", xlab="Coefficients", pch=20) 
segments(1:7,cis[,1], 1:7, cis[,2] ) 
lassocoefs<- coef(lasso)[-2, 
                         which(lasso$lambda==cvlasso$lambda.min)] 
points(1:7, lassocoefs, col=2, pch=20) 
axis(1, at=1:7, names(coef(glmfitOD3))) 

#5. 
plot(1:7, coef(glmfitOD3), ylim=c(range(cis[-c(1,3),])), 
     xaxt="n", xlab="Coefficients", pch=20) 
segments(1:7,cis[,1], 1:7, cis[,2] ) 
points(1:7, lassocoefs, col=2, pch=20) 
axis(1, at=1:7, names(coef(glmfitOD3)))

#6. 
plot(1:7, coef(glmfitOD3), ylim=c(range(cis[-c(1,3,5),])), 
     xaxt="n", xlab="Coefficients", pch=20) 
segments(1:7,cis[,1], 1:7, cis[,2] ) 
points(1:7, lassocoefs, col=2, pch=20) 
axis(1, at=1:7, names(coef(glmfitOD3)))

plot(pairs(pract1mt5764[,1:7], pch=20))
#7. 
lassocoefs
ifelse(cis[,1]<lassocoefs & lassocoefs < cis[,2],1,0)

######### ELASTIC NET REGRESSION ########
#1.
require(glmnet) 
xmatrix<- model.matrix(glmfitOD3) 
enet<- glmnet(xmatrix, pract1mt5764$Nhat, family="poisson", 
              offset=log(area), alpha=0.5)
cvenet<- cv.glmnet(xmatrix, pract1mt5764$Nhat, family="poisson",
                   offset=log(area), nfolds=10, alpha=0.5)

#2.
par(mfrow=c(1,2)) 
plot(enet, xvar="lambda") 
abline(v=log(cvenet$lambda.min)) 
plot(cvenet) 
abline(v=log(cvenet$lambda.min))

#3. 
par(mfrow=c(1,1)) 
plot(enet, xvar="lambda", ylim=c(-0.00001, 0.00001))

#4. 
par(mfrow=c(1,1)) 
plot(1:7, coef(glmfitOD3), ylim=c(range(cis)), 
     xaxt="n", xlab="Coefficients", pch=20) 
segments(1:7,cis[,1], 1:7, cis[,2] ) 
enetcoefs<- coef(enet)[-2, 
                       which(enet$lambda==cvenet$lambda.min)] 
points(1:7, enetcoefs, col=2, pch=20) 
axis(1, at=1:7, names(coef(glmfitOD3)))

#5.
plot(1:7, coef(glmfitOD3), ylim=c(range(cis[-c(1,3),])), 
     xaxt="n", xlab="Coefficients", pch=20) 
segments(1:7,cis[,1], 1:7, cis[,2] ) 
points(1:7, enetcoefs, col=2, pch=20) 
axis(1, at=1:7, names(coef(glmfitOD3)))

#6. 
plot(1:7, coef(glmfitOD3), ylim=c(range(cis[-c(1,3,5),])), 
     xaxt="n", xlab="Coefficients", pch=20) 
segments(1:7,cis[,1], 1:7, cis[,2] ) 
points(1:7, enetcoefs, col=2, pch=20) 
axis(1, at=1:7, names(coef(glmfitOD3)))

#7. 
ifelse(cis[,1]<enetcoefs & enetcoefs< cis[,2],1,0)

#8.
lmcoefs<- coef(glmfitOD3) 
for(i in 1:length(lmcoefs)){ 
  plot(1:4, c(lmcoefs[i], ridgecoefs[i], lassocoefs[i], enetcoefs[i]), 
       ylim=range(c(lmcoefs[i], ridgecoefs[i], lassocoefs[i], enetcoefs[i], cis[i,])), 
       pch=20, col=c(1:4),xlab="Modelling method", 
       xaxt="n", ylab="Coefficients", main=names(coef(glmfitOD3))[i]) 
  axis(1, at=1:4, c("LM", "Ridge", "LASSO", "Enet"))
  segments(1,cis[i,1],1 , cis[i,2]) 
  abline(h=c(cis[i,1], cis[i,2]))}


################# PRACTICAL 3 #######################

#1.
testVals<- rnorm(50)

#2. 
sign(testVals) 
plot(sign(testVals), type="l") 

#3.
#plot residuals for the working model in order: 
plot(sign(residuals(glmfitOD3, type="pearson")[1:800]), type="l")

#4.
require(lawstat) 
runs.test(residuals(glmfitOD3, type="pearson"))

#5. 
require(car) 
residualPlots(glmfitOD3, quadratic=T, type = "pearson", 
              ylim=c(-20,20),terms = ~ Depth, fitted=FALSE) 

#6. 
pract1mt5764$Depth[order(pract1mt5764$Depth)][1:100] 
par(mfrow=c(1,2)) 
plot(sign(residuals(glmfitOD3, type="pearson")[order(pract1mt5764$Depth)])[1:100], 
     type="l", main="100 Residuals in Depth order", ylab="Pearsons Residuals") 
plot(sign(rnorm(100)), type="l", main="Random values") 

#7. 
require(lawstat)
runs.test(residuals(glmfitOD3, type="pearson")[order(pract1mt5764$Depth)])

#8. 
residualPlots(glmfitOD3, quadratic=T, type = "pearson", 
              ylim=c(-20,20),terms = ~ x.pos, fitted=FALSE) 

#9. 
require(lawstat) 
runs.test(residuals(glmfitOD3, type="pearson")[order(pract1mt5764$x.pos)])

#10.
residualPlots(glmfitOD3, quadratic=T, type = "pearson", 
              ylim=c(-20,20),terms = ~ y.pos, fitted=FALSE)

#11.
require(lawstat) 
runs.test(residuals(glmfitOD3, type="pearson")[order(pract1mt5764$y.pos)])

