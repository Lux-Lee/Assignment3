##Picloram##
Picloram_Data$Proportion <- Picloram_Data$Number.Died/Picloram_Data$Total.Exposed
plot(Picloram_Data$Dose,Picloram_Data$Proportion, xlab="Dose(unit)", 
     ylab="Proportion of death", main="Proportion of death vs Dose")
Picloram.model1 <-lm(Proportion~Dose, data=Picloram_Data)
summary(Picloram.model1)
abline(Picloram.model1, col="blue")

Picloram.model2 <-glm(cbind(Number.Died, Number.Survived)~Dose, data=Picloram_Data,
                     family=binomial)
summary(Picloram.model2)
Dose.value <- seq(-1,5,0.005)
Pro.value <- ilogit(Picloram.model2$coefficients[1]+Picloram.model2$coefficients[2]*Dose.value)

dose.p(Picloram.model2,p=0.5)
dose.p(Picloram.model2,p=0.9)
points(1.36479, 0.5, pch=15, col="red")
points(2.54078, 0.9, pch=15, col="red")

predict(Picloram.model1, newdata=data.frame(Dose=(2.2)))
ilogit(Picloram.model2$coefficients[1]+Picloram.model2$coefficients[2]*2.2)
predict(Picloram.model3, newdata=data.frame(Dose=2.2))
pnorm(0.7132237)

##Picloram Probit##
plot(Picloram_Data$Dose,Picloram_Data$Proportion, xlab="Dose(unit)", 
     ylab="Proportion of death", main="Proportion of death vs Dose")
Picloram.model3 <-glm(cbind(Number.Died, Number.Survived)~Dose, data=Picloram_Data,
                      family=binomial(link=probit))
summary(Picloram.model3)
Pro.value2 <- pnorm(Picloram.model3$coefficients[1]+Picloram.model3$coefficients[2]*Dose.value)
lines(Dose.value, Pro.value2, col="green")
dose.p(Picloram.model3, p=0.5)
dose.p(Picloram.model3, p=0.9)

points(1.429759,0.5, col="green")
points(2.811824,0.9, col="green")

legend("bottomright", legend=c("Logit","Probit"), pch=c(15,1), lty=c(1,1), 
       col=c("red","green"))

##West Nile Virus##
WNV.model1 <-glm(No.HorseCase~offset(log(No.Farms))+No.BirdCase+Area+Population, family=poisson, data=WNV_Data)
summary(WNV.model1)
WNV.model2 <-glm(No.HorseCase~offset(log(No.Farms))+No.BirdCase+Area, family=poisson, data=WNV_Data)
summary(WNV.model2)
anova(WNV.model2,WNV.model1, test="Chisq")
###Metribuzin###
##SSD##
Metribuzin_fit<-ssd_fit_dists(Metribuzin_Species,dists=c("lnorm", "gamma", "invpareto",
                                                   "llogis", "lgumbel", "weibull", "gompertz"))
Metribuzin_gof<-ssd_gof(Metribuzin_fit)
Metribuzin_gof[order(Metribuzin_gof$delta),]
Metribuzin_fits<-ssd_fit_dists(Metribuzin_Species,dists=c("lgumbel"))
Metribuzin_fits
set.seed(99)
Metribuzin_pred<-predict(Metribuzin_fits,ci=TRUE)
Metribuzin_pred
Metplot<-ssd_plot(Metribuzin_Species,Metribuzin_pred,xlab="Concentration (mg/L)",
                  ylab="Centile",ribbon=TRUE, col="Group")
print(Metplot)
autoplot(Metribuzin_fit, delta=100)
##EED##
Metribuzin_dist<- ssd_fit_dists(Metribuzin_Exposure,dists=c("lnorm", "gamma",
                                                            "invpareto", "llogis", "lgumbel", "weibull", "gompertz"))
autoplot(Metribuzin_dist, delta=10000)
Metribuzin_E_pred<-predict(Metribuzin_Exposure,ci=TRUE, nboot = 1000,
                            parametric = FALSE)
Metplot2<- ggplot(Metribuzin_pred, aes_string(x = "est")) +
  geom_xribbon(aes_string(xmin = "lcl", xmax = "ucl", y = "percent/100"), alpha = 0.2) +
  geom_line(aes_string(y = "percent/100")) +
  geom_ssdpoint(data = Metribuzin_Species, aes_string(x = "Conc",color="Group")) +
  scale_y_continuous("Centile", labels = scales::percent) +
  expand_limits(y = c(0, 1)) +
  xlab("Concentration (mg/L)")
Metplot2
