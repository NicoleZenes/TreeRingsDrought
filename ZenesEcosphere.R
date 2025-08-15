#Ecosphere
#Article
#ECS25-0062 - Complex trajectories of tree growth in the southwestern United States after severe drought
#Nicole Zenes, Leander D. L. Anderegg, Kiona Ogle, Drew M. P. Peltier, and William R. L.Anderegg

#Libraries for data analysis
library(MuMIn)
library(dplyr)
library(lme4)
library(tidyr)
library(reshape2)
library(data.table)

#Libraries for visualization
library(ggpattern)
library(viridis)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(plotrix)

#Internal functions for easy of analysis
#function to remove outliers from dataset
outlierReplace = function(dataframe, cols, rows, newValue = NA) {
  if (any(rows)) {
    set(dataframe, rows, cols, newValue)
  }
}
#number extract function for plot names (AZ435 vs. AZ_PIED_435)
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 


###########################################################################################
period <- as.data.frame(c("aet","def","soil","swe","PDSI","vpd","tmax","ppt"))

###-----single pixel file-----#
single <- read.csv("/Users/nicolezenes/Library/Mobile Documents/com~apple~CloudDocs/FIA/PlotLocationFuzzed_NDVI_singlePixel.csv")
zscores <- read.csv("/Users/nicolezenes/Library/Mobile Documents/com~apple~CloudDocs/FIA/NDVI_singlePixel_Zscore_test.csv")

###########################################################################################
#Doing nested models by species with individual within site 
#source("/Users/nicolezenes/Library/Mobile Documents/com~apple~CloudDocs/FIA/Terraclim/terraclimR.R")
#change above to read in files instead of pulling off the web all the time 
##Read in climate files to average per year
for (p in (1:8)) {
  var <- as.character(period[p,1])
  
  #Assign variable to variable name
  assign(paste0(var, "climate"), read.csv(paste(
    "/Users/nicolezenes/Library/Mobile Documents/com~apple~CloudDocs/FIA/Terraclim/", var,"full.csv")))
  
  ##average per year
  average <- as.data.frame(matrix(0, ncol=57, nrow=34))
  j <- 3
  for (i in (seq(3, 674, by=12))) {
    x <- as.data.frame(assign(paste0(var, "climate"), read.csv(paste(
      "/Users/nicolezenes/Library/Mobile Documents/com~apple~CloudDocs/FIA/Terraclim/", var,"full.csv"))))
    average[,j] <- rowMeans(x[,i:(i+11)])
    j <- j + 1
  }
  average[,1:2] <- x[,1:2]
  assign(paste0(var, "average"), average)
  
}


setwd("/Users/nicolezenes/Library/Mobile Documents/com~apple~CloudDocs/FIA/AllRW_data_bill/")
fullrrwdrew <- read.csv("ring_dataID.csv") 
corelist <- read.csv("corelist.csv", header=F)

xx <- as.data.frame(unique(fullrrwdrew$file)) #get unique plot names from master datasheet 
period <- as.data.frame(c("aet","def","soil","swe","PDSI","vpd","tmax","ppt")) #climate variables 
state <- as.data.frame(c("AZ_", "UT_", "CO_")) #states
species <- as.data.frame(c('POTR','PIED', 'PIEN', 'JUOS', 'PIPO')) #species

#subset data from master sheet for each plot by matching state and plot number with climate sheets
for (i in 1:nrow(xx)) {
  num <- as.character(xx[i,1])
  uniplot <- fullrrwdrew[(str_detect(fullrrwdrew[,2], num)),]
  #find combination of species, state, and plot number that is correct 
  uniplot1 <- uniplot
  
  #cast uniplot to sepearate unique cores
  uniplot1
  cast <- dcast(uniplot1, year ~ UniqueID)
  #separate into 1960-1989, and 2000-2013
  try( castbaseline <- cast[ which(grepl(1960, cast$year)):which(grepl(1989, cast$year)),])
  cast2013 <- cast[ which(grepl(2000, cast$year)):which(grepl(2013, cast$year)),]
  
  #subset climate variables for plot for the 60's, 70's, 80's baseline 
  for (p in 1:nrow(period)) {
    var <- as.character(period[p,1])
    #get dataframe for assigned climate variable 
    fullvar <- get(paste0(var,"average"))
    fullvar <- cbind(corelist[1:34,],fullvar)
    
    #extract row that matches plot above
    plotvar <- fullvar[(str_detect(fullvar[,1], num)),]
    #extract 1960-1989 and 2000-2013 climate data and assign to variable
    assign(paste0("plot",var, "baseline"), t(plotvar[,6:35]))
    assign(paste0("plot",var, "2013"), t(plotvar[,46:59]))
  }
  
  #do long form with climate variables for 1960-1989
  castbaseline <- cbind(castbaseline,(plotaetbaseline),(plotdefbaseline),(plotsoilbaseline),
                        (plotswebaseline),(plotPDSIbaseline),(plotvpdbaseline),
                        (plottmaxbaseline),
                        (plotpptbaseline))
  #Add column names for each of the variables (hard coded for these names and number of variables - fix with period above)
  colnames(castbaseline)[(length(castbaseline)-7):length(castbaseline)] <- c("aet","def","soil","swe","PDSI","vpd","tmax","ppt")
  meltbaseline <- melt(data = castbaseline, id.vars = c("year","aet","def","soil","swe","PDSI","vpd","tmax","ppt"))
  
  #add column with plot number and species repeated
  plotID <- rep((paste0(num)), nrow(meltbaseline))
  meltbaseline <- cbind(meltbaseline, plotID)
  
  #assign name to file with state, plot number, and species ID
  assign(paste0(num), meltbaseline)
  
  #do long form with climate variables for 2000-2013
  cast2013 <- cbind(cast2013,(plotaet2013),(plotdef2013),(plotsoil2013),
                    (plotswe2013),(plotPDSI2013),(plotvpd2013),
                    (plottmax2013),
                    (plotppt2013))
  
  #Add column names for each of the variables (hard coded for these names and number of variables - fix with period above)
  colnames(cast2013)[(length(cast2013)-7):length(cast2013)] <- c("aet","def","soil","swe","PDSI","vpd","tmax","ppt")
  melt2013 <- melt(data = cast2013, id.vars = c("year","aet","def","soil","swe","PDSI","vpd","tmax","ppt"))
  
  #add column with plot number and species repeated
  plotID <- rep((paste0(num)), nrow(melt2013))
  melt2013 <- cbind(melt2013, plotID)
  
  #assign name to file with state, plot number, and species ID 
  assign(paste0(num, "2013"), melt2013)
}


###########################POTR#####################
#####combine into species full 1960-1989
#initialize dataframe for potr
POTRclimateall <- get(paste0(corelist[23,]))
#Hard coded - which rows are potr
list <- as.data.frame(c(24:34))
#get climate files and row bind long form
for (c in 1:nrow(list)){
  POTRclimateall <- rbind(POTRclimateall,get(paste0(corelist[(list[c,]),])))
}
POTRclimateall <- na.omit(POTRclimateall)

model = lmer(value ~ year + ppt + vpd + def + tmax + aet + soil + PDSI + swe + (1|variable) + (1|plotID), 
             data = POTRclimateall, REML=TRUE)
summary(model)
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)  2.145e+01  1.404e+00  4.620e+03  15.281  < 2e-16 ***
#   year        -1.027e-02  7.074e-04  4.649e+03 -14.523  < 2e-16 ***
#   ppt         -4.370e-03  1.019e-03  4.635e+03  -4.290 1.83e-05 ***
#   vpd         -1.548e-01  2.616e-01  4.262e+03  -0.592  0.55417    
# def         -1.058e-03  2.143e-03  4.630e+03  -0.493  0.62172    
# tmax        -2.919e-02  1.191e-02  3.364e+03  -2.452  0.01427 *  
#   aet          7.682e-03  2.776e-03  4.635e+03   2.767  0.00567 ** 
#   soil        -1.959e-03  9.467e-04  4.306e+03  -2.069  0.03856 *  
#   PDSI         1.290e-02  6.050e-03  3.434e+03   2.133  0.03299 *  
#   swe          8.739e-04  4.092e-04  4.376e+03   2.136  0.03274 *  
AIC(model) #4416
r.squaredGLMM(model)
# R2m       R2c
# [1,] 0.06536307 0.5707313

model1 = lmer(value ~ year + ppt + tmax + aet + PDSI + soil + swe + (1|variable) + (1|plotID), data = POTRclimateall, REML=FALSE)
summary(model1)
qqnorm(resid(model1))
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)  2.106e+01  1.364e+00  4.554e+03  15.440  < 2e-16 ***
#   year        -1.012e-02  6.955e-04  4.608e+03 -14.557  < 2e-16 ***
#   ppt         -4.070e-03  9.480e-04  4.644e+03  -4.293 1.80e-05 ***
#   tmax        -3.669e-02  9.830e-03  1.793e+03  -3.732 0.000196 ***
#   aet          8.891e-03  1.556e-03  4.630e+03   5.715 1.16e-08 ***
#   PDSI         1.389e-02  5.969e-03  3.603e+03   2.328 0.019987 *  
#   soil        -1.970e-03  9.438e-04  4.303e+03  -2.087 0.036919 *  
#   swe          7.839e-04  4.010e-04  4.505e+03   1.955 0.050695 .  

AIC(model1) #4666
r.squaredGLMM(model1)
# R2m       R2c
#[1,] 0.06032489 0.5520439

rand(model1)
# ANOVA-like table for random-effects: Single term deletions
# 
# Model:
#   value ~ year + ppt + aet + PDSI + swe + (1 | plotID:variable) + 
#   (1 | plotID)
# npar  logLik    AIC     LRT Df Pr(>Chisq)    
# <none>                   9 -3113.7 6245.4                          
# (1 | plotID:variable)    8 -4063.9 8143.8 1900.38  1  < 2.2e-16 ***
#   (1 | plotID)             8 -3143.7 6303.4   59.99  1  9.556e-15 ***


#initialize dataframe for potr for 2000-2013
POTRclimate2013 <- get(paste0(corelist[23,],"2013"))
#Hard coded - which rows are potr
list <- as.data.frame(c(24:34))
#get climate files and row bind long form
for (c in 1:nrow(list)){
  POTRclimate2013 <- rbind(POTRclimate2013,get(paste0(corelist[(list[c,]),],"2013")))
}

#get predicted values for 2000-2013 year + ppt + tmax + aet + PDSI + soil + swe
pred <- summary(model1)$coefficients[1,1] + summary(model1)$coefficients[2,1] * POTRclimate2013$year + 
  summary(model1)$coefficients[3,1] * POTRclimate2013$ppt +  
  summary(model1)$coefficients[4,1] * POTRclimate2013$tmax + 
  summary(model1)$coefficients[5,1] * POTRclimate2013$soil +
  summary(model1)$coefficients[6,1] * POTRclimate2013$PDSI +
  summary(model1)$coefficients[7,1] * POTRclimate2013$swe

POTRclimate2013 <- cbind(POTRclimate2013, pred)
diff <- POTRclimate2013$value - POTRclimate2013$pred
POTRclimate2013 <- cbind(POTRclimate2013, diff)

#plots looking at the data, remove outliers 
boxplot(POTRclimate2013$value, POTRclimate2013$diff)
plot(POTRclimate2013$value, POTRclimate2013$pred)
plot(model1)

boxplot(diff ~ plotID*year, data=POTRclimate2013)
abline(h=0)
boxplot(POTRclimate2013$value)$out
min(boxplot(POTRclimate2013$value)$out)
qplot(data=POTRclimate2013, x=value)

outlierReplace(POTRclimate2013, "value", which(POTRclimate2013$value > 1.558), NA)
qplot(data=POTRclimate2013, x=value)
withoutpotr <- na.omit(POTRclimate2013)


#line plots with means 
means <- aggregate(withoutpotr[, 14], list(withoutpotr$year, withoutpotr$plotID), mean)
se <- aggregate(diff ~ year + plotID, data = withoutpotr,FUN=function(x) c(mean=mean(x),se=std.error(x)))
means <- cbind(means, (se$diff[,1]-se$diff[,2]), (se$diff[,1]+se$diff[,2]))

#################
####---2024 final paper edits potr SI figures (A)---####
#line plots with means 
means <- aggregate(withoutpotr[, 14], list(withoutpotr$year, withoutpotr$plotID), mean)
se <- aggregate(diff ~ year + plotID, data = withoutpotr,FUN=function(x) c(mean=mean(x),se=std.error(x)))
means <- cbind(means, (se$diff[,1]-se$diff[,2]), (se$diff[,1]+se$diff[,2]))

potrrrwse <- ggplot(data = means, aes(x=Group.1, y=x, color=Group.2)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin=(means[,4]), ymax=(means[,5])), width=.2) + ylab("Observed - predicted RRW") +
  scale_color_hue(name = "PlotID", labels = c("AZ1060", "AZ 940", "CO 166", "CO 25", "CO 25A", "CO 3014", "CO 37",
                                              "CO 55", "CO 63", "UT 74", "UT 297", "UT 73")) + geom_hline(yintercept = 0) +
  xlab("Year") + labs(tag = "A")
potrrrwse
ggsave(
  "sifigure_potr_rrw_indiv",
  plot = potrrrwse,
  device = png)

plots <- data.frame (name = c("AZ1060", "AZ 940", "CO 166", "CO 25", "CO 25A", "CO 3014", "CO 37",
                              "CO 55", "CO 63", "UT 74", "UT 297", "UT 73" ), 
                     value = c(0.161,0.0755, 0.0917,0.111,0.2038,0.619,0.0685, 0.01815, 0.124, 0.1889, 0.1579, 0.0188))

potrmort <- ggplot(data = plots, mapping = aes(x = name, y =value)) +
  geom_bar(stat="identity", aes(color=name, fill=name), width=.5) + ylab("Fraction of mortality reported 2014")+
  xlab("Plot") + ylim(c(0,1)) +  theme(legend.title = element_blank()) + labs(tag="A") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position="none")
potrmort
ggsave("sifigure_potr_mort",
       plot = potrmort,
       device = png)


#point plot of means with se bars
dodge <- position_dodge(.3)
pp <- ggplot(data = means, aes(x=Group.1, y=x, color=Group.2)) + geom_point(position=dodge)
pp <- pp +
  geom_errorbar(aes(ymin=(means[,4]), ymax=(means[,5])), width=1, position=dodge) + geom_hline(yintercept=0)

#Single pixel 
potrsingle <- single[c(21,22,23,24,25,26,27,28,29,30,31,32), c(1, 4:18)]
potrsingle <- gather(potrsingle, year, measurement, X2000:X2014, factor_key = TRUE)
potrndvi <- ggplot(data = potrsingle, aes(x=year, y=measurement, color=ID, group=ID)) + geom_point() + geom_line() +
  xlab("Year") + ylab("NDVI") + scale_x_discrete(breaks=c("X2000", "X2005", "X2010"), labels = c("2000", "2005",
                                                                                                 "2010"))
potrzscore <- zscores[c(21,22,23,24,25,26,27,28,29,30,31,32), c(1, 4:18)]
potrzscore <- gather(potrzscore, year, measurement, X2000:X2014, factor_key = TRUE)
potrndviz <- ggplot(data = potrzscore, aes(x=year, y=measurement, color=ID, group=ID)) + geom_point() + geom_line() +
  xlab("Year") + ylab("NDVI Z-Score") + scale_x_discrete(breaks=c("X2000", "X2005", "X2010"), labels = c("2000", "2005",
                                                                                                         "2010")) +
  ggtitle("Potr")+
  scale_colour_manual(values = c("#440154FF", "#482173FF", "#433E85FF",
                                 "#38598CFF", "#2D708EFF", "#25858EFF",
                                 "#1E9B8AFF", "#2BB07FFF", "#51C56AFF",
                                 "#85D54AFF", "#C2DF23FF", "#FDE725FF")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())




###########################JUOS#####################
#####combine into species full 1960-1989
#initialize dataframe for juos - 5, 7, 11, 14
juosclimateall <- get(paste0(corelist[5,]))
#Hard coded - which rows are juos
list <- as.data.frame(c(7,11,14))
#get climate files and row bind long form
for (c in 1:nrow(list)){
  juosclimateall <- rbind(juosclimateall,get(paste0(corelist[(list[c,]),])))
}

#random intercept by plot 
model = lmer(value ~ year + ppt + vpd + def + tmax + aet + soil + PDSI + swe + (1|variable) + (1|plotID), 
             data = juosclimateall, REML=TRUE)
summary(model)

# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)  2.017e+01  2.296e+00  1.582e+03   8.786  < 2e-16 ***
#   year        -1.021e-02  1.156e-03  1.647e+03  -8.835  < 2e-16 ***
#   ppt          3.749e-03  1.911e-03  1.619e+03   1.962  0.04990 *  
#   vpd         -6.003e-01  4.934e-01  1.192e+03  -1.217  0.22397    
# def          4.317e-03  3.912e-03  1.641e+03   1.103  0.27004    
# tmax         5.625e-02  2.217e-02  1.571e+03   2.537  0.01128 *  
#   aet         -7.364e-03  4.724e-03  1.647e+03  -1.559  0.11920    
# soil         6.094e-03  2.068e-03  1.563e+03   2.947  0.00326 ** 
#   PDSI         7.125e-02  1.002e-02  1.529e+03   7.107 1.81e-12 ***
#   swe         -2.789e-03  9.190e-04  1.587e+03  -3.034  0.00245 ** 
AIC(model) #1223

#ppt, tmax, soil, PDSIm and swe year significant
model1 = lmer(value ~ year + ppt + tmax + soil + PDSI + swe + (1|variable) + (1|plotID), 
              data = juosclimateall, REML=TRUE)
summary(model1)
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)  1.769e+01  1.968e+00  1.594e+03   8.990  < 2e-16 ***
#   year        -9.043e-03  1.021e-03  1.649e+03  -8.860  < 2e-16 ***
#   ppt          3.564e-04  1.387e-03  1.648e+03   0.257 0.797286    
# tmax         4.322e-02  1.452e-02  1.260e+03   2.978 0.002962 ** 
#   soil         7.092e-03  1.988e-03  1.633e+03   3.567 0.000372 ***
#   PDSI         4.960e-02  8.179e-03  1.598e+03   6.065 1.65e-09 ***
#   swe         -1.157e-03  7.840e-04  1.650e+03  -1.476 0.140166   
AIC(model1) #1215

model2 = lmer(value ~ year + tmax + soil + PDSI + (1|variable) + (1|plotID), 
              data = juosclimateall, REML=TRUE)
r.squaredGLMM(model2)

summary(model2)
# 
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)  1.723e+01  1.946e+00  1.595e+03   8.856  < 2e-16 ***
#   year        -8.840e-03  1.013e-03  1.651e+03  -8.726  < 2e-16 ***
#   tmax         4.582e-02  1.431e-02  1.241e+03   3.202  0.00140 ** 
#   soil         6.975e-03  1.987e-03  1.635e+03   3.510  0.00046 ***
#   PDSI         4.605e-02  6.985e-03  1.622e+03   6.592 5.84e-11 ***
AIC(model2) #1190


#initialize dataframe for juos for 2000-2013
juosclimate2013 <- get(paste0(corelist[5,],"2013"))
#Hard coded - which rows are juos
list <- as.data.frame(c(7,11,14))
#get climate files and row bind long form
for (c in 1:nrow(list)){
  juosclimate2013 <- rbind(juosclimate2013,get(paste0(corelist[(list[c,]),],"2013")))
}

#get predicted values for 2000-2013
pred <- summary(model2)$coefficients[1,1] + summary(model2)$coefficients[2,1] * juosclimate2013$year + 
  summary(model2)$coefficients[3,1] * juosclimate2013$tmax +  
  summary(model2)$coefficients[4,1] * juosclimate2013$soil +  
  summary(model2)$coefficients[5,1] * juosclimate2013$PDSI 
juosclimate2013 <- cbind(juosclimate2013, pred)
diff <- juosclimate2013$value - juosclimate2013$pred
juosclimate2013 <- cbind(juosclimate2013, diff)

#plots looking at the data, remove outliers
boxplot(juosclimate2013$value, juosclimate2013$diff)
plot(juosclimate2013$value, juosclimate2013$pred)
plot(model1)

boxplot(diff ~ year, data=juosclimate2013)
abline(h=0)
min(boxplot(juosclimate2013$value)$out)

outlierReplace(juosclimate2013, "value", which(juosclimate2013$value > 1.312), NA)
qplot(data=juosclimate2013, x=value)

boxplot(juosclimate2013$value)$out
#[1] 1.283 1.293 1.270
outlierReplace(juosclimate2013, "value", which(juosclimate2013$value > 1.26), NA)

boxplot(juosclimate2013$diff)$out
boxplot(juosclimate2013$pred)$out
withoutjuos <- na.omit(juosclimate2013)

#line plots with means 
means <- aggregate(withoutjuos[, 14], list(withoutjuos$year, withoutjuos$plotID), mean)
se <- aggregate(diff ~ year + plotID, data = withoutjuos,FUN=function(x) c(mean=mean(x),se=std.error(x)))
means <- cbind(means, (se$diff[,1]-se$diff[,2]), (se$diff[,1]+se$diff[,2]))

#Single pixel 
juossingle <- single[c(5,7,11,14), c(1, 4:18)]
juossingle <- gather(juossingle, year, measurement, X2000:X2014, factor_key = TRUE)
juosndvi <- ggplot(data = juossingle, aes(x=year, y=measurement, color=ID, group=ID)) + geom_point() + geom_line() +
  xlab("Year") + ylab("NDVI Z-score") + scale_x_discrete(breaks=c("X2000", "X2005", "X2010"), labels = c("2000", "2005","2010"))
                                                                                                         
                                                                                                         
#################
####---2024 final paper edits juos SI figures (B)---####
#line plots with means 
means <- aggregate(withoutjuos[, 14], list(withoutjuos$year, withoutjuos$plotID), mean)
se <- aggregate(diff ~ year + plotID, data = withoutjuos,FUN=function(x) c(mean=mean(x),se=std.error(x)))
means <- cbind(means, (se$diff[,1]-se$diff[,2]), (se$diff[,1]+se$diff[,2]))

juosrrwse <- ggplot(data = means, aes(x=Group.1, y=x, color=Group.2)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin=(means[,4]), ymax=(means[,5])), width=.2) + xlab("Year") + ylab("Observed - Predicted RRW") +
  scale_color_hue(name = "PlotID", labels = c("AZ 940", "CO 137","CO 31", "CO 95")) +  guides(color=guide_legend("PlotID")) +
  geom_hline(yintercept=0) + labs(tag="B") + geom_hline(yintercept = 0)
juosrrwse

ggsave(
  "sifigure_juos_rrw_indiv",
  plot = juosrrwse,
  device = png)

juosplots <- data.frame (name = c("AZ 940", "CO 137","CO 31", "CO 95"), 
                         value = c(.73, .66,0, 0))

juosmort <- ggplot(data = juosplots, mapping = aes(x = name, y =value)) +
  geom_bar(stat="identity", aes(color=name, fill=name), width=.5) + ylab("Fraction of FIA mortality reported 2014") +
  xlab("Plot") + ylim(c(0,1)) +theme(legend.title = element_blank()) + labs(tag="B") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position="none")
juosmort

ggsave(
  "sifigure_juos_mort",
  plot = juosmort,
  device = png)

#point plot of means with se bars
dodge <- position_dodge(.3)
pp <- ggplot(data = means, aes(x=Group.1, y=x, color=Group.2)) + geom_point(position=dodge)
pp
pp +
  geom_errorbar(aes(ymin=(means[,4]), ymax=(means[,5])), width=1, position=dodge)

#separate boxplot per plot across all years 
p1 <- ggplot(without, aes(x=year, y=diff, fill=plotID, group=year)) + geom_boxplot() +facet_wrap(~plotID)
p1
#boxplots for each plot aggregated over years
p2 <- ggplot(without, aes(x=plotID, y=diff, fill=plotID)) + geom_boxplot()
p2


###########################PIED#####################
#####combine into species full 1960-1989
#initialize dataframe for pied - 5, 7, 11, 14
piedclimateall <- get(paste0(corelist[5,]))
#Hard coded - which rows are juos plus two correct pied ones - error in the pied ones for same plots
list <- as.data.frame(c(7,11,14,4,18))
#get climate files and row bind long form
for (c in 1:nrow(list)){
  piedclimateall <- rbind(piedclimateall,get(paste0(corelist[(list[c,]),])))
}
#rename plots
piedclimateall$plotID <- gsub('JUOS', 'PIED', piedclimateall$plotID)
model = lmer(value ~ year + ppt + vpd + def + tmax + aet + soil + PDSI + swe + (1|variable) + (1|plotID), 
             data = piedclimateall, REML=TRUE)
summary(model)

# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)  1.779e+01  1.716e+00  2.338e+03  10.372  < 2e-16 ***
#   year        -8.902e-03  8.740e-04  2.472e+03 -10.185  < 2e-16 ***
#   ppt          3.410e-03  1.331e-03  2.487e+03   2.561  0.01049 *  
#   vpd         -7.248e-01  3.435e-01  1.731e+03  -2.110  0.03498 *  
#   def          1.574e-03  2.743e-03  2.400e+03   0.574  0.56608    
# tmax         5.968e-02  1.827e-02  1.898e+03   3.266  0.00111 ** 
#   aet         -4.648e-03  3.273e-03  2.465e+03  -1.420  0.15574    
# soil         4.923e-03  1.582e-03  2.304e+03   3.113  0.00187 ** 
#   PDSI         5.752e-02  7.620e-03  2.132e+03   7.548 6.52e-14 ***
#   swe         -1.788e-03  7.165e-04  2.418e+03  -2.496  0.01264 *  
AIC(model) #1469

model1 = lmer(value ~ year + ppt + vpd + tmax + soil + PDSI + swe + (1|variable) + (1|plotID), 
              data = piedclimateall, REML=TRUE)
summary(model1)

# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)  1.654e+01  1.628e+00  2.363e+03  10.160  < 2e-16 ***
#   year        -8.308e-03  8.391e-04  2.481e+03  -9.901  < 2e-16 ***
#   ppt          1.058e-03  1.041e-03  2.354e+03   1.017 0.309322    
# vpd         -4.511e-01  2.946e-01  9.275e+02  -1.531 0.126062    
# tmax         5.080e-02  1.765e-02  1.947e+03   2.879 0.004037 ** 
#   soil         5.511e-03  1.531e-03  2.414e+03   3.599 0.000325 ***
#   PDSI         4.615e-02  6.253e-03  2.301e+03   7.380  2.2e-13 ***
#   swe         -8.593e-04  6.518e-04  2.383e+03  -1.318 0.187504   
AIC(model1) #1454

model2 = lmer(value ~ year + tmax + soil + PDSI + (1|variable) + (1|plotID), 
              data = piedclimateall, REML=TRUE)
r.squaredGLMM(model2)
summary(model2)
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)  1.505e+01  1.476e+00  2.409e+03  10.194  < 2e-16 ***
#   year        -7.565e-03  7.725e-04  2.487e+03  -9.793  < 2e-16 ***
#   tmax         3.098e-02  1.076e-02  1.412e+03   2.879  0.00405 ** 
#   soil         4.673e-03  1.502e-03  2.485e+03   3.110  0.00189 ** 
#   PDSI         4.748e-02  5.456e-03  2.353e+03   8.702  < 2e-16 ***
AIC(model2) #1434


#initialize dataframe for pied for 2000-2013
piedclimate2013 <- get(paste0(corelist[5,],"2013"))
#Hard coded - which rows are juos plus two correct pied ones - error in the pied ones for same plots
list <- as.data.frame(c(7,11,14,4,18))
#get climate files and row bind long form
for (c in 1:nrow(list)){
  piedclimate2013 <- rbind(piedclimate2013,get(paste0(corelist[(list[c,]),],"2013")))
}
piedclimate2013$plotID <- gsub('JUOS', 'PIED', piedclimate2013$plotID)
pred <- summary(model2)$coefficients[1,1] + summary(model2)$coefficients[2,1] * piedclimate2013$year + 
  summary(model2)$coefficients[3,1] * piedclimate2013$tmax +  
  summary(model2)$coefficients[4,1] * piedclimate2013$soil +
  summary(model2)$coefficients[5,1] * piedclimate2013$PDSI 


piedclimate2013 <- cbind(piedclimate2013, pred)
diff <- piedclimate2013$value - piedclimate2013$pred
piedclimate2013 <- cbind(piedclimate2013, diff)

#plots looking at the data, remove major outliers 
boxplot(piedclimate2013$value, piedclimate2013$diff)
plot(piedclimate2013$value, piedclimate2013$pred)
plot(model2)

boxplot(diff ~ year, data=piedclimate2013)
abline(h=0)

min(boxplot(piedclimate2013$value)$out) #1.442
qplot(data=piedclimate2013, x=value)

outlierReplace(piedclimate2013, "value", which(piedclimate2013$value > 1.441), NA)
qplot(data=piedclimate2013, x=value)

boxplot(piedclimate2013$value)$out
min(boxplot(piedclimate2013$value)$out) 

outlierReplace(piedclimate2013, "value", which(piedclimate2013$value > 1.327), NA)



##########################
###---2024 final paper edits pied SI figures (C)---####
#line plots with means 
means <- aggregate(withoutpied[, 14], list(withoutpied$year, withoutpied$plotID), mean)
se <- aggregate(diff ~ year + plotID, data = withoutpied,FUN=function(x) c(mean=mean(x),se=std.error(x)))
means <- cbind(means, (se$diff[,1]-se$diff[,2]), (se$diff[,1]+se$diff[,2]))

peidrrwse <- ggplot(data = means, aes(x=Group.1, y=x, color=Group.2)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin=(means[,4]), ymax=(means[,5])), width=.2) + xlab("Year") + ylab("Observed - Predicted RRW") +
  scale_color_hue(labels =  c("AZ 435", "AZ 940", "CO 137","CO 31", "CO 95", "UT 228")) +  guides(color=guide_legend("PlotID")) +
  geom_hline(yintercept=0) + labs(tag="C")
peidrrwse
ggsave("sifigure_pied_rrw_indv",plot = peidrrwse,device = png)

plots <- data.frame (name = c("AZ 435", "UT 228", "AZ 940", "CO 137","CO 31", "CO 95"), 
                     value = c(.55, .2395, .73, .66,0, 0))

#plots <- cbind(plots, piedtrader$deadID)
peidmort <- ggplot(data = plots, mapping = aes(x = name, y =value), pattern = deadID) +
  geom_bar(stat="identity", aes(color=name, fill=name), width=.5) + ylab("Fraction of mortality reported 2014") +
  xlab("Plot") + ylim(c(0,1)) + theme(legend.title = element_blank()) + labs(tag="C") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position="none")
ggsave("sifigure_pied_mort",plot = peidmort,device = png)

#jitter plot points over time with mean plotted as a line
p <- ggplot(withoutpied, aes(x=year, y=diff,  color=plotID)) + 
  geom_jitter(position=position_jitter(0.2)) + 
  geom_line(data = means, aes(x=Group.1, y=x, color=Group.2)) 

#line plot of means with se bars
pdf(file="Pied_MortRRW_ndvi.pdf", width=8.5,height=6, paper='special')

piedsingle <- single[c(4,5,7,11,14,16), c(1, 4:18)]
piedsingle <- gather(piedsingle, year, measurement, X2000:X2014, factor_key = TRUE)
piedndvi <- ggplot(data = piedsingle, aes(x=year, y=measurement, color=ID, group=ID)) + geom_point() + geom_line() +
  xlab("Year") + ylab("NDVI") + scale_x_discrete(breaks=c("X2000", "X2005", "X2010"), labels = c("2000", "2005",
                                                                                                 "2010"))
piedzscore <- zscores[c(4,5,7,11,14,16), c(1, 4:18)]
piedzscore <- gather(piedzscore, year, measurement, X2000:X2014, factor_key = TRUE)
piedndviz <- ggplot(data = piedzscore, aes(x=year, y=measurement, color=ID, group=ID)) + geom_point() + geom_line() +
  xlab("Year") + ylab("NDVI Z-Score") + scale_x_discrete(breaks=c("X2000", "X2005", "X2010"), labels = c("2000", "2005",
                                                                                                         "2010")) +
  ggtitle("Pied") +
  scale_colour_manual(values = c("#440154FF", "#414487FF", "#2A788EFF",
                                 "#22A884FF", "#7AD151FF", "#FDE725FF")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

piedtrader <- trader[c(1,3,8,9,12,19), ]
piedtrader <- gather(piedtrader, year, trader, X2000:X2014, factor_key = TRUE)
piedtraderplot <- ggplot(data = piedtrader, aes(x=year, y=trader, fill=plot)) + geom_bar(position="dodge", stat="identity") +
  xlab("Year") + ylab("Trader detected disturbance event") + scale_x_discrete(breaks=c("X2000", "X2005", "X2010"), labels = c("2000", "2005",
                                                                                                                              "2010"))
pe1 <- ggplot(data = means, aes(x=Group.1, y=x, color=Group.2)) + geom_point() + geom_line()
pe1 <- pe1 +  geom_errorbar(aes(ymin=(means[,4]), ymax=(means[,5])), width=.2) + xlab("Year") + ylab("Observed - Predicted RRW") +
  scale_color_hue(labels =  c("AZ 435", "AZ 940", "CO 137","CO 31", "CO 95", "UT 228")) +  guides(color=guide_legend("Plot")) +
  geom_hline(yintercept=0)

plots <- data.frame (name = c("AZ 435", "UT 228", "AZ 940", "CO 137","CO 31", "CO 95"), 
                     value = c(.55, .2395, .73, .66,0, 0))

mort1 <- ggplot(data = plots, mapping = aes(x = name, y =value), pattern = deadID) +
  geom_bar(stat="identity", aes(color=name, fill=name), width=.5) + ylab("Fraction of mortality reported 2014") +
  xlab("Plot") + ylim(c(0,1)) + ggtitle("Pied plots") + theme(legend.title = element_blank())
grid.arrange(mort1, pe1, piedndviz, piedtraderplot, nrow=2,ncol=2)
dev.off()

dodge <- position_dodge(.3)
pp <- ggplot(data = means, aes(x=Group.1, y=x, color=Group.2)) + geom_point(position=dodge)
pp <- pp + geom_errorbar(aes(ymin=(means[,4]), ymax=(means[,5])), width=1, position=dodge)

#separate boxplot per plot across all years 
p1 <- ggplot(without, aes(x=year, y=diff, fill=plotID, group=year)) + geom_boxplot() +facet_wrap(~plotID)

#boxplots for each plot aggregated over years
p2 <- ggplot(without, aes(x=plotID, y=diff, fill=plotID)) + geom_boxplot()
grid.arrange(p, p2, pe1, pp, nrow=3, ncol=2)


###########################PIEN#####################
#####combine into species full 1960-1989
#Remove some CO_SF_49_PIEN
#initialize dataframe for pien - 1, 6, 9, 12, 13, 20, 21
pienclimateall <- get(paste0(corelist[1,]))
#Hard coded - which rows are potr
list <- as.data.frame(c(6, 9, 12, 13, 20, 21))
#get climate files and row bind long form
for (c in 1:nrow(list)){
  pienclimateall <- rbind(pienclimateall,get(paste0(corelist[(list[c,]),])))
}

model = lmer(value ~ year + ppt + vpd + def + tmax + aet + soil + PDSI + swe + (1|variable) + (1|plotID), 
             data = pienclimateall, REML=TRUE)
summary(model)
# 
# 
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)  1.434e+01  2.146e+00  2.670e+03   6.681 2.88e-11 ***
#   year        -6.053e-03  1.086e-03  2.930e+03  -5.573 2.73e-08 ***
#   ppt         -5.607e-03  1.698e-03  2.925e+03  -3.303 0.000969 ***
#   vpd         -1.477e+00  4.397e-01  2.926e+03  -3.358 0.000795 ***
#   def         -4.921e-03  3.322e-03  2.926e+03  -1.481 0.138669    
# tmax         3.846e-02  1.796e-02  2.916e+03   2.141 0.032356 *  
#   aet          3.825e-04  4.508e-03  2.928e+03   0.085 0.932376    
# soil        -6.516e-03  2.655e-03  2.911e+03  -2.454 0.014187 *  
#   PDSI         2.809e-02  9.567e-03  2.876e+03   2.936 0.003350 ** 
#   swe          2.883e-04  6.285e-04  2.928e+03   0.459 0.646444    

rand(model)

model1 = lmer(value ~ year  + ppt + vpd + tmax + soil + PDSI + (1|variable) + (1|plotID), 
              data = pienclimateall, REML=TRUE)
r.squaredGLMM(model1)

summary(model1)

# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)  1.498e+01  2.100e+00  2.674e+03   7.131 1.28e-12 ***
#   year        -6.402e-03  1.071e-03  2.933e+03  -5.976 2.57e-09 ***
#   ppt         -3.405e-03  1.430e-03  2.931e+03  -2.381   0.0173 *  
#   vpd         -1.853e+00  3.631e-01  2.907e+03  -5.102 3.58e-07 ***
#   tmax         4.073e-02  1.708e-02  2.883e+03   2.385   0.0172 *  
#   soil        -5.818e-03  2.611e-03  2.918e+03  -2.229   0.0259 *  
#   PDSI         3.308e-02  7.570e-03  2.901e+03   4.369 1.29e-05 ***


#initialize dataframe for pien for 2000-2013
pienclimate2013 <- get(paste0(corelist[1,],"2013"))
#Hard coded - which rows are pien
list <- as.data.frame(c(6, 9, 12, 13, 20, 21))
#get climate files and row bind long form
for (c in 1:nrow(list)){
  pienclimate2013 <- rbind(pienclimate2013,get(paste0(corelist[(list[c,]),],"2013")))
}

#get predicted values for 2005-2013
pred <- summary(model1)$coefficients[1,1] + summary(model1)$coefficients[2,1] * pienclimate2013$year + 
  summary(model1)$coefficients[3,1] * pienclimate2013$ppt +  
  summary(model1)$coefficients[4,1] * pienclimate2013$vpd +  
  summary(model1)$coefficients[5,1] * pienclimate2013$tmax +
  summary(model1)$coefficients[6,1] * pienclimate2013$soil +
  summary(model1)$coefficients[7,1] * pienclimate2013$PDSI 

pienclimate2013 <- cbind(pienclimate2013, pred)
diff <- pienclimate2013$value - pienclimate2013$pred
pienclimate2013 <- cbind(pienclimate2013, diff)

#plots looking at the data, remove outliers
boxplot(pienclimate2013$value, pienclimate2013$diff)
plot(pienclimate2013$value, pienclimate2013$pred)
plot(model1)

boxplot(diff ~ year, data=pienclimate2013)
abline(h=0)
qplot(data=pienclimate2013, x=value)
outlierReplace(pienclimate2013, "value", which(pienclimate2013$value > 2.201), NA)
withoutpien <- na.omit(pienclimate2013)

#line plots with means 
means <- aggregate(withoutpien[, 14], list(withoutpien$year, withoutpien$plotID), mean)
se <- aggregate(diff ~ year + plotID, data = withoutpien,FUN=function(x) c(mean=mean(x),se=std.error(x)))
means <- cbind(means, (se$diff[,1]-se$diff[,2]), (se$diff[,1]+se$diff[,2]))

################
####---2024 final paper edits pien SI figures (E)----####
#line plot of means with se bars!
means <- aggregate(withoutpien[, 14], list(withoutpien$year, withoutpien$plotID), mean)
se <- aggregate(diff ~ year + plotID, data = withoutpien,FUN=function(x) c(mean=mean(x),se=std.error(x)))
means <- cbind(means, (se$diff[,1]-se$diff[,2]), (se$diff[,1]+se$diff[,2]))

pienrrwse <- ggplot(data = means, aes(x=Group.1, y=x, color=Group.2)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin=(means[,4]), ymax=(means[,5])), width=.2) + xlab("Year") +
  ylab("Observed - predicted RRW")+
  scale_color_hue(labels =  c("AZ 1007","CO 108", "CO 23", "CO 49", "CO 8", "UT 57", "UT 63")) +  guides(color=guide_legend("PlotID")) +
  geom_hline(yintercept=0) + labs(tag="E")
pienrrwse
ggsave("sifigure_pien_rrw_indv",plot = pienrrwse, device = png)

plots <- data.frame (name = c("AZ 1007","CO 108", "CO 23", "CO 49", "CO 8", "UT 57", "UT 63"), value = c(.21,0,.31,.30,0,0,0))

pienmort <- ggplot(data = plots, mapping = aes(x = name, y =value)) +
  geom_bar(stat="identity", aes(color=name, fill=name)) + ylab("Fraction of mortality reported 2014") +
  theme(legend.title = element_blank()) + labs(tag="E") + xlab("Plot") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  scale_y_continuous(limits = c(0, 1))

pienmort
ggsave("sifigure_pien_mort",plot = pienrrwse, device = png)


#Single pixel 
piensingle <- single[c(1,6,9,12,13,18,19), c(1, 4:18)]
piensingle <- gather(piensingle, year, measurement, X2000:X2014, factor_key = TRUE)
pienndvi <- ggplot(data = piensingle, aes(x=year, y=measurement, color=ID, group=ID)) + geom_point() + geom_line() +
  xlab("Year") + ylab("NDVI") + scale_x_discrete(breaks=c("X2000", "X2005", "X2010"), labels = c("2000", "2005",
                                                                                                 "2010"))

pienzscore <- zscores[c(1,6,9,12,13,18,19), c(1, 4:18)]
pienzscore <- gather(pienzscore, year, measurement, X2000:X2014, factor_key = TRUE)
pienndviz <- ggplot(data = pienzscore, aes(x=year, y=measurement, color=ID, group=ID)) + geom_point() + geom_line() +
  xlab("Year") + ylab("NDVI Z-score") + scale_x_discrete(breaks=c("X2000", "X2005", "X2010"), labels = c("2000", "2005",
                                                                                                         "2010")) +
  ggtitle("Pien")+
  scale_colour_manual(values = c("#440154FF",
                                 "#443A83FF", "#31688EFF", "#21908CFF",
                                 "#35B779FF", "#8FD744FF", "#FDE725FF")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

pientrader <- trader[c(6,15,16,17,18,23,24), ]
pientrader <- gather(pientrader, year, trader, X2000:X2014, factor_key = TRUE)
pientraderplot <- ggplot(data = pientrader, aes(x=year, y=trader, fill=plot)) + geom_bar(position="dodge", stat="identity") +
  xlab("Year") + ylab("Trader detected disturbance event") + scale_x_discrete(breaks=c("X2000", "X2005", "X2010"), labels = c("2000", "2005",
                                                                                                                              "2010"))
grid.arrange(mort, pe, pienndviz, pientraderplot, nrow=2,ncol=2)
dev.off()



###########################pipo#####################
#####combine into species full 1960-1989
#initialize dataframe for pipo - 3, 8, 10, 15:17, 19, 22, 
#Hard coded - which rows are pipo
pipoclimateall <- get(paste0(corelist[3,])) 
list <- as.data.frame(c(8, 10, 15:17, 19, 22))

#get climate files and row bind long form
for (c in 1:nrow(list)){
  pipoclimateall <- rbind(pipoclimateall,get(paste0(corelist[(list[c,]),])))
}

model = lmer(value ~ year + ppt + vpd + def + tmax + aet + soil + PDSI + swe + (1|variable) + (1|plotID), 
             data = pipoclimateall, REML=TRUE)
summary(model)
# 
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept) -1.226e+01  2.135e+00  3.099e+03  -5.742 1.03e-08 ***
#   year         6.387e-03  1.097e-03  2.708e+03   5.821 6.53e-09 ***
#   ppt         -1.261e-03  1.903e-03  2.943e+03  -0.663  0.50758    
# vpd         -1.552e-01  4.203e-01  1.360e+03  -0.369  0.71196    
# def          1.397e-03  3.330e-03  2.996e+03   0.420  0.67485    
# tmax         3.222e-02  1.952e-02  8.024e+02   1.651  0.09920 .  
# aet          1.395e-02  4.246e-03  3.252e+03   3.286  0.00103 ** 
#   soil         2.931e-03  1.629e-03  1.771e+03   1.800  0.07210 .  
# PDSI         7.385e-02  9.346e-03  1.971e+03   7.902 4.53e-15 ***
#   swe         -4.782e-03  1.185e-03  3.170e+03  -4.034 5.60e-05 ***

model1 = lmer(value ~ year + aet + PDSI + swe + (1|variable) + (1|plotID), 
              data = pipoclimateall, REML=TRUE)
r.squaredGLMM(model1)

summary(model1)
# # 
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept) -1.358e+01  2.014e+00  3.279e+03  -6.744 1.81e-11 ***
#   year         7.306e-03  1.012e-03  3.276e+03   7.221 6.38e-13 ***
#   aet          1.170e-02  2.173e-03  3.142e+03   5.384 7.84e-08 ***
#   PDSI         7.345e-02  7.909e-03  3.246e+03   9.286  < 2e-16 ***
#   swe         -5.022e-03  1.149e-03  3.224e+03  -4.370 1.28e-05 ***

#initialize dataframe for pipo for 2005-2013
pipoclimate2013 <- get(paste0(corelist[3,], "2013"))
#Hard coded - which rows are pipo
list <- as.data.frame(c(8, 10, 15:17, 19, 22))

#get climate files and row bind long form
for (c in 1:nrow(list)){
  pipoclimate2013 <- rbind(pipoclimate2013,get(paste0(corelist[(list[c,]),],"2013")))
}

#get predicted values for 2005-2013
pred <- summary(model1)$coefficients[1,1] + summary(model1)$coefficients[2,1] * pipoclimate2013$year + 
  summary(model1)$coefficients[3,1] * pipoclimate2013$aet +  
  summary(model1)$coefficients[4,1] * pipoclimate2013$PDSI +
  summary(model1)$coefficients[5,1] * pipoclimate2013$swe 

pipoclimate2013 <- cbind(pipoclimate2013, pred)
diff <- pipoclimate2013$value - pipoclimate2013$pred
pipoclimate2013 <- cbind(pipoclimate2013, diff)

#plots looking at the data, remove outliers
boxplot(pipoclimate2013$value, pipoclimate2013$diff)
plot(pipoclimate2013$value, pipoclimate2013$pred)
plot(model1)
min(boxplot(pipoclimate2013$value)$out) #1.876
qplot(data=pipoclimate2013, x=value)
min(boxplot(pipoclimate2013$value)$out) #1.824

outlierReplace(pipoclimate2013, "value", which(pipoclimate2013$value > 1.824), NA)
qplot(data=pipoclimate2013, x=value)

#line plots with means 
means <- aggregate(withoutpipo[, 14], list(withoutpipo$year, withoutpipo$plotID), mean)
se <- aggregate(diff ~ year + plotID, data = withoutpipo,FUN=function(x) c(mean=mean(x),se=std.error(x)))
means <- cbind(means, (se$diff[,1]-se$diff[,2]), (se$diff[,1]+se$diff[,2]))


###########################
####---2024 final paper edits pipo SI figures (D)---####
means <- aggregate(withoutpipo[, 14], list(withoutpipo$year, withoutpipo$plotID), mean)
se <- aggregate(diff ~ year + plotID, data = withoutpipo,FUN=function(x) c(mean=mean(x),se=std.error(x)))
means <- cbind(means, (se$diff[,1]-se$diff[,2]), (se$diff[,1]+se$diff[,2]))

piporrwse <- ggplot(data = means, aes(x=Group.1, y=x, color=Group.2)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin=(means[,4]), ymax=(means[,5])), width=.2) + xlab("Year") +
  ylab("Observed - predicted RRW")+
  scale_color_hue(labels =  c("AZ 11","AZ 23", "CO 165", "CO 3010", "NM 350", "NM 425", "UT 183", "UT 3009", "UT 82")) +  
  guides(color=guide_legend("PlotID")) +
  geom_hline(yintercept=0) +
  xlab("Year") + ylab("Observed - predicted RRW") + labs(tag="D")
piporrwse
ggsave("sifigure_pipo_rrw_indv",plot = piporrwse, device = png)


#Mortality plot
plots <- data.frame (name = c("AZ 11","AZ 23", "CO 165", "CO 3010", "UT 183", "UT 3009", "UT 82"), 
                     value = c(.17,0,.16,0,0,0,.01))

pipomort <- ggplot(data = plots, mapping = aes(x = name, y =value)) +
  geom_bar(stat="identity", aes(color=name, fill=name)) + ylab("Fraction of mortality reported 2014")+ 
  theme(legend.title = element_blank()) + xlab("Plot") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position="none") + labs(tag ="D") +
  scale_y_continuous(limits = c(0, 1))

pipomort
ggsave("sifigure_pipo_mort",plot = pipomort, device = png)

#Single pixel 
piposingle <- single[c(2,3,8,10,17,19,22), c(1, 4:18)]
piposingle <- gather(piposingle, year, measurement, X2000:X2014, factor_key = TRUE)
pipondvi <- ggplot(data = piposingle, aes(x=year, y=measurement, color=ID, group=ID)) + geom_point() + geom_line() +
  xlab("Year") + ylab("NDVI") + scale_x_discrete(breaks=c("X2000", "X2005", "X2010"), labels = c("2000", "2005",
                                                                                                 "2010"))

pipozscore <- zscores[c(2,3,8,10,17,15,20), c(1, 4:18)]
pipozscore <- gather(pipozscore, year, measurement, X2000:X2014, factor_key = TRUE)
pipondviz <- ggplot(data = pipozscore, aes(x=year, y=measurement, color=ID, group=ID)) + geom_point() + geom_line() +
  xlab("Year") + ylab("NDVI Z-score") + scale_x_discrete(breaks=c("X2000", "X2005", "X2010"), labels = c("2000", "2005",
                                                                                                         "2010")) +
  ggtitle("Pipo") + scale_colour_manual(values = c("#440154FF",
                                                   "#443A83FF", "#31688EFF", "#21908CFF",
                                                   "#35B779FF", "#8FD744FF", "#FDE725FF")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

pipotrader <- trader[c(4,5,13,14,20,21,22), ]
pipotrader <- gather(pipotrader, year, trader, X2000:X2014, factor_key = TRUE)
pipotraderplot <- ggplot(data = pipotrader, aes(x=year, y=trader, fill=plot)) + geom_bar(position="dodge", stat="identity") +
  xlab("Year") + ylab("Trader detected disturbance event") + 
  scale_x_discrete(breaks=c("X2000", "X2005", "X2010"), labels = c("2000", "2005", "2010"))
grid.arrange(mort, pe, pipondviz, pipotraderplot, nrow=2,ncol=2)
dev.off()


##############
####----2024 combined SI figures----####
#RRW plot
ggarrange(potrrrwse, juosrrwse, peidrrwse, piporrwse, pienrrwse)
mortcombo <- ggarrange(potrmort, juosmort, peidmort, pipomort, pienmort)
mortcombo
ggsave("s7_sifigure_all_mort.png",plot = mortcombo, device = png, width = 7, height = 8, units = "in")


######Figure 2 - combined plot with all five species renamed#######
juosndviz <- juosndviz + ggtitle("Juniper")
piedndviz <- piedndviz + ggtitle("Pinyon")
pienndviz <- pienndviz + ggtitle("Spruce")
potrndviz <- potrndviz + ggtitle("Aspen")
pipondviz <- pipondviz + ggtitle("Ponderosa")

pdf(file="NDVIallTGIF.pdf", width=13,height=7, paper='special')
grid.arrange(potrndviz, juosndviz, piedndviz, pipondviz, pienndviz, nrow=2,ncol=3)

dev.off()


#boxplots of residuals 1-3 years, 1-5 years, 1-10 years after drought (2002)
three <- c("2003", "2004", "2005")
five <- c("2003", "2004", "2005", "2006", "2007")
ten <- c("2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013")

#Potr
withoutpotr <- withoutpotr
threeyearspotr <- withoutpotr[which(withoutpotr$year==three),]
threeyearspotr$group <- "3 year"
fiveyearspotr <- withoutpotr[which(withoutpotr$year==five),]
fiveyearspotr$group <- "5 year"
tenyearspotr <- withoutpotr[which(withoutpotr$year==ten),]
tenyearspotr$group <- "10 year"
#group 3, 5, and 10 year subsets together
fullpotryears <- rbind(threeyearspotr, fiveyearspotr, tenyearspotr)
fullpotryears$group <- factor(fullpotryears$group , levels=c("3 year", "5 year", "10 year"))

#Boxplot with points
potrboxjitter <- ggplot( data = fullpotryears, aes(x=group, y=diff, group=group, color = group)) +
  geom_boxplot() +
  geom_jitter(size=0.4, alpha=0.9
  ) +
  ggtitle("Potr") +
  xlab("Number of years post drought") +
  ylab("Observed - Predicted RRW") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position="none",
                     plot.title = element_text(size=11)) +
  scale_color_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"))


#T-test: H0 true mean is equal to 0, HA true mean is not equal to 0 
t.test(threeyearspotr$diff, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, 
       paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(fiveyearspotr$diff, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, 
       paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(tenyearspotr$diff, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, 
       paired = FALSE, var.equal = FALSE, conf.level = 0.95)

##Anova
fullpotrbins <- rbind(threeyearspotr, fiveyearspotr, tenyearspotr)
potraov <- aov(diff ~ group, data = fullpotrbins)
summary(potraov) #p = 0.0904
#TukeyHSD(potraov) 

##Add p-values to graph 
prot2 = distinct(fullpotryears, group) %>%
  arrange(group)
prot2$yloc = 1.5
prot2$label = c("p < 0.001", "p < 0.001", "p < 0.001")
letter <- data.frame(x=.5,y = -1.5, label="A")

potrboxjitter <- potrboxjitter + 
  ylim(-1.5,1.5) +
  geom_text(data = prot2, aes(y = yloc, label = label), 
            position = position_dodge(width = .75)) + geom_hline(yintercept=0)+ 
  geom_text(data = letter, aes(x = x, y = y, label = label), inherit.aes = FALSE)

#######juos
withoutjuos <- withoutjuos
threeyearsjuos <- withoutjuos[which(withoutjuos$year==three),]
threeyearsjuos$group <- "3 year"
fiveyearsjuos <- withoutjuos[which(withoutjuos$year==five),]
fiveyearsjuos$group <- "5 year"
tenyearsjuos <- withoutjuos[which(withoutjuos$year==ten),]
tenyearsjuos$group <- "10 year"

#group 3, 5, and 10 year subsets together
fulljuosyears <- rbind(threeyearsjuos, fiveyearsjuos, tenyearsjuos)
fulljuosyears$group <- factor(fulljuosyears$group , levels=c("3 year", "5 year", "10 year"))

#Boxplot with points
juosboxjitter <- ggplot( data = fulljuosyears, aes(x=group, y=diff, group=group, color = group)) +
  geom_boxplot() +
  geom_jitter(size=0.4, alpha=0.9) +
  ggtitle("Juos") +
  xlab("Number of years post drought") +
  ylab("Observed - Predicted RRW")+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position="none",
                     plot.title = element_text(size=11)) +
  scale_color_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"))


##T-tests
t.test(threeyearsjuos$diff, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, 
       paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(fiveyearsjuos$diff, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, 
       paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(tenyearsjuos$diff, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, 
       paired = FALSE, var.equal = FALSE, conf.level = 0.95)

##Anova
fulljuosbins <- rbind(threeyearsjuos, fiveyearsjuos, tenyearsjuos)
juosaov <- aov(diff ~ group, data = fulljuosbins)
summary(juosaov) #p = 0.085
#TukeyHSD(juosaov) 

##Add p-values to graph 
prot2 = distinct(fulljuosyears, group) %>%
  arrange(group)
prot2$yloc = 1.5
prot2$label = c("p = 0.005", "p = 0.003", "p = 0.96")
letter <- data.frame(x=.5,y = -1.5, label="B")

juosboxjitter <- juosboxjitter + 
  ylim(-1.5,1.5) +
  geom_text(data = prot2, aes(y = yloc, label = label), 
            position = position_dodge(width = .75)) + geom_hline(yintercept=0)+ 
  geom_text(data = letter, aes(x = x, y = y, label = label), inherit.aes = FALSE)



#####pipo
withoutpipo <- withoutpipo
threeyearspipo <- withoutpipo[which(withoutpipo$year==three),]
threeyearspipo$group <- "3 year"
fiveyearspipo <- withoutpipo[which(withoutpipo$year==five),]
fiveyearspipo$group <- "5 year"
tenyearspipo <- withoutpipo[which(withoutpipo$year==ten),]
tenyearspipo$group <- "10 year"

#group 3, 5, and 10 year subsets together
fullpipoyears <- rbind(threeyearspipo, fiveyearspipo, tenyearspipo)
fullpipoyears$group <- factor(fullpipoyears$group , levels=c("3 year", "5 year", "10 year"))

#Boxplot with points
pipoboxjitter <- ggplot( data = fullpipoyears, aes(x=group, y=diff, group=group, color = group)) +
  geom_boxplot() +
  geom_jitter( size=0.4, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Pipo") +
  xlab("Number of years post drought") +
  ylab("Observed - Predicted RRW")+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position="none",
                     plot.title = element_text(size=11)) +
  scale_color_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"))

##T-tests
t.test(threeyearspipo$diff, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, 
       paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(fiveyearspipo$diff, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, 
       paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(tenyearspipo$diff, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, 
       paired = FALSE, var.equal = FALSE, conf.level = 0.95)


##Anova
fullpipobins <- rbind(threeyearspipo, fiveyearspipo, tenyearspipo)
pipoaov <- aov(diff ~ group, data = fullpipobins)
summary(pipoaov) #p = 0.00126 **
TukeyHSD(pipoaov) 
# diff        lwr         upr     p adj
# 3 year-10 year -0.14728115 -0.2637126 -0.03084973 0.0087190
# 5 year-10 year -0.16582698 -0.2797439 -0.05191009 0.0019783
# 5 year-3 year  -0.01854582 -0.1302579  0.09316621 0.9192874

##Add p-values to graph 
prot2 = distinct(fullpipoyears, group) %>%
  arrange(group)
prot2$yloc = 1.5
prot2$label = c("p < 0.001", "p < 0.001", "p < 0.001")
letter <- data.frame(x=.5,y = -1.5, label="D")

pipoboxjitter <- pipoboxjitter + 
  ylim(-1.5,1.5) +
  geom_text(data = prot2, aes(y = yloc, label = label), 
            position = position_dodge(width = .75)) + geom_hline(yintercept=0)+ 
  geom_text(data = letter, aes(x = x, y = y, label = label), inherit.aes = FALSE)




#pien
withoutpien <- withoutpien
threeyearspien <- withoutpien[which(withoutpien$year==three),]
threeyearspien$group <- "3 year"
fiveyearspien <- withoutpien[which(withoutpien$year==five),]
fiveyearspien$group <- "5 year"
tenyearspien <- withoutpien[which(withoutpien$year==ten),]
tenyearspien$group <- "10 year"

#group 3, 5, and 10 year subsets together
fullpienyears <- rbind(threeyearspien, fiveyearspien, tenyearspien)
fullpienyears$group <- factor(fullpienyears$group , levels=c("3 year", "5 year", "10 year"))

#Boxplot with points
pienboxjitter <- ggplot( data = fullpienyears, aes(x=group, y=diff, group=group, color = group)) +
  geom_boxplot() +
  geom_jitter( size=0.4, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Pien") +
  xlab("Number of years post drought") +
  ylab("Observed - Predicted RRW") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position="none",
                     plot.title = element_text(size=11)) +
  scale_color_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"))


##T-tests
t.test(threeyearspien$diff, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, 
       paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(fiveyearspien$diff, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, 
       paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(tenyearspien$diff, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, 
       paired = FALSE, var.equal = FALSE, conf.level = 0.95)

##Anova
fullpienbins <- rbind(threeyearspien, fiveyearspien, tenyearspien)
pienaov <- aov(diff ~ group, data = fullpienbins)
summary(pienaov) #p = 0.000801
TukeyHSD(pienaov) 
# $group
# diff        lwr         upr     p adj
# 3 year-10 year -0.2289376 -0.3855941 -0.07228106 0.0019072
# 5 year-10 year -0.2186741 -0.3756752 -0.06167298 0.0033187
# 5 year-3 year   0.0102635 -0.1385849  0.15911191 0.9855571

##Add p-values to graph 
prot2 = distinct(fullpienyears, group) %>%
  arrange(group)
prot2$yloc = 1.5
prot2$label = c("p < 0.001", "p < 0.001", "p < 0.001")
letter <- data.frame(x=.5,y = -1.5, label="E")

pienboxjitter <- pienboxjitter + 
  ylim(-1.5,1.5) +
  geom_text(data = prot2, aes(y = yloc, label = label), 
            position = position_dodge(width = .75)) + geom_hline(yintercept=0)+ 
  geom_text(data = letter, aes(x = x, y = y, label = label), inherit.aes = FALSE)


###pied
withoutpied <- withoutpied
threeyearspied <- withoutpied[which(withoutpied$year==three),]
threeyearspied$group <- "3 year"
fiveyearspied <- withoutpied[which(withoutpied$year==five),]
fiveyearspied$group <- "5 year"
tenyearspied <- withoutpied[which(withoutpied$year==ten),]
tenyearspied$group <- "10 year"

#group 3, 5, and 10 year subsets together
fullpiedyears <- rbind(threeyearspied, fiveyearspied, tenyearspied)
fullpiedyears$group <- factor(fullpiedyears$group , levels=c("3 year", "5 year", "10 year"))

#Boxplot with points
piedboxjitter <- ggplot( data = fullpiedyears, aes(x=group, y=diff, group=group, color = group)) +
  geom_boxplot() +
  geom_jitter(size=0.4, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Pied") +
  xlab("Number of years post drought") +
  ylab("Observed - Predicted RRW") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position="none",
                     plot.title = element_text(size=11)) +
  scale_color_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"))


##T-tests
t.test(threeyearspied$diff, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, 
       paired = FALSE, var.equal = FALSE, conf.level = 0.95) #.651
t.test(fiveyearspied$diff, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, 
       paired = FALSE, var.equal = FALSE, conf.level = 0.95) #.8897
t.test(tenyearspied$diff, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, 
       paired = FALSE, var.equal = FALSE, conf.level = 0.95) #.3598

##Anova
fullpiedbins <- rbind(threeyearspied, fiveyearspied, tenyearspied)
piedaov <- aov(diff ~ group, data = fullpiedbins)
summary(piedaov) #p = 0.875


##Add p-values to graph 
prot2 = distinct(fullpiedyears, group) %>%
  arrange(group)
prot2$yloc = 1.5
prot2$label = c("p = 0.65", "p = 0.89", "p = 0.36")
letter <- data.frame(x=.5,y = -1.5, label="C")
piedboxjitter <- piedboxjitter + 
  ylim(-1.5,1.5) +
  geom_text(data = prot2, aes(y = yloc, label = label), 
            position = position_dodge(width = .75)) + geom_hline(yintercept=0) + 
  geom_text(data = letter, aes(x = x, y = y, label = label), inherit.aes = FALSE)


######Pdf exportation of figure 3#####  ######Pdf exportation of figure 3#####  
juosboxjitter <- juosboxjitter + ggtitle("Juniper")
piedboxjitter <- piedboxjitter + ggtitle("Pinyon")
pienboxjitter <- pienboxjitter + ggtitle("Spruce")
potrboxjitter <- potrboxjitter + ggtitle("Aspen")
pipoboxjitter <- pipoboxjitter + ggtitle("Ponderosa")

pdf(file="BoxplotJitterFull.pdf", width=9,height=7, paper='special')
grid.arrange(potrboxjitter, (juosboxjitter + geom_hline(yintercept=0)), piedboxjitter, pipoboxjitter, pienboxjitter, nrow=2,ncol=3)
dev.off()

#####Combined means of each plot for the five year bins#####
combinedfiveyears <- rbind(fiveyearsjuos,fiveyearspied,fiveyearspien,fiveyearspipo,fiveyearspotr)
fullmeandiff <- as.data.frame(sapply(split(combinedfiveyears$diff, combinedfiveyears$plotID),mean))
colnames(fullmeandiff) <- c("diff")
fullmeandiff <- cbind(fullmeandiff, unique(combinedfiveyears$plotID))
colnames(fullmeandiff) <- c("diff", "plotID")

#Extract climate variables for each plot by removing individual tree cores (same climate data)
piedclimateplot <- piedclimateall[!duplicated(piedclimateall[,c(1,12)]),]
juosclimateplot <- juosclimateall[!duplicated(juosclimateall[,c(1,12)]),]
pienclimateplot <- pienclimateall[!duplicated(pienclimateall[,c(1,12)]),]
potrclimateplot <- POTRclimateall[!duplicated(POTRclimateall[,c(1,12)]),]
pipoclimateplot <- pipoclimateall[!duplicated(pipoclimateall[,c(1,12)]),]

#combine into one dataframe 
combinedclimateplot <- rbind(juosclimateplot, piedclimateplot, pienclimateplot, potrclimateplot, pipoclimateplot)
combinedclimateplot <- na.omit(combinedclimateplot)

#####Subset climate variable - PDSI 
pdsiclimateplot <- as.data.frame(sapply(split(combinedclimateplot$PDSI, combinedclimateplot$plotID),mean))
pdsiclimateplot <- cbind(pdsiclimateplot, unique(combinedclimateplot$plotID))
names(pdsiclimateplot)[1] <- "pdsi"
names(pdsiclimateplot)[2] <- "plotID"

#merge RRW difference and psdi climate mean for 1960-1989
mergedpdsi  <- left_join(fullmeandiff, pdsiclimateplot, by = 'plotID')

#Hard coding adding species column - edited to remove AZ11 Pipo
mergedpdsi[1:4,4] <- 'juos'
mergedpdsi[5:9,4] <- 'pied'
mergedpdsi[10:16,4] <- 'pien'
mergedpdsi[17:24,4] <- 'pipo'
mergedpdsi[25:36,4] <- 'potr'
names(mergedpdsi)[4] <- "species"

#mixed effects model (need to separate by species, not plot)
pdsilm <- lm(mergedpdsi$diff ~ mergedpdsi$pdsi + (1/mergedpdsi$species))
summary(pdsilm)
# Coefficients:
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)     -0.10587    0.08041  -1.317    0.197
# mergedpdsi$pdsi -0.73464    0.98556  -0.745    0.461 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4343 on 36 degrees of freedom
# Multiple R-squared:  0.001745,	Adjusted R-squared:  -0.02598 
# F-statistic: 0.06293 on 1 and 36 DF,  p-value: 0.8034
AIC(pdsilm) #53.62

pdsilm1 <- lm(mergedpdsi$diff ~ mergedpdsi$pdsi + (mergedpdsi$pdsi/mergedpdsi$species))
summary(pdsilm1)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                            -0.10655    0.08885  -1.199    0.240
# mergedpdsi$pdsi                        -2.78923    5.42398  -0.514    0.611
# mergedpdsi$pdsi:mergedpdsi$speciespied  1.76884    6.12235   0.289    0.775
# mergedpdsi$pdsi:mergedpdsi$speciespien  2.95120    5.84998   0.504    0.618
# mergedpdsi$pdsi:mergedpdsi$speciespipo  1.57229    5.71330   0.275    0.785
# mergedpdsi$pdsi:mergedpdsi$speciespotr  2.71862    5.89235   0.461    0.648
# 
# Residual standard error: 0.5095 on 30 degrees of freedom
# Multiple R-squared:  0.03178,	Adjusted R-squared:  -0.1296 
# F-statistic: 0.197 on 5 and 30 DF,  p-value: 0.9612
AIC(pdsilm1) #61.19


###6.4.2021 Calculating RSS from model
pdsilm2 <- lm(mergedpdsi$diff ~ mergedpdsi$pdsi + (mergedpdsi$species))
drop1(pdsilm2, test = "F")
# 
# Model:
#   mergedpdsi$diff ~ mergedpdsi$pdsi + (mergedpdsi$species)
#                   Df Sum of Sq    RSS      AIC F value    Pr(>F)    
# <none>                          1.3391 -106.494                      
# mergedpdsi$pdsi     1    0.0166 1.3558 -108.049  0.3725    0.5462    
# mergedpdsi$species  4    6.5735 7.9126  -50.542 36.8155 3.603e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Linear regression with 95% confidence intervals
pdsici <- ggplot(data = mergedpdsi, aes(x=pdsi, y=diff, group=species, color=species)) + geom_point() + 
  geom_smooth(method = "lm") + xlab("PDSI") + ylab("Observed - Predicted RRW")
##Add R2 values to graph 
pdsici <- pdsici + 
  annotate(geom="text",x=-.05, y=1, label="Multiple R2 = .03")


#############Anova###########
#One way ANOVA for difference in five year post drought RRW across species
#mergedpdsi <- mergedpdsi[1:38,]
mergedpdsi$species <- as.factor(mergedpdsi$species)
mergedpdsi$diff <- as.numeric(as.character(mergedpdsi$diff))
levels(mergedpdsi$species)
pdsiaov <- aov(mergedpdsi$diff ~ mergedpdsi$species)
summary(pdsiaov)                   
# Df Sum Sq Mean Sq F value  Pr(>F)    
# mergedpdsi$species  4  5.563   1.391    17.4 1.4e-07 ***
#   Residuals          31  2.478   0.080   

TukeyHSD(pdsiaov)
# $`mergedpdsi$species`
# $`mergedpdsi$species`
# diff         lwr         upr     p adj
# pied-juos -0.08155052 -0.63063550  0.46753446 0.9925161
# pien-juos -0.64528757 -1.15832687 -0.13224828 0.0080834*
# pipo-juos -0.58534253 -1.07721638 -0.09346868 0.0133422*
# potr-juos  0.28707836 -0.19083877  0.76499549 0.4261124
# pien-pied -0.56373705 -1.04301771 -0.08445639 0.0147522*
# pipo-pied -0.50379201 -0.96034487 -0.04723915 0.0247765*
# potr-pied  0.36862889 -0.07285205  0.81010982 0.1375329
# pipo-pien  0.05994504 -0.35255407  0.47244416 0.9931125
# potr-pien  0.93236594  0.53661295  1.32811892 0.0000012*
# potr-pipo  0.87242089  0.50452012  1.24032167 0.0000010*


#####Subset climate variable - ppt 
pptclimateplot <- as.data.frame(sapply(split(combinedclimateplot$ppt, combinedclimateplot$plotID),mean))
pptclimateplot <- cbind(pptclimateplot, unique(combinedclimateplot$plotID))
names(pptclimateplot)[1] <- "ppt"
names(pptclimateplot)[2] <- "plotID"

#merge RRW difference and psdi climate mean for 1960-1989
mergedppt  <- left_join(fullmeandiff, pptclimateplot, by = 'plotID')

#Hard coding adding species column
mergedppt[1:4,4] <- 'juos'
mergedppt[5:9,4] <- 'pied'
mergedppt[10:16,4] <- 'pien'
mergedppt[17:24,4] <- 'pipo'
mergedppt[25:36,4] <- 'potr'
names(mergedppt)[4] <- "species"

#mixed effects model (need to separate by species, not plot)
pptlm <- lm(mergedppt$diff ~ mergedppt$ppt + (1/mergedppt$species))
summary(pptlm)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)   -0.743950   0.403313  -1.845   0.0738 .
# mergedppt$ppt  0.015716   0.009735   1.614   0.1157  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4687 on 34 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.0712,	Adjusted R-squared:  0.04388 
# F-statistic: 2.606 on 1 and 34 DF,  p-value: 0.1157
AIC(pptlm) #51.54

pptlm1 <- lm(mergedppt$diff ~ mergedppt$ppt + (mergedppt$ppt/mergedppt$species))
summary(pptlm1)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         -0.165435   0.225747  -0.733 0.469346    
# mergedppt$ppt                        0.006852   0.006832   1.003 0.323892    
# mergedppt$ppt:mergedppt$speciespied -0.002179   0.003971  -0.549 0.587325    
# mergedppt$ppt:mergedppt$speciespien -0.015523   0.003717  -4.176 0.000235 ***
#   mergedppt$ppt:mergedppt$speciespipo -0.019601   0.003734  -5.249 1.15e-05 ***
#   mergedppt$ppt:mergedppt$speciespotr  0.005166   0.003556   1.453 0.156604    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2197 on 30 degrees of freedom
# Multiple R-squared:  0.8199,	Adjusted R-squared:  0.7899 
# F-statistic: 27.31 on 5 and 30 DF,  p-value: 2.542e-10

AIC(pptlm1) #49

pptci <- ggplot(data = mergedppt, aes(x=ppt, y=diff, group=species, color=species)) + geom_point() + geom_smooth(method = "lm") + 
  geom_smooth(method = "lm") + xlab("Precipitation") + ylab("Observed - Predicted RRW")
##Add R2 values to graph 
pptci <- pptci + 
  annotate(geom="text",x=10, y=1, label="Multiple R2 = .82")


#####Subset climate variable - def 
defclimateplot <- as.data.frame(sapply(split(combinedclimateplot$def, combinedclimateplot$plotID),mean))
defclimateplot <- cbind(defclimateplot, unique(combinedclimateplot$plotID))
names(defclimateplot)[1] <- "def"
names(defclimateplot)[2] <- "plotID"

#merge RRW difference and psdi climate mean for 1960-1989
mergeddef  <- left_join(fullmeandiff, defclimateplot, by = 'plotID')

#Hard coding adding species column
mergeddef[1:4,4] <- 'juos'
mergeddef[5:9,4] <- 'pied'
mergeddef[10:16,4] <- 'pien'
mergeddef[17:24,4] <- 'pipo'
mergeddef[25:36,4] <- 'potr'
names(mergeddef)[4] <- "species"

#mixed effects model (need to separate by species, not plot)
deflm <- lm(mergeddef$diff ~ mergeddef$def + (1/mergeddef$species))
summary(deflm)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)   -0.016718   0.210674  -0.079    0.937
# mergeddef$def -0.001990   0.004376  -0.455    0.652
# 
# Residual standard error: 0.4849 on 34 degrees of freedom
# Multiple R-squared:  0.006044,	Adjusted R-squared:  -0.02319 
# F-statistic: 0.2067 on 1 and 34 DF,  p-value: 0.6522
AIC(deflm) #54

deflm1 <- lm(mergeddef$diff ~ mergeddef$def + (mergeddef$def/mergeddef$species))
summary(deflm1)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         -0.135515   0.155202  -0.873  0.38952    
# mergeddef$def                        0.003321   0.003224   1.030  0.31131    
# mergeddef$def:mergeddef$speciespied -0.001088   0.002685  -0.405  0.68826    
# mergeddef$def:mergeddef$speciespien -0.014614   0.004337  -3.369  0.00209 ** 
#   mergeddef$def:mergeddef$speciespipo -0.012247   0.002692  -4.549  8.3e-05 ***
#   mergeddef$def:mergeddef$speciespotr  0.011217   0.003460   3.242  0.00291 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.255 on 30 degrees of freedom
# Multiple R-squared:  0.7574,	Adjusted R-squared:  0.717 
# F-statistic: 18.73 on 5 and 30 DF,  p-value: 1.986e-08
AIC(deflm1) #11


##########Bill stats#######
r.squaredGLMM(deflm) #0.0059
r.squaredGLMM(deflm1) #0.73

########Manual stats 
deflm2 <- lm(mergeddef$diff ~ mergeddef$def + (mergeddef$species))
drop1(deflm2, test='F')

#find sse
sse <- sum((fitted(deflm2) - mergeddef$diff)^2)
sse #1.33

#find ssr
ssr <- sum((fitted(deflm2) - mean(mergeddef$diff))^2)
ssr #6.71

#find sst
sst <- ssr + sse
sst #8.04
ssr/sst #0.834 R2

########Manual stats 
deflm3 <- lm(mergeddef$diff ~ mergeddef$def + (1/mergeddef$species))
drop1(deflm3, test='F')

#find sse
sse <- sum((fitted(deflm3) - mergeddef$diff)^2)
sse #8.00

#find ssr
ssr <- sum((fitted(deflm3) - mean(mergeddef$diff))^2)
ssr #0.049

#find sst
sst <- ssr + sse
sst #8.04
ssr/sst #0.006 R2

#########If we want species names
mergeddef$species <- gsub('juos', 'Juniper', mergeddef$species)
mergeddef$species <- gsub('pien', 'Spruce', mergeddef$species)
mergeddef$species <- gsub('pipo', 'Ponderosa', mergeddef$species)
mergeddef$species <- gsub('potr', 'Aspen', mergeddef$species)
mergeddef$species <- gsub('pied', 'Pinyon', mergeddef$species)

#Rename column with capital letter for legend 
names(mergeddef)[4] <- "Species"
#get R2 superscript
R2.expdef <- expression(paste("Marginal ",R^2 ,"= 0.73"))

defci <- ggplot(data = mergeddef, aes(x=def, y=diff, group=Species, color=Species, linetype=Species)) + geom_point() +  
  geom_smooth(method = "lm") + 
  xlab("Climatic water deficit") + ylab("Observed - Predicted RRW") +
  scale_linetype_manual(values=c("dashed", "dashed", "solid", "solid", "solid")) +
  annotate(geom="text",x=20, y=1, label=R2.expdef) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#A6761D")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


#####Subset climate variable - swe 
sweclimateplot <- as.data.frame(sapply(split(combinedclimateplot$swe, combinedclimateplot$plotID),mean))
sweclimateplot <- cbind(sweclimateplot, unique(combinedclimateplot$plotID))
names(sweclimateplot)[1] <- "swe"
names(sweclimateplot)[2] <- "plotID"

#merge RRW difference and psdi climate mean for 1960-1989
mergedswe  <- left_join(fullmeandiff, sweclimateplot, by = 'plotID')

#Hard coding adding species column
mergedswe[1:4,4] <- 'juos'
mergedswe[5:9,4] <- 'pied'
mergedswe[10:16,4] <- 'pien'
mergedswe[17:24,4] <- 'pipo'
mergedswe[25:36,4] <- 'potr'
names(mergedswe)[4] <- "species"

#mixed effects model (need to separate by species, not plot)
swelm <- lm(mergedswe$diff ~ mergedswe$swe + (1/mergedswe$species))
summary(swelm)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)   -0.269998   0.147677  -1.828   0.0763 .
# mergedswe$swe  0.003478   0.002632   1.321   0.1952  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4743 on 34 degrees of freedom
# Multiple R-squared:  0.04884,	Adjusted R-squared:  0.02087 
# F-statistic: 1.746 on 1 and 34 DF,  p-value: 0.1952
AIC(swelm) #52

swelm1 <- lm(mergedswe$diff ~ mergedswe$swe + (mergedswe$swe/mergedswe$species))
summary(swelm1)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         -0.0465028  0.1185363  -0.392  0.69760    
# mergedswe$swe                        0.0046284  0.0061281   0.755  0.45598    
# mergedswe$swe:mergedswe$speciespied -0.0042820  0.0068601  -0.624  0.53722    
# mergedswe$swe:mergedswe$speciespien -0.0116412  0.0055101  -2.113  0.04306 *  
#   mergedswe$swe:mergedswe$speciespipo -0.0236079  0.0061460  -3.841  0.00059 ***
#   mergedswe$swe:mergedswe$speciespotr  0.0008543  0.0054590   0.157  0.87669    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2805 on 30 degrees of freedom
# Multiple R-squared:  0.7065,	Adjusted R-squared:  0.6575 
# F-statistic: 14.44 on 5 and 30 DF,  p-value: 3.146e-07

AIC(swelm1) #18

sweci <- ggplot(data = mergedswe, aes(x=swe, y=diff, group=species, color=species)) + geom_point() + geom_smooth(method = "lm")  + 
  geom_smooth(method = "lm") + xlab("Climatic water sweicit") + ylab("Observed - Predicted RRW")
##Add R2 values to graph 
sweci <- sweci + 
  annotate(geom="text",x=25, y=1, label="Multiple R2 = .71")


#####Subset climate variable - swe 
sweclimateplot <- as.data.frame(sapply(split(combinedclimateplot$swe, combinedclimateplot$plotID),mean))
sweclimateplot <- cbind(sweclimateplot, unique(combinedclimateplot$plotID))
names(sweclimateplot)[1] <- "swe"
names(sweclimateplot)[2] <- "plotID"

#merge RRW difference and psdi climate mean for 1960-1989
mergedswe  <- left_join(fullmeandiff, sweclimateplot, by = 'plotID')

#Hard coding adding species column
mergedswe[1:4,4] <- 'juos'
mergedswe[5:9,4] <- 'pied'
mergedswe[10:16,4] <- 'pien'
mergedswe[17:24,4] <- 'pipo'
mergedswe[25:36,4] <- 'potr'
names(mergedswe)[4] <- "species"

#mixed effects model (need to separate by species, not plot)
swelm <- lm(mergedswe$diff ~ mergedswe$swe + (1/mergedswe$species))
summary(swelm)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.91393 -0.38183  0.03032  0.36513  0.77209 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)     -0.128329   0.295323  -0.435    0.667
# mergedswe$swe  0.001842   0.022599   0.082    0.936
# 
# Residual standard error: 0.4863 on 34 degrees of freedom
# Multiple R-squared:  0.0001954,	Adjusted R-squared:  -0.02921 
# F-statistic: 0.006644 on 1 and 34 DF,  p-value: 0.9355
AIC(swelm) #54

swelm1 <- lm(mergedswe$diff ~ mergedswe$swe + (mergedswe$swe/mergedswe$species))
summary(swelm1)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            -0.268620   0.200315  -1.341  0.18999    
# mergedswe$swe                         0.021755   0.014535   1.497  0.14492    
# mergedswe$swe:mergedswe$speciespied -0.004908   0.009798  -0.501  0.62006    
# mergedswe$swe:mergedswe$speciespien -0.047810   0.014386  -3.323  0.00235 ** 
#   mergedswe$swe:mergedswe$speciespipo -0.046487   0.009609  -4.838 3.68e-05 ***
#   mergedswe$swe:mergedswe$speciespotr  0.036760   0.011027   3.334  0.00229 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2376 on 30 degrees of freedom
# Multiple R-squared:  0.7893,	Adjusted R-squared:  0.7542 
# F-statistic: 22.48 on 5 and 30 DF,  p-value: 2.53e-09

AIC(swelm1) #6.14

sweci <- ggplot(data = mergedswe, aes(x=swe, y=diff, group=species, color=species)) + geom_point() + geom_smooth(method = "lm")  + 
  geom_smooth(method = "lm") + xlab("Snow water equivalent") + ylab("Observed - Predicted RRW")
##Add R2 values to graph 
sweci <- sweci + 
  annotate(geom="text",x=4, y=1, label="Multiple R2 = .79")





#####Subset climate variable - aet 
aetclimateplot <- as.data.frame(sapply(split(combinedclimateplot$aet, combinedclimateplot$plotID),mean))
aetclimateplot <- cbind(aetclimateplot, unique(combinedclimateplot$plotID))
names(aetclimateplot)[1] <- "aet"
names(aetclimateplot)[2] <- "plotID"

#merge RRW difference and psdi climate mean for 1960-1989
mergedaet  <- left_join(fullmeandiff, aetclimateplot, by = 'plotID')

#Hard coding adding species column
mergedaet[1:4,4] <- 'juos'
mergedaet[5:9,4] <- 'pied'
mergedaet[10:16,4] <- 'pien'
mergedaet[17:24,4] <- 'pipo'
mergedaet[25:36,4] <- 'potr'
names(mergedaet)[4] <- "species"

#mixed effects model (need to separate by species, not plot)
aetlm <- lm(mergedaet$diff ~ mergedaet$aet + (1/mergedaet$species))
summary(aetlm)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)   -0.92807    0.55589  -1.670    0.104
# mergedaet$aet  0.02437    0.01630   1.495    0.144
# 
# Residual standard error: 0.4711 on 34 degrees of freedom
# Multiple R-squared:  0.06171,	Adjusted R-squared:  0.03411 
# F-statistic: 2.236 on 1 and 34 DF,  p-value: 0.1441
AIC(aetlm) #52

aetlm1 <- lm(mergedaet$diff ~ mergedaet$aet + (mergedaet$aet/mergedaet$species))
summary(aetlm1)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         -0.309760   0.273413  -1.133   0.2662    
# mergedaet$aet                        0.012442   0.009207   1.351   0.1867    
# mergedaet$aet:mergedaet$speciespied -0.002803   0.004516  -0.621   0.5395    
# mergedaet$aet:mergedaet$speciespien -0.019326   0.004236  -4.562 8.00e-05 ***
#   mergedaet$aet:mergedaet$speciespipo -0.022187   0.004197  -5.287 1.04e-05 ***
#   mergedaet$aet:mergedaet$speciespotr  0.006807   0.003983   1.709   0.0978 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2181 on 30 degrees of freedom
# Multiple R-squared:  0.8226,	Adjusted R-squared:  0.793 
# F-statistic: 27.81 on 5 and 30 DF,  p-value: 2.041e-10

AIC(aetlm1) #-0.04

aetci <- ggplot(data = mergedaet, aes(x=aet, y=diff, group=species, color=species)) + geom_point() + geom_smooth(method = "lm")  + 
  geom_smooth(method = "lm") + xlab("Actual evapotranspiration") + ylab("Observed - Predicted RRW")
##Add R2 values to graph 
aetci <- aetci + 
  annotate(geom="text",x=30, y=1, label="Multiple R2 = .82")


#####Subset climate variable - aet 
aetclimateplot <- as.data.frame(sapply(split(combinedclimateplot$aet, combinedclimateplot$plotID),mean))
aetclimateplot <- cbind(aetclimateplot, unique(combinedclimateplot$plotID))
names(aetclimateplot)[1] <- "aet"
names(aetclimateplot)[2] <- "plotID"

#merge RRW difference and psdi climate mean for 1960-1989
mergedaet  <- left_join(fullmeandiff, aetclimateplot, by = 'plotID')

#Hard coding adding species column
mergedaet[1:4,4] <- 'juos'
mergedaet[5:9,4] <- 'pied'
mergedaet[10:16,4] <- 'pien'
mergedaet[17:24,4] <- 'pipo'
mergedaet[25:36,4] <- 'potr'
names(mergedaet)[4] <- "species"

#mixed effects model (need to separate by species, not plot)
aetlm <- lm(mergedaet$diff ~ mergedaet$aet + (1/mergedaet$species))
summary(aetlm)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.91393 -0.38183  0.03032  0.36513  0.77209 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)     -0.128329   0.295323  -0.435    0.667
# mergedaet$aet  0.001842   0.022599   0.082    0.936
# 
# Residual standard error: 0.4863 on 34 degrees of freedom
# Multiple R-squared:  0.0001954,	Adjusted R-squared:  -0.02921 
# F-statistic: 0.006644 on 1 and 34 DF,  p-value: 0.9355
AIC(aetlm) #54

aetlm1 <- lm(mergedaet$diff ~ mergedaet$aet + (mergedaet$aet/mergedaet$species))
summary(aetlm1)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            -0.268620   0.200315  -1.341  0.18999    
# mergedaet$aet                         0.021755   0.014535   1.497  0.14492    
# mergedaet$aet:mergedaet$speciespied -0.004908   0.009798  -0.501  0.62006    
# mergedaet$aet:mergedaet$speciespien -0.047810   0.014386  -3.323  0.00235 ** 
#   mergedaet$aet:mergedaet$speciespipo -0.046487   0.009609  -4.838 3.68e-05 ***
#   mergedaet$aet:mergedaet$speciespotr  0.036760   0.011027   3.334  0.00229 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2376 on 30 degrees of freedom
# Multiple R-squared:  0.7893,	Adjusted R-squared:  0.7542 
# F-statistic: 22.48 on 5 and 30 DF,  p-value: 2.53e-09

AIC(aetlm1) #6.14

aetci <- ggplot(data = mergedaet, aes(x=aet, y=diff, group=species, color=species)) + geom_point() + geom_smooth(method = "lm")  + 
  geom_smooth(method = "lm") + xlab("Actual evapotranspiration") + ylab("Observed - Predicted RRW")
##Add R2 values to graph 
aetci <- aetci + 
  annotate(geom="text",x=4, y=1, label="Multiple R2 = .79")


#Combine into one graph and export
pdf(file="regressCI.pdf", width=10,height=9, paper='special')
grid.arrange(pdsici, pptci, tmaxci, defci, sweci, aetci, nrow=3,ncol=2)
dev.off()




#######-----FINAL GRAPHS for FIG 4

#####Subset climate variable - def 
defclimateplot <- as.data.frame(sapply(split(combinedclimateplot$def, combinedclimateplot$plotID),mean))
defclimateplot <- cbind(defclimateplot, unique(combinedclimateplot$plotID))
names(defclimateplot)[1] <- "def"
names(defclimateplot)[2] <- "plotID"

#merge RRW difference and psdi climate mean for 1960-1989
mergeddef  <- left_join(fullmeandiff, defclimateplot, by = 'plotID')

#Hard coding adding species column
mergeddef[1:4,4] <- 'juos'
mergeddef[5:9,4] <- 'pied'
mergeddef[10:16,4] <- 'pien'
mergeddef[17:24,4] <- 'pipo'
mergeddef[25:36,4] <- 'potr'
names(mergeddef)[4] <- "species"

#mixed effects model (need to separate by species, not plot)
deflm <- lmer(mergeddef$diff ~ mergeddef$def + (1|mergeddef$species))
deflim1 <- lmer(mergeddef$diff ~ mergeddef$def + (mergeddef$def|mergeddef$species))

AIC(deflm)    #26
#AIC(deflm1)  #NA
r.squaredGLMM(deflm) 
# R2m       R2c
#  0.004134383 0.8773579
#r.squaredGLMM(deflm1)  #NA



#########If we want species names
mergeddef$species <- gsub('juos', 'Juniper', mergeddef$species)
mergeddef$species <- gsub('pien', 'Spruce', mergeddef$species)
mergeddef$species <- gsub('pipo', 'Ponderosa', mergeddef$species)
mergeddef$species <- gsub('potr', 'Aspen', mergeddef$species)
mergeddef$species <- gsub('pied', 'Pinyon', mergeddef$species)

#Rename column with capital letter for legend 
names(mergeddef)[4] <- "Species"

defci <- ggplot(data = mergeddef, aes(x=def, y=diff, group=Species, color=Species, linetype=Species)) + geom_point(size=2) +  
  geom_smooth(method = "lm", size=2) + 
  xlab("Climatic water deficit") + ylab("Observed - Predicted RRW") +
  scale_linetype_manual(values=c("dashed", "dashed", "solid", "solid", "solid")) +
  annotate(geom="text",x=80, y=1, label="A", size=8) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#A6761D")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 20))
defci


#############stacked bar plots for age/density 
setwd("/Users/nicolezenes/Library/Mobile Documents/com~apple~CloudDocs/FIA/")
#write.csv(mergedtmax, "mergedtmax.csv") get 5 year diff into master variable document
mastervar <- read.csv("MasterVariables1.csv")

#########If we want species names
mastervar$spcies <- gsub('JUOS', 'Juniper', mastervar$spcies)
mastervar$spcies <- gsub('PIEN', 'Spruce', mastervar$spcies)
mastervar$spcies <- gsub('PIPO', 'Ponderosa', mastervar$spcies)
mastervar$spcies <- gsub('POTR', 'Aspen', mastervar$spcies)
mastervar$spcies <- gsub('PIED', 'Pinyon', mastervar$spcies)

######Linear regression for RRW and stand density 2004#####
BA2m <- lmer(diff ~ totBA2004 + (1|spcies), data=mastervar)                            # Correct mixed effects model specification
BA2m1 <- lmer(diff ~ totBA2004 + (totBA2004 |spcies), data=mastervar)

AIC(BA2m)    # AIC better here: 35.9
AIC(BA2m1)    # AIC 32.2
r.squaredGLMM(BA2m) 
# R2m       R2c
# [1,] 0.004816051 0.7283478
r.squaredGLMM(BA2m1) 
# R2m       R2c
# [1,] 0.06214326 0.7942435

#Rename column with capital letter for legend 
names(mastervar)[3] <- "Species"

baplot <- ggplot(data = mastervar, aes(x = totBA2004, y=diff, color=Species, linetype=Species)) + geom_point(size=2)+ 
  geom_smooth(method = "lm", size = 2) +
  xlab("Total basal area 2004") +
  xlim(0,50) + ylab("Observed-Predicted RRW") + ylim(-1,1) +
  scale_linetype_manual(values=c("solid", "dashed", "solid", "solid", "solid")) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#A6761D")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 20)) + 
  annotate(geom="text",x=47, y=.95, label="B", size=8)
baplot

#-------------------Bill check 
age2m <- lmer(diff ~ avgage + (1|spcies), data=mastervar)    # Correct mixed effects model specification
age2m1 <- lmer(diff ~ avgage + (avgage |spcies), data=mastervar)

AIC(age2m)    # AIC better here: 36.6
AIC(age2m1)    # AIC 40.4
r.squaredGLMM(age2m) 
# R2m       R2c
# 0.002723935 0.7231787
r.squaredGLMM(age2m1) 
# [1,] 0.00341694 0.7248958

#Rename column with capital letter for legend 
names(mastervar)[3] <- "Species"

ageplot <- ggplot(data = mastervar, aes(x = avgage, y=diff, color=Species, linetype=Species)) + geom_point(size=2)+ 
  geom_smooth(method = "lm", size = 2) +
  xlab("Average stand age") + ylab("Observed - Predicted RRW") +
  scale_linetype_manual(values=c("dashed", "dashed", "dashed", "solid", "solid")) +
  annotate(geom="text",x=115, y=1.1, label="C", size=8) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#A6761D")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 20)) + xlim(30,120)
ageplot

grid.arrange(baplot, ageplot, nrow=2,ncol=1)


#########Dead ID subs
mastervar$deadID <- gsub('both', 'Both', mastervar$deadID)
mastervar$deadID <- gsub('diff', 'Different', mastervar$deadID)
mastervar$deadID <- gsub('none', 'None', mastervar$deadID)
mastervar$deadID <- gsub('self', 'Same', mastervar$deadID)

#Column name change
names(mastervar)[11] <- "DeadID"

stacked <- ggplot(mastervar, aes(fill=DeadID, y=change, x=Species)) + 
  geom_bar(position="stack", stat="identity") + ylab("Average percent difference") + 
  annotate(geom="text",x=5, y=600, label="D", size=8) +
  labs(color = "Species that died") + xlab("Species") +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 20)) 
stacked


##############Figure 4 - defci, basal area, age, deadID##############
pdf(file="combined4(def, BA, age, deadID.pdf", width=15,height=10, paper='special')

grid.arrange(arrangeGrob(defci, baplot, ageplot, stacked, nrow=2, ncol=2))
dev.off()        

#########Map of site locations for Figure 1##########
setwd("/Users/nicolezenes/Library/Mobile Documents/com~apple~CloudDocs/FIA/")
latlong <- read.csv("PlotLocationFuzzed_Grayson_NONM.csv")
points(latlong[,3], latlong[,2])

latlong$Species <- gsub('JUOS/PIED', 'Juniper/Pinyon', latlong$Species)
latlong$Species <- gsub('PIEN', 'Spruce', latlong$Species)
latlong$Species <- gsub('PIPO', 'Ponderosa', latlong$Species)
latlong$Species <- gsub('POTR', 'Aspen', latlong$Species)
latlong$Species <- gsub('PIED', 'Pinyon', latlong$Species)

####Figure 1 map
library(ggplot2)
library(ggspatial)
library(sf)
library(maps)

gusa <- map_data("state")
gusa <- subset(gusa, long < -100)
gusa <- subset(gusa, lat < 45)

gusa_sf <- st_as_sf(maps::map("state", regions = c("arizona", "california", "colorado", "idaho", "kansas", "montana", "nebraska", "nevada", "new mexico", "north dakota", "oklahoma", "oregon", "south dakota", "texas", "utah", "washington", "wyoming"), 
                              fill=TRUE, plot=FALSE))

gusa_proj <- st_transform(gusa_sf, crs = 5070)

latlong_sf <- st_as_sf(latlong, coords = c("LON_FUZZED", "LAT_FUZZED"), crs = 4326)

latlong_proj <- st_transform(latlong_sf, crs = 5070)

states_to_label <- c("utah", "colorado", "arizona", "new mexico")
label_data_manual <- data.frame(state_label = c("AZ", "CO", "NM", "UT"), long=c(-111,-107.5,-107.5,-111), lat=c(33,39.5,33,39.5))
label_sf <- st_as_sf(label_data_manual, coords=c("long", "lat"),crs=4326)
label_proj <- st_transform(label_sf,crs=5070)
mapfigure <- ggplot() +
  geom_sf(data = gusa_proj, fill = NA, color = "darkgrey") +
  geom_sf(data = latlong_proj, aes(color = Species)) +
  coord_sf() + # This ensures the plot uses the projected CRS
  geom_sf_text(data = label_proj, aes(label = state_label), color = "black", size = 3.5) +
  labs(x=NULL,y=NULL) +

    annotation_scale(
    location = "br", 
    style = "ticks", 
    width_hint = 0.2, 
    unit_category = "metric"
  ) +
  
  annotation_north_arrow(
    location = "bl", 
    pad_x = unit(0.2, "in"), 
    pad_y = unit(0.2, "in"), 
    style = north_arrow_fancy_orienteering
  ) +
  theme_minimal() +
  labs(color = "Species")

mapfigure

ggsave(file="Fig1MapFIAPlots.png", plot=mapfigure,units="in", width=6,height=4,dpi=600)

pdf(file="Fig1MapFIAPlots.pdf", width=6,height=4, paper='special')
mapfigure
dev.off()
