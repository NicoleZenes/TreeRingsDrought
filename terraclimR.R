#Ecosphere
#Article
#ECS25-0062 - Complex trajectories of tree growth in the southwestern United States after severe drought
#Nicole Zenes, Leander D. L. Anderegg, Kiona Ogle, Drew M. P. Peltier, and William R. L.Anderegg

#Abatzoglou, J., Dobrowski, S., Parks, S. et al. TerraClimate, a high-resolution global dataset of monthly climate 
#and climatic water balance from 1958–2015. Sci Data 5, 170191 (2018). https://doi.org/10.1038/sdata.2017.191

#Libraries used
library(ncdf4)
library(MASS)
library(doParallel)
library(ncdf4)


###########################################################################################
setwd("/Users/nicolezenes/Library/Mobile Documents/com~apple~CloudDocs/FIA/")
location <- read.csv("PlotLocationFuzzed_Grayson.csv")
plot=as.data.frame(seq(1, 34, by=1))
variable=as.data.frame(seq(1,8, by=1))

#Climate variables to extract - actual evapotrans, cwd, soil moisture, snow water equivalent
#PDSI, VPD, max monthly temp, precip accumulation
period <- as.data.frame(c("aet","def","soil","swe","PDSI","vpd","tmax","ppt","srad"))

foreach (p = iter(variable, by ="row")) %:%
  foreach(i = iter(plot, by = "row"), .packages='MASS') %do% {
    # enter in variable you want to download see: http://thredds.northwestknowledge.net:8080/thredds/terraclimate_aggregated.html
    var = as.character(period[p,1])
    
    # enter in longitude, latitude here
    lat = location[i,2]
    long = location[i,3]
    x<-c(long, lat)

    baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1958_CurrentYear_GLOBE.nc")
    
    nc <- nc_open(baseurlagg)
    lon <- ncvar_get(nc, "lon")
    lat <- ncvar_get(nc, "lat")
    flat = match(abs(lat - x[2]) < 1/48, 1)
    latindex = which(flat %in% 1)
    flon = match(abs(lon - x[1]) < 1/48, 1)
    lonindex = which(flon %in% 1)
    start <- c(lonindex, latindex, 1)
    count <- c(1, 1, -1)
    
    
    # read in the full period of record using aggregated files
    data <- as.numeric(ncvar_get(nc, varid = var,start = start, count))
    
    #Write to csv with var and plot ID
    write.csv(data, paste("/Users/nicolezenes/Library/Mobile Documents/com~apple~CloudDocs/FIA/Terraclim/", var,location[i,1],".csv"))
    
  }

#### Combine all climate variables per plot into one spreadsheet 
setwd("/Users/nicolezenes/Documents/ONEFlux/terraclim/")
#setwd("/Users/nicolezenes/Library/Mobile Documents/com~apple~CloudDocs/FIA/")
location <- read.csv("latlong.csv")

#Climate variables to extract - actual evapotrans, cwd, soil moisture, snow water equivalent
#PDSI, VPD, max monthly temp, precip accumulation
period <- as.data.frame(c("aet","def","soil","swe","PDSI","vpd","tmax","ppt"))
fullvar <- as.data.frame(matrix(0, ncol=745, nrow=25))
fullvar[,1] <- location[,1]

#For loops to go through every plot for each climate variable 
for (p in (1:8)) {
  for(i in (1:34)) {
    # enter in variable 
    var <- as.character(period[p,1])
    # plot ID
    plot <- location[i,1]
    
    #read in csv
    climate <- read.csv(paste(
      "/Users/nicolezenes/Library/Mobile Documents/com~apple~CloudDocs/FIA/Terraclim/", var,location[i,1],".csv"))
    fullvar[i, 2:745] <- t(climate[,2])
    
    #If all plots have been read in, write to csv
    if (i > 33) {
      ##Remove 2014-2019
      fullvar <- fullvar[,1:673]
      
      write.csv(fullvar, paste(
        "/Users/nicolezenes/Library/Mobile Documents/com~apple~CloudDocs/FIA/Terraclim/", var,"full.csv"))
    }
  }
}


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



  
  
  
  
