library(ncdf4)
library(tidyverse)
library(oxcAAR)
library(here)
library(lipdR)
library(purrr)

lldb <- nc_open("lake_level_database.nc")
vn <- attributes(lldb$var)$names

lldbl <- list()
#pull all the data and put it in a list
for(v in vn){
  lldbl[[v]] <- ncvar_get(lldb, v)
}

#these are pretty metadata sparse, so they will be slim LiPD files.

#let's also load in some relevant data from the text files.
lakedata <- read_lines("dataFromWDS/lakedata.txt")
numrefs <- read_lines("dataFromWDS/numrefs.txt")

pubindex <- read_delim("dataFromWDS/pubindex.txt",skip = 6,delim = "    ",col_names = c("Index","CitationIndex"))

#oof. This structure is a little rough.
sn <- lldbl$basin_name
D <- list()#initialize list of LiPD files

#Let's calibrate those radiocarbon dates. 
getCalMedian <- function(c14,unc){
  calout <- oxcAAR::oxcalCalibrate(c14,unc,)
  cp <- cumsum(calout$`1`$raw_probabilities$probabilities)
  cp <- cp/max(cp)
  calAge <- geoChronR::convertAD2BP(calout$`1`$raw_probabilities$dates[which(abs(cp-.5)==min(abs(cp-.5)))])
  return(calAge)
}

unc <- 200
calBinYoung <- purrr::map_dbl(lldbl$age_ka*1000,getCalMedian,unc)
calBinMid <- purrr::map_dbl(lldbl$age_ka*1000+500,getCalMedian,unc)
calBinOld <- purrr::map_dbl(lldbl$age_ka*1000+1000,getCalMedian,unc)

ageStarter <- list()
ageStarter$variableName <- "age"
ageStarter$units <- "BP"
ageStarter$variableType <- "inferred"
ageStarter$inferredVariableType <- "age"
ageStarter$description <- "estimated calibrated bin mid-point age"
ageStarter$notes <- "These data are estimated from what were 1000 14C year intervals, calibrated using OxCal v4.4.2 and IntCal20"
ageStarter$values <- as.matrix(calBinMid)

ageStarterYoung <- ageStarter
ageStarter$variableName <- "ageYoung"
ageStarter$description <- "estimated calibrated bin young age"
ageStarter$values <- as.matrix(calBinYoung)

ageStarterOld <- ageStarter
ageStarter$variableName <- "ageOld"
ageStarter$description <- "estimated calibrated bin old age"
ageStarter$values <- as.matrix(calBinOld)



for(i in 1:length(sn)){
  L <- list()#initialize new LiPD file
  
  bmn <- lldbl$brown_number[i]
  
  #geo
  geo <- list()
  geo$latitude <- lldbl$latitude[i]
  geo$longitude <- lldbl$longitude[i]
  geo$description <- sn[i]
  geo$siteName <- stringr::str_extract(sn[i],pattern = "^([^,])+") %>% 
    tolower() %>% 
    tools::toTitleCase() %>% 
    str_remove_all("[^A-Za-z0-9]")
  
  L$geo <- geo
  
  #pub
  #find site in pubindex
  tp <- pubindex$CitationIndex[which(pubindex$Index == bmn)] %>% 
    str_trim() %>% 
    str_split(",") %>% 
    unlist()
  
  if(!is.null(tp)){
    tp <- 0
  }
  
  pub <- vector(mode = "list",length = length(tp)+1)
  
  pub[[1]]$citation <- "Street-Perrott, F.A., D.S. Marchand, N. Roberts, and S.P Harisson, 1989. Global Lake-Level Variations from 18,000 to 0 Years Ago: A Paleoclimatic Analysis.   U.S. Department of Energy Technical Report 46, Washington, D.C. 20545.   Distributed by National Technical Information Service, Springfield, VA 22161."
  
  
  if(!all(tp==0)){
    for(p in 2:length(pub)){
      wp <- which(startsWith(numrefs,tp[p-1]))
      pub[[p]]$citation <- str_remove(numrefs[wp],str_c(tp[p-1]," "))
    }
  }
  
  L$pub <- pub 
  
  
  #get primary reference
  pfi <- which(startsWith(lakedata,sn[i]))
  pr <- str_extract(lakedata[pfi+1],pattern = "(?<=>)(.*?)(?=<)") %>% 
    str_trim()
  if(is.na(pr)){#try for only one author
    pr <- str_extract(lakedata[pfi+1],pattern = "[>].*") %>% 
      str_remove_all(">") %>% 
      str_trim()
  }
  
  
  spl <- str_split(pr,",")
  firstAuthor <- stringr::str_extract(spl[[1]][1],pattern = "^([^ ])+") %>% 
    tolower() %>% 
    tools::toTitleCase() %>% 
    str_remove_all("[^A-Za-z0-9]")
  
  pubYear <- spl[[1]][2]
  
  #root
  L$dataSetName <- str_c(geo$siteName,".",firstAuthor,".",pubYear) %>% str_remove_all(" ")
  print(L$dataSetName)
  L$lipdVersion <- 1.3
  L$datasetId <- lipdverseR::createDatasetId()
  L$originalDataUrl <- "ftp://ftp.ncdc.noaa.gov/pub/data/paleo/paleolimnology/lakelevels/oxford/"
  L$createdBy <- "http://github.com/nickmckay/oxforLakeStatus2Lipd"
  L$archiveType <- "LakeDeposits"
  #Paleodata
  mt <- vector(mode = "list",length = 1)
  #ages first
  mt[[1]]$age <- ageStarter
  mt[[1]]$age$TSid <- str_c("OLS_",lipdR::createTSid())
  
  mt[[1]]$ageYoung <- ageStarterYoung
  mt[[1]]$ageYoung$TSid <- str_c("OLS_",lipdR::createTSid())
  
  mt[[1]]$ageOld <- ageStarterOld
  mt[[1]]$ageOld$TSid <- str_c("OLS_",lipdR::createTSid())
  
  mt[[1]]$lakeStatus <- list()
  mt[[1]]$lakeStatus$variableName <- "lakeStatus"
  mt[[1]]$lakeStatus$TSid <- str_c("OLS_",lipdR::createTSid())
  mt[[1]]$lakeStatus$units <- "unitless index"
  mt[[1]]$lakeStatus$variableType <- "inferred"
  mt[[1]]$inferredVariableType <- "lakeLevel"
  mt[[1]]$lakeStatus$values <- as.matrix(lldbl$lake_status[,i])
  mt[[1]]$lakeStatus$description <- "This is a millennial-scale estimate of lake status in relative units. 3 corresponds to LOW: 0-15% of the total altitudinal range of fluctuation, including dry lakes, 2 corresponds to INTERMEDIATE, or 15-70% of the total range, and 1 corresponds to HIGH: 70-100% of the total range of fluctuation, including overflowing lakes. 0 is uncodable."
  mt[[1]]$lakeStatus$interpretation <- vector(mode = "list",length = 1)
  mt[[1]]$lakeStatus$interpretation[[1]]$variable <- "P-E"
  mt[[1]]$lakeStatus$interpretation[[1]]$variableDetail <- "lakeLevel@surface"
  mt[[1]]$lakeStatus$interpretation[[1]]$scope <- "climate"
  mt[[1]]$lakeStatus$interpretation[[1]]$direction <- "negative"
  
  icb <- list(compilationName = "OxfordLSDB-LiPD",compilationVersion = "1_0_0")
  mt[[1]]$lakeStatus$inCompilationBeta <- list(icb)
  
  
  #lake trend
  mt[[1]]$lakeTrend <- list()
  mt[[1]]$lakeTrend$variableName <- "lakeTrend"
  mt[[1]]$lakeTrend$TSid <- str_c("OLS_",lipdR::createTSid())
  mt[[1]]$lakeTrend$units <- "unitless index"
  mt[[1]]$lakeTrend$values <- as.matrix(lldbl$lake_trend[,i])
  mt[[1]]$lakeTrend$variableType <- "inferred"
  mt[[1]]$lakeTrend$description <- "This is a millennial-scale estimate of lake trend in relative units. 3 corresponds to falling, 2 corresponds to stable, and 1 corresponds to rising. 0 is uncodable."
  
  L$paleoData <- vector(mode = "list",length = 1)
  L$paleoData[[1]]$measurementTable <- mt
  
  
  #initialize the changelog
  L <- lipdverseR::initializeChangelog(L,notes = "Created from the Oxford Lake Status Database. See github.com/nickmckay/oxfordLakeStatus2Lipd")
  
  
  #chronData. There are extensive chronData here, but it's quite messy and so I'm putting this off for now
  writeLipd(L,here("lipds"))
  
  
}
