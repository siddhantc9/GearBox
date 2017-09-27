#install.packages(c("e1071","RPostgreSQL"))

library(e1071)
# Create a connection to the database
library(RPostgreSQL)
library(lubridate)


globVariables <- function(){
  
  basedir <<- "C:/Users/A672257/Siddhant/Workspace/Gearbox/dataset/";
  
  # Bearing data
  Bd <<- 0.331; # ball diameter, in inches
  Pd <<- 2.815; # pitch diameter, in inches
  Nb <<- 16; # number of rolling elements
  a <<- 15.17*pi/180; # contact angle, in radians
  s <<- 2000/60; # rotational frequency, in Hz
  
  ratio <<- Bd/Pd * cos(a);
  ftf <<- s/2 * (1 - ratio);
  bpfi <<- Nb/2 * s * (1 + ratio);
  bpfo <<- Nb/2 * s * (1 - ratio);
  bsf <<- Pd/Bd * s/2 * (1 - ratio**2);
}

# Helper functions
fft.spectrum <- function (d)
{
  fft.data <- fft(d)
  # Ignore the 2nd half, which are complex conjugates of the 1st half,
  # and calculate the Mod (magnitude of each complex number)
  return (Mod(fft.data[1:(length(fft.data)/2)]))
}

freq2index <- function(freq)
{
  step <- 10000/10240 # 10kHz over 10240 bins
  return (floor(freq/step))
}


all.features <- function(d)
{
  # Statistical features
  features <- c(quantile(d, names=FALSE), mean(d), sd(d), skewness(d), kurtosis(d))
  
  # RMS
  features <- append(features, sqrt(mean(d**2)))
  
  # Key frequencies
  fft.amps <- fft.spectrum(d)
  
  features <- append(features, fft.amps[freq2index(ftf)])
  features <- append(features, fft.amps[freq2index(bpfi)])
  features <- append(features, fft.amps[freq2index(bpfo)])
  features <- append(features, fft.amps[freq2index(bsf)])
  
  # Strongest frequencies
  n <- 5
  frequencies <- seq(0, 10000, length.out=length(fft.amps))
  sorted <- sort.int(fft.amps, decreasing=TRUE, index.return=TRUE)
  top.ind <- sorted$ix[1:n] # indexes of the largest n components
  features <- append(features, frequencies[top.ind]) # convert indexes to frequencies
  
  # Power in frequency bands
  vhf <- freq2index(6000):length(fft.amps)    # 6kHz plus
  hf <- freq2index(2600):(freq2index(6000)-1) # 2.6kHz to 6kHz
  mf <- freq2index(1250):(freq2index(2600)-1) # 1.25kHz to 2.6kHz
  lf <- 0:(freq2index(1250)-1)                # forcing frequency band
  
  powers <- c(sum(fft.amps[vhf]), sum(fft.amps[hf]), sum(fft.amps[mf]), sum(fft.amps[lf]))
  features <- append(features, powers)
  
  return(features)
}

vibrationsAsPerTimestamp <- function(){
  
  # Set up storage for bearing-grouped data
  b1m <- matrix(nrow=0, ncol=(2*23))
  b2m <- matrix(nrow=0, ncol=(2*23))
  b3m <- matrix(nrow=0, ncol=(2*23))
  b4m <- matrix(nrow=0, ncol=(2*23))
  # and for timestamps
  timestamp <- vector()
  
  i = 1
  for (filename in list.files(basedir)){
    cat(i, filename, "\n")
    i = i + 1
    
    ts <- as.character(strptime(filename, format="%Y.%m.%d.%H.%M.%S"))
    
    data <- read.table(paste0(basedir, filename), header=FALSE, sep="\t")
    colnames(data) <- c("b1.x", "b1.y", "b2.x", "b2.y", "b3.x", "b3.y", "b4.x", "b4.y")
    
    # Bind the new rows to the bearing matrices
    b1m <- rbind(b1m, c(all.features(data$b1.x), all.features(data$b1.y)))
    b2m <- rbind(b2m, c(all.features(data$b2.x), all.features(data$b2.y)))
    b3m <- rbind(b3m, c(all.features(data$b3.x), all.features(data$b3.y)))
    b4m <- rbind(b4m, c(all.features(data$b4.x), all.features(data$b4.y)))
    
    timestamp <- c(timestamp, ts)
  }
  
  timestamp <- suppressWarnings(changeDateFormat(timestamp))
  
  timestamp <- suppressWarnings(generateTimeSeries(timestamp))
  
  # Bind the new rows to the bearing matrices
  b1m <- data.frame(b1m, timestamp,"b1")
  b2m <- data.frame(b2m, timestamp,"b2")
  b3m <- data.frame(b3m, timestamp,"b3")
  b4m <- data.frame(b4m, timestamp,"b4")
  
  cnames <- c("Min.x", "Qu.1.x", "Median.x", "Qu.3.x", "Max.x", "Mean.x", "SD.x",
              "Skew.x", "Kurt.x", "RMS.x", "FTF.x", "BPFI.x", "BPFO.x", "BSF.x",
              "F1.x", "F2.x", "F3.x", "F4.x", "F5.x", "VHF.pow.x", "HF.pow.x",
              "MF.pow.x", "LF.pow.x", "Min.y", "Qu.1.y", "Median.y", "Qu.3.y", 
              "Max.y", "Mean.y", "SD.y", "Skew.y", "Kurt.y", "RMS.y", "FTF.y", 
              "BPFI.y", "BPFO.y", "BSF.y", "F1.y", "F2.y", "F3.y", "F4.y", 
              "F5.y", "VHF.pow.y", "HF.pow.y", "MF.pow.y", "LF.pow.y", 
              "Timestamp","Bearing")
  
  colnames(b1m) <- cnames
  colnames(b2m) <- cnames
  colnames(b3m) <- cnames
  colnames(b4m) <- cnames
  
  allData <- rbind(b1m,b2m,b3m,b4m)
    
  return(allData)
}

insertInDatabase <- function(data = NULL){
  
  ## Loading required package: DBI
  pg = dbDriver("PostgreSQL")
  
  con = dbConnect(pg, user="postgres", password="postgres",
                  host="localhost", port=5432, dbname="gearBox")
  
  dbWriteTable(con,'gear_box',data, row.names=FALSE)
  
  # commit the change
  dbCommit(con)
  
  # disconnect from the database
  dbDisconnect(con)
  
}

dbSafeNames = function(names) {
  names = gsub('[^a-z0-9]+','_',tolower(names))
  names = make.names(names, unique=TRUE, allow_=TRUE)
  names = gsub('.','_',names, fixed=TRUE)
  return(names)
}

#Function to identify timestamp format
changeDateFormat <- function(timestamp = NULL){
  if(is.na(dmy_hms(timestamp)) & is.na(mdy_hms(timestamp)))
  {
    formatt ="ymd_hms"
    timestamp <- gsub("-","/",ymd_hms(timestamp))
  }else if(is.na(mdy_hms(timestamp)) & is.na(ymd_hms(timestamp)))
  {
    formatt = "dmy_hms"
    timestamp <- gsub("-","/",dmy_hms(timestamp))
  }else
  {
    formatt = "mdy_hms"
    timestamp <- gsub("-","/",mdy_hms(timestamp))
  }
}

### Timedifference,FINDMIN,FINDMAX
minmaxdatetime<-function(timestamp = NULL){
  
  dataset <- strptime(timestamp, "%Y/%m/%d %H:%M:%S")
  #dataset <- lapply((timestamp), strptime, "%Y/%m/%d %H:%M:%S")[1]
  vecNa <- !is.na(timestamp)
  dataset <- dataset[vecNa]
  #minimumdate
  mindate<-min(as.POSIXlt(dataset))
  #maxdate
  maxdate<-max(as.POSIXlt(dataset))
  
  
  #code for time difference
  N<-length(dataset)
  
  time1<-as.POSIXct(dataset[2:N])
  time2<-as.POSIXct(dataset[1:(N-1)])
  timedata<-data.frame(difftime(time1 = time1,time2 = time2,units="secs"))
  names(timedata)<-"timediff"
  
  temp<-data.frame(timedata$timediff[timedata$timediff!=0])
  
  names(temp)<-"timediffsort"
  print(sort(table((temp$timediffsort)),decreasing = TRUE)[1:10])
  timediff<-names(sort(table((temp$timediffsort)),decreasing = TRUE)[1:4])[1]
  timediff<-as.numeric(timediff)
  #single dataset for minimum date,maximum date,time difference
  mindate <- as.POSIXct(gsub("-","/",mindate))
  maxdate <- as.POSIXct(gsub("-","/",maxdate))
  minmaxtime<-data.frame(mindate,maxdate,timediff)
  return(minmaxtime)
  
}

generateTimeSeries <- function(timestamp = NULL){
  
  if(!is.null(timestamp)){
    minmaxdifftime <- minmaxdatetime(timestamp)
    startTS <- minmaxdifftime$mindate
    endTS <- minmaxdifftime$maxdate
    timeDifferenceInSeconds <- minmaxdifftime$timediff
    MainTSList <- list();
    i=1
    while(length(timestamp)!=(i-1)) {
      MainTSList[[i]] <- startTS
      i <- i + 1
      startTS=startTS+as.numeric(timeDifferenceInSeconds)
    }
    TSDf=do.call(rbind, lapply(MainTSList, data.frame, stringsAsFactors=FALSE))
    colnames(TSDf)=c("timestamp")
    
    return(TSDf)
  }
  
}

mainFunction <- function(){
  globVariables()
  finalData <- vibrationsAsPerTimestamp() #Localize the scope
  colnames(finalData) <- dbSafeNames(colnames(finalData))
  finalData <<- finalData;
  
  # #Invoke identify_format function
  # finalData$timestamp <- suppressWarnings(changeDateFormat(finalData))
  # 
  # timeSeries <- suppressWarnings(generateTimeSeries(data = finalData))
  # 
  # finalData2 <<- timeSeries;
  
  #insertInDatabase(finalData)
}

