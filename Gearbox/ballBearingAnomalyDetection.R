basedir <- "C:/Users/A672257/Downloads/1st_test/"

library(e1071)
library("nnet")

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

# Bearing data
Bd <- 0.331 # ball diameter, in inches
Pd <- 2.815 # pitch diameter, in inches
Nb <- 16 # number of rolling elements
a <- 15.17*pi/180 # contact angle, in radians
s <- 2000/60 # rotational frequency, in Hz

ratio <- Bd/Pd * cos(a)
ftf <- s/2 * (1 - ratio)
bpfi <- Nb/2 * s * (1 + ratio)
bpfo <- Nb/2 * s * (1 - ratio)
bsf <- Pd/Bd * s/2 * (1 - ratio**2)


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


# Set up storage for bearing-grouped data
b1m <- matrix(nrow=0, ncol=(2*23))
b2m <- matrix(nrow=0, ncol=(2*23))
b3m <- matrix(nrow=0, ncol=(2*23))
b4m <- matrix(nrow=0, ncol=(2*23))
# and for timestamps
timestamp <- vector()

i = 1
for (filename in head(list.files(basedir))){
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

# Bind the new rows to the bearing matrices
b1m <- cbind(b1m, timestamp)
b2m <- cbind(b2m, timestamp)
b3m <- cbind(b3m, timestamp)
b4m <- cbind(b4m, timestamp)

cnames <- c("Min.x", "Qu.1.x", "Median.x", "Qu.3.x", "Max.x", "Mean.x", "SD.x", "Skew.x", "Kurt.x", "RMS.x", "FTF.x", "BPFI.x", "BPFO.x", "BSF.x", "F1.x", "F2.x", "F3.x", "F4.x", "F5.x", "VHF.pow.x", "HF.pow.x", "MF.pow.x", "LF.pow.x", "Min.y", "Qu.1.y", "Median.y", "Qu.3.y", "Max.y", "Mean.y", "SD.y", "Skew.y", "Kurt.y", "RMS.y", "FTF.y", "BPFI.y", "BPFO.y", "BSF.y", "F1.y", "F2.y", "F3.y", "F4.y", "F5.y", "VHF.pow.y", "HF.pow.y", "MF.pow.y", "LF.pow.y", "Timestamp")

colnames(b1m) <- cnames
colnames(b2m) <- cnames
colnames(b3m) <- cnames
colnames(b4m) <- cnames
b1 <- data.frame(b1m)
b2 <- data.frame(b2m)
b3 <- data.frame(b3m)
b4 <- data.frame(b4m)

### Add an extra column for classes of bearing state

# Bearing state classes (determined by eye)
b1.labels <- c(rep("early", times=150), rep("unknown", times=450), rep("normal", times=899), rep("suspect", times=600), rep("failure.b1", times=57))
b2.labels <- c(rep("early", times=499), rep("normal", times=1500), rep("suspect", times=120), rep("failure.b2", times=37))
b3.labels <- c(rep("early", times=499), rep("normal", times=1290), rep("suspect", times=330), rep("failure.inner", times=37))
b4.labels <- c(rep("early", times=199), rep("normal", times=800), rep("suspect", times=435), rep("failure.roller", times=405), rep("stage2", times=317))

b1 <- cbind(b1, State=b1.labels)
b2 <- cbind(b2, State=b2.labels)
b3 <- cbind(b3, State=b3.labels)
b4 <- cbind(b4, State=b4.labels)

# Now make one large "normal" dataset
# Ignore columns 47 (timestamp), 16 and 39 (F1.x and F1.y), and 48 (state), which are not useful for correlation
norm <- rbind(
  cbind((b1[(b1$State == "normal"),-c(47, 16, 39, 48)]), bearing="b1"),
  cbind((b2[(b2$State == "normal"),-c(47, 16, 39, 48)]), bearing="b2"),
  cbind((b3[(b3$State == "normal"),-c(47, 16, 39, 48)]), bearing="b3"),
  cbind((b4[(b4$State == "normal"),-c(47, 16, 39, 48)]), bearing="b4")
)

# Export
write.table(norm, file=paste0(basedir, "../normal_bearings.csv"), sep=",", row.names=FALSE)

library(corrgram)
corrgram(norm)

# Find high (absolute) correlations (ignoring non-numeric bearing name col)
cor <- cov2cor(cov(norm[, c(-44,-45)]))
alikes <- apply(cor, 2, function(col) { names(Filter(function (val) { val > 0.9 }, sort.int(abs(col), decreasing=TRUE))) } )

cat(str(alikes, vec.len=10))
# As an offline exercise, minimise the set of features by removing those with high correlations

# I found offline analysis to give the following minimal set
best <- c("Min.x", "Median.x", "Max.x", "Mean.x", "Skew.x", "Kurt.x", "FTF.x", "BPFI.x", "BPFO.x", "BSF.x", "F2.x", "F3.x", "F4.x", "F5.x", "Min.y", "Max.y", "Skew.y", "Kurt.y", "FTF.y", "BPFI.y", "BPFO.y", "BSF.y", "F2.y", "F3.y", "F4.y", "F5.y", "Qu.1.x", "VHF.pow.x", "Qu.1.y", "Median.y", "HF.pow.y")

# Select all rows from all bearings, only the best feature cols plus State, and bind on bearing
data <- rbind(
  cbind(bearing="b1", (b1[,c("State", best)])),
  cbind(bearing="b2", (b2[,c("State", best)])),
  cbind(bearing="b3", (b3[,c("State", best)])),
  cbind(bearing="b4", (b4[,c("State", best)]))
)

write.table(data, file=paste0(basedir, "../all_bearings.csv"), sep=",", row.names=FALSE)


library("entropy")
# MI = H(x) + H(y) - H(x, y)
H.x <- entropy(table(data$State))
mi <- apply(data[, -c(1, 2)], 2, function(col) { H.x + entropy(table(col)) - entropy(table(data$State, col))})
sort(mi, decreasing=TRUE)

# Now select the correct size of FV by testing different classifiers
library("rpart")
library("caret")

train.rows <- seq(1, length(data$State), by=4) # every fourth index
sorted <- names(sort(mi, decreasing=TRUE))
accuracies <- vector()

for (i in 2:length(sorted)){
  form <- paste0("data$State ~ data$", Reduce(function (r, name) { paste(r, paste0("data$", name), sep=" + ") }, sorted[1:i]))
  model <- rpart(as.formula(form), data, subset=train.rows)
  pred <- predict(model, data, type="class") # contains train and test set predictions	
  cm <- confusionMatrix(pred[-train.rows], data[-train.rows, 2]) # only the non-training rows
  # pull out the answer
  accuracies <- c(accuracies, cm["overall"][[1]]["Accuracy"][[1]])
}

plot(accuracies, xlab="Number of features", ylab="Accuracy")

# 14 features looks best. save those plus bearing and state labels
write.table(data[c("bearing", "State", sorted[1:14])], file=paste0(basedir, "../all_bearings_best_fv.csv"), sep=",", row.names=FALSE)


calculate.confusion <- function(states, clusters)
{
  # generate a confusion matrix of cols C versus states S
  d <- data.frame(state = states, cluster = clusters)
  td <- as.data.frame(table(d))
  # convert from raw counts to percentage of each label
  pc <- matrix(ncol=max(clusters),nrow=0) # k cols
  for (i in 1:9) # 9 labels
  {
    total <- sum(td[td$state==td$state[i],3])
    pc <- rbind(pc, td[td$state==td$state[i],3]/total)
  }
  rownames(pc) <- td[1:9,1]
  return(pc)
}

assign.cluster.labels <- function(cm, k)
{
  # take the cluster label from the highest percentage in that column
  cluster.labels <- list()
  for (i in 1:k)
  {
    cluster.labels <- rbind(cluster.labels, row.names(cm)[match(max(cm[,i]), cm[,i])])
  }
  
  # this may still miss some labels, so make sure all labels are included
  for (l in rownames(cm)) 
  { 
    if (!(l %in% cluster.labels)) 
    { 
      cluster.number <- match(max(cm[l,]), cm[l,])
      cluster.labels[[cluster.number]] <- c(cluster.labels[[cluster.number]], l)
    } 
  }
  return(cluster.labels)
}

calculate.accuracy <- function(states, clabels)
{
  # For each means$cluster c, if cluster.labels[[c]] contains data$State, it's correct
  matching <- Map(function(state, labels) { state %in% labels }, states, clabels)
  tf <- unlist(matching, use.names=FALSE)
  return (sum(tf)/length(tf))
}



# Read in best features
basedir <- "C:/Users/A672257/Downloads/1st_test/"
data <- read.table(file=paste0(basedir, "../all_bearings_best_fv.csv"), sep=",", header=TRUE)

# Test multiple k values
results <- matrix(ncol=2, nrow=0)
models <- list()

for (k in 6:18){
  # Don't cluster columns for bearing or State  
  means <- kmeans(data[,-c(1, 2)], k)
  
  # generate a confusion matrix of cols C versus states S
  conf.mat <- calculate.confusion(data$State, means$cluster)
  cluster.labels <- assign.cluster.labels(conf.mat, k)
  
  # Now calculate accuracy, using states and groups of labels for each cluster
  accuracy <- calculate.accuracy(data$State, cluster.labels[means$cluster])
  results <- rbind(results, c(k, accuracy))
  models[[(length(models)+1)]] <- means
}


# Can run the loop above multiple times to get different random starts
# Then pick out the most accurate
best.row <- match(max(results[,2]), results[,2])
best.kmeans <- models[[best.row]]
save(best.kmeans, file=paste0(basedir, "../../models/kmeans.obj"))


# Now visualise
k <- length(best.kmeans$size)
conf.mat <- calculate.confusion(data$State, best.kmeans$cluster)
cluster.labels <- assign.cluster.labels(conf.mat, k)
acc <- calculate.accuracy(data$State, cluster.labels[best.kmeans$cluster])
cat("For", k, "means with accuracy", acc, ", labels are assigned as:\n")
cat(str(cluster.labels))

# use the same colours for states as before
cols <- do.call(rbind, Map(function(s)
{
  if (s=="early") "green" 
  else if (s == "normal") "blue" 
  else if (s == "suspect") "darkgoldenrod" 
  else if (s == "stage2") "salmon"
  else if (s == "unknown") "black"
  else "red"
}, data$State))


# plot each bearing changing state
par(mfrow=c(2,2))
for (i in 1:4){
  s <- (i-1)*2156 + 1 # 2156 datapoints per bearing
  e <- i*2156
  plot(best.kmeans$cluster[s:e], col=cols[s:e], ylim=c(1,k), main=paste0("Bearing ", i), ylab="Cluster")
}

# plot the data with cluster centres superimposed
par(mfrow=c(2,3))
plot(Kurt.x ~ Skew.x, data=data, col=cols)
points(Kurt.x ~ Skew.x, data=best.kmeans$centers, pch=10)
plot(FTF.x ~ BPFI.x, data=data, col=cols)
points(FTF.x ~ BPFI.x, data=best.kmeans$centers, pch=10)
plot(BPFO.x ~ BSF.x, data=data, col=cols)
points(BPFO.x~ BSF.x, data=best.kmeans$centers, pch=10)
plot(FTF.y ~ BPFI.y, data=data, col=cols)
points(FTF.y ~ BPFI.y, data=best.kmeans$centers, pch=10)
plot(BPFO.y ~ BSF.y, data=data, col=cols)
points(BPFO.y ~ BSF.y, data=best.kmeans$centers, pch=10)
plot(VHF.pow.x ~ HF.pow.y, data=data, col=cols)
points(VHF.pow.x ~ HF.pow.y, data=best.kmeans$centers, pch=10)

data <- read.table(file=paste0(basedir, "../all_bearings_best_fv.csv"), sep=",", header=TRUE)

# Read in the best kmeans model. Straight after running kmeans.R the 
# filename will be "kmeans.obj", but this is the best model I found.
load(paste0(basedir, "../../models/best-kmeans-12-0.612.obj"))

## Adjust the class labels as a result of k-means
# Cluster 2 should be labelled "suspect"
data[best.kmeans$cluster==2,2] <- "suspect"

# Cluster 3 should be labelled "normal"
data[best.kmeans$cluster==3,2] <- "normal"

# Cluster 9 datapoints labelled "early" should be "normal"
data[((best.kmeans$cluster==9)&(data$State=="early")),2] <- "normal"

# b1 failure looks like a rolling element failure
data[data$State=="failure.b1", 2] <- "failure.roller"
data[((best.kmeans$cluster==11)&(data$State=="unknown")),2] <- "failure.roller"
data[((best.kmeans$cluster==11)&(data$State=="normal")),2] <- "suspect"
data[((best.kmeans$cluster==10)&(data$State=="unknown")),2] <- "early"

# Cluster 6 should all be "normal"
data[best.kmeans$cluster==6,2] <- "normal"

## Now plot to check the result
# use the same colours for states as before
cols <- do.call(rbind, Map(function(s)
{
  if (s=="early") "green" 
  else if (s == "normal") "blue" 
  else if (s == "suspect") "darkgoldenrod" 
  else if (s == "stage2") "salmon"
  else if (s == "unknown") "black"
  else "red"
}, data$State))

# plot each bearing changing state
par(mfrow=c(2,2))
for (i in 1:4){
  s <- (i-1)*2156 + 1 # 2156 datapoints per bearing
  e <- i*2156
  plot(best.kmeans$cluster[s:e], col=cols[s:e], ylim=c(1,k), main=paste0("Bearing ", i), ylab="Cluster")
}

# Now save
write.table(data, file=paste0(basedir, "../all_bearings_relabelled.csv"), sep=",", row.names=FALSE)

# Read in the relabelled best features
basedir <- "C:/Users/A672257/Downloads/1st_test/"
data <- read.table(file=paste0(basedir, "../all_bearings_relabelled.csv"), sep=",", header=TRUE)

# Split into train and test sets, preserving percentage across states
train.pc <- 0.7
train <- vector()
test <- vector()
for (state in unique(data$State)){
  all.samples <- data[data$State==state,]
  len <- length(all.samples[,1])
  rownums <- sample(len, len*train.pc, replace=FALSE)
  train <- c(train, as.integer(row.names(all.samples)[rownums]))
  test <- c(test, as.integer(row.names(all.samples)[-rownums]))
}

# Calculate accuracy weighted by counts per class
weighted.acc <- function(predictions, actual)
{
  freqs <- as.data.frame(table(actual))
  tmp <- t(mapply(function (p, a) { c(a, p==a) }, predictions, actual, USE.NAMES=FALSE)) # map over both together
  tab <- as.data.frame(table(tmp[,1], tmp[,2])[,2]) # gives rows of [F,T] counts, where each row is a state
  acc.pc <- tab[,1]/freqs[,2]
  return(sum(acc.pc)/length(acc.pc))
}

# Set up class weights to penalise the minority classes more
cw1 <- rep(1, 8) # all equal
cw2 <- c(10, 100, 100, 10, 1, 10, 1,10) # 1/order of count

freqs <- as.data.frame(table(data$State))
cw3 <- cbind(freqs[1], apply(freqs, 1, function(s) { length(data[,1])/as.integer(s[2])})) # 1/weight

class.weights <- rbind(cw1, cw2, cw3[,2])
colnames(class.weights) <- c("early", "failure.b2", "failure.inner", "failure.roller", "normal", "stage2", "suspect","unknown")

# Also normalise the data for comparison
normed <- cbind(data[,1:2], as.data.frame(lapply(data[,-c(1,2)], function(col) { col / max(abs(col)) })))

results <- matrix(ncol=6, nrow=0)
models <- list()

# Run three iterations of each
for (i in 1:3)
{
  for (c in 1:length(class.weights[,1]))
  {
    data.weights <- do.call(rbind, Map(function(s)
    {
      class.weights[c,s]
    }, data$State))
    
    for (h in 2:30)
    {
      cat("Run", i, "for c", c, "and h", h, "\n")
      # With range
      ann <- nnet(State ~ ., data=data[train,-1], weights=data.weights[train], size=h, decay=5e-4, rang=(1/max(data[,-c(1,2)])), maxit=200)
      pred <- predict(ann, data[,-1], type="class")
      tacc <- weighted.acc(pred[train], data[train,2])
      wacc <- weighted.acc(pred[-train], data[-train,2])
      pacc <- sum(pred[-train]==data[-train,2])/length(pred[-train])
      
      results <- rbind(results, c(h, tacc, wacc, pacc, c, 1))
      models[[(length(models)+1)]] <- ann
      
      # With normalised data (no need for range now)
      ann <- nnet(State ~ ., data=normed[train,-1], weights=data.weights[train], size=h, decay=5e-4, maxit=200)
      pred <- predict(ann, normed[,-1], type="class")
      tacc <- weighted.acc(pred[train], normed[train,2])
      wacc <- weighted.acc(pred[-train], normed[-train,2])
      pacc <- sum(pred[-train]==normed[-train,2])/length(pred[-train])
      
      results <- rbind(results, c(h, tacc, wacc, pacc, c, 2))
      models[[(length(models)+1)]] <- ann
      
      # With neither range nor normalisation
      ann <- nnet(State ~ ., data=data[train,-1], weights=data.weights[train], size=h, decay=5e-4, maxit=200)
      pred <- predict(ann, data[,-1], type="class")
      tacc <- weighted.acc(pred[train], data[train,2])
      wacc <- weighted.acc(pred[-train], data[-train,2])
      pacc <- sum(pred[-train]==data[-train,2])/length(pred[-train])
      
      results <- rbind(results, c(h, tacc, wacc, pacc, c, 3))
      models[[(length(models)+1)]] <- ann
      
    }
  }
}

# Visualise results
cols <- do.call(rbind, Map(function(c)
{
  if (c==1) "green" 
  else if (c == 2) "blue" 
  else if (c == 3) "red" 
  else "black"
}, results[,5]))

pts <- do.call(rbind, Map(function(v)
{
  if (v==1) "r" # range
  else if (v == 2) "n" # normalised input
  else if (v == 3) "p" # 
  else "x"
}, results[,6]))


plot(results[,3] ~ results[,1], ylim=c(0,1), col=cols, pch=pts, xlab="Hidden neurons", ylab="Weighted accuracy")


# Save everything
save(results, file=paste0(basedir, "../ann.results.obj"))
save(models, file=paste0(basedir, "../ann.models.obj"))

write.table(results, file=paste0(basedir, "../ann.results.csv"), sep=",")

best.row <- match(max(results[,3]), results[,3])
best.ann <- models[[best.row]]
save(best.ann, file=paste0(basedir, "../best.ann.obj"))

finalPred <- predict(best.ann, normed[,-1], type="class")
table(normed[test,2], finalPred[test])
(288+1+12+49+1743+93+223+2)/length(test) #0.9301698 accuracy


data$ProjectedState <- finalPred
projected <- read.csv(file = "projected.csv")
projected$X <- NULL

ball1 <- subset(projected, projected$bearing == "b1")
ball2 <- subset(projected, projected$bearing == "b2")
ball3 <- subset(projected, projected$bearing == "b3")
ball4 <- subset(projected, projected$bearing == "b4")


vec <- vector()
for(rowNum in 11:nrow(ball1)){
  temp = data.frame("Time"= ball1$Time[c(1:rowNum)], "Freq"= ball1$HF.pow.y[c(1:rowNum)])
  x = xts(x=temp$Freq, order.by=as.POSIXct(temp$Time))
  x.ts = ts(x, freq=31536000, start=c(2017, 200))
  ARIMAfit <- auto.arima(x.ts, approximation=FALSE,trace=FALSE)
  pred <- forecast(ARIMAfit,10)
  vec[rowNum] <- pred[[4]][1]
  plot(pred,type="l",xlab = "Year",ylab = "Balance",lwd = 2,col = 'red',main="Forecasting using ARIMA Model")
}

vec[c(1:10)] = 0
ball1$ProjectedHFPow <- vec


vec <- vector()
for(rowNum in 11:nrow(ball2)){
  temp = data.frame("Time"= ball2$Time[c(1:rowNum)], "Freq"= ball2$HF.pow.y[c(1:rowNum)])
  x = xts(x=temp$Freq, order.by=as.POSIXct(temp$Time))
  x.ts = ts(x, freq=31536000, start=c(2017, 200))
  ARIMAfit <- auto.arima(x.ts, approximation=FALSE,trace=FALSE)
  pred <- forecast(ARIMAfit,10)
  vec[rowNum] <- pred[[4]][1]
  plot(pred,type="l",xlab = "Year",ylab = "Balance",lwd = 2,col = 'red',main="Forecasting using ARIMA Model")
}

vec[c(1:10)] = 0
ball2$ProjectedHFPow <- vec


vec <- vector()
for(rowNum in 11:nrow(ball3)){
  temp = data.frame("Time"= ball3$Time[c(1:rowNum)], "Freq"= ball3$HF.pow.y[c(1:rowNum)])
  x = xts(x=temp$Freq, order.by=as.POSIXct(temp$Time))
  x.ts = ts(x, freq=31536000, start=c(2017, 200))
  ARIMAfit <- auto.arima(x.ts, approximation=FALSE,trace=FALSE)
  pred <- forecast(ARIMAfit,10)
  vec[rowNum] <- pred[[4]][1]
  plot(pred,type="l",xlab = "Year",ylab = "Balance",lwd = 2,col = 'red',main="Forecasting using ARIMA Model")
}

vec[c(1:10)] = 0
ball3$ProjectedHFPow <- vec

vec <- vector()
for(rowNum in 11:nrow(ball4)){
  temp = data.frame("Time"= ball4$Time[c(1:rowNum)], "Freq"= ball4$HF.pow.y[c(1:rowNum)])
  x = xts(x=temp$Freq, order.by=as.POSIXct(temp$Time))
  x.ts = ts(x, freq=31536000, start=c(2017, 200))
  ARIMAfit <- auto.arima(x.ts, approximation=FALSE,trace=FALSE)
  pred <- forecast(ARIMAfit,10)
  vec[rowNum] <- pred[[4]][1]
  plot(pred,type="l",xlab = "Year",ylab = "Balance",lwd = 2,col = 'red',main="Forecasting using ARIMA Model")
}

vec[c(1:10)] = 0
ball4$ProjectedHFPow <- vec

write.csv(file= "ball1.csv",ball1)
write.csv(file= "ball2.csv",ball2)
write.csv(file= "ball3.csv",ball3)
write.csv(file= "ball4.csv",ball4)

ball11 <- read.csv("ball1.csv")
ball21 <- read.csv("ball2.csv")
ball31 <- read.csv("ball3.csv")
ball41 <- read.csv("ball4.csv")
ball11$X <- NULL
ball21$X <- NULL
ball31$X <- NULL
ball41$X <- NULL

finalData <- rbind(ball1,ball2,ball3,ball4)
write.csv(file= "finalData.csv",finalData)
