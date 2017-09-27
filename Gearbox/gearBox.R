
finalData$timestamp <- as.factor(finalData$timestamp)
finalData$bearing <- as.factor(finalData$bearing)

b1 <- subset(finalData, bearing == "b1")
b2 <- subset(finalData, bearing == "b2")
b3 <- subset(finalData, bearing == "b3")
b4 <- subset(finalData, bearing == "b4")


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
  cbind((b1[(b1$State == "normal"),-c(47, 48, 49)]), bearing="b1"),
  cbind((b2[(b2$State == "normal"),-c(47, 48, 49)]), bearing="b2"),
  cbind((b3[(b3$State == "normal"),-c(47, 48, 49)]), bearing="b3"),
  cbind((b4[(b4$State == "normal"),-c(47, 48, 49)]), bearing="b4")
)
library(corrgram)
#corrgram(norm)

# Find high (absolute) correlations (ignoring non-numeric bearing name col)
cor <- cov2cor(cov(norm[, c(-47)]))
alikes <- apply(cor, 2, function(col) { names(Filter(function (val) { val > 0.9 }, sort.int(abs(col), decreasing=TRUE))) } )

cat(str(alikes, vec.len=10))

# As an offline exercise, minimise the set of features by removing those with high correlations

# I found offline analysis to give the following minimal set
best <- c("min_x","median_x","max_x","mean_x","skew_x","kurt_x",  
          "ftf_x","bpfi_x","bpfo_x","bsf_x","f2_x","f3_x",
          "f4_x", "f5_x", "min_y","max_y","skew_y","kurt_y",
          "ftf_y","bpfi_y","bpfo_y","bsf_y","f2_y","f3_y",
          "f4_y", "f5_y","qu_1_x","vhf_pow_x","qu_1_y","median_y","hf_pow_y")

# Select all rows from all bearings, only the best feature cols plus State, and bind on bearing
data <- rbind(
  cbind(bearing="b1", (b1[,c("State", best)])),
  cbind(bearing="b2", (b2[,c("State", best)])),
  cbind(bearing="b3", (b3[,c("State", best)])),
  cbind(bearing="b4", (b4[,c("State", best)]))
)

#write.table(data, file=paste0(basedir, "../all_bearings.csv"), sep=",", row.names=FALSE)


library("entropy")
# MI = H(x) + H(y) - H(x, y)
H.x <- entropy(table(data$State))
mi <- apply(data[, -c(1,2)], 2, function(col) { H.x + entropy(table(col)) - entropy(table(data$State, col))})
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
  cm <- confusionMatrix(pred[-train.rows], data[-train.rows, "State"]) # only the non-training rows
  # pull out the answer
  accuracies <- c(accuracies, cm["overall"][[1]]["Accuracy"][[1]])
}

plot(accuracies, xlab="Number of features", ylab="Accuracy")

data <- data[c("bearing", "State", sorted[1:14])]

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
save(best.kmeans, file=paste0(basedir, "models/kmeans.obj"))


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
print(table(data.frame(data$State, best.kmeans$cluster)))

# Read in the best kmeans model. Straight after running kmeans.R the 
# filename will be "kmeans.obj", but this is the best model I found.
load(paste0(basedir, "models/kmeans.obj"))

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
cw1 <- rep(1, 9) # all equal
cw2 <- c(10, 100, 100, 10, 1, 10, 1,10,1) # 1/order of count

freqs <- as.data.frame(table(data$State))
cw3 <- cbind(freqs[1], apply(freqs, 1, function(s) { length(data[,1])/as.integer(s[2])})) # 1/weight

class.weights <- rbind(cw1, cw2, cw3[,2])
colnames(class.weights) <- c("early", "failure.b2", "failure.inner", "failure.roller", "normal", "stage2", "suspect","unknown","failure.b1")

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
             