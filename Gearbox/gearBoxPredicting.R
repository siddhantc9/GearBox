

globVariables <- function(){
  
  basedir <<- "C:/Users/A672257/Siddhant/Workspace/Gearbox/";
  
}

x <- read.table(paste0(basedir, "models/bestColnames.txt"))

columnsToKeep <- function(data = NULL, KeepOnlyCols = NULL){
  filterData=data[,KeepOnlyCols]
  return(filterData)
}

xNew <- as.vector(x[-1,1])
y <- columnsToKeep(data = finalData, KeepOnlyCols = xNew)

x <- read.csv(paste0(basedir, "models/maxAbsColumns.csv"))

for(i in colnames(y)){
  
  maxAbs <- x[,i]
  y[,i] <- y[,i]/maxAbs
}

load(paste0(basedir, "models/best.ann.obj"))

predictions <- predict(best.ann, newdata = y, type="class")

x <- predictions[test]
y <- normed[test,2]
l <- union(x, y)
Table2 <- table(factor(x, l), factor(y, l))
print(confusionMatrix(Table2))

