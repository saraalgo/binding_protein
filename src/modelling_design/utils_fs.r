# Preprocessing
# ===================================
## FCBF (Fast Correlation Based Filter for Feature Selection)
fast.cor.FS = function(data, thresh, normalize = T){

  stopifnot('target' %in% names(data))
  require(FCBF)

  y = as.factor(data$target)
  x = subset(data, select = - c(target))

  if (normalize == T) {
    x = apply(x, 2, function(x) log2(x + 1))
  }

  dis = discretize_exprs(t(x))

  fcbf = fcbf(dis, y, verbose = T, thresh = thresh)
  xx = x[,fcbf$index]

  xx = as.data.frame(cbind(xx, target = y))

  return(xx)
}

## Kruskal
kruskal = function(data, fs.type = 'kruskal.test', perc, output.dir, filename){

  require(mlr)
  stopifnot('target' %in% names(data))

  task = makeClassifTask(data = data, target = 'target')

  tasks = filterFeatures(task, method = fs.type, perc = perc)

  tasks$task.desc$id = paste(fs.type, ncol(tasks$env$data) -1, sep = "_")

  tasks = tasks$env$data

  return(tasks)
}

## PCA
pca = function(data){

  y = as.factor(data$target)
  x = subset(data, select = - c(target))

  res = prcomp(as.matrix(x), scale = FALSE)

  res = cbind.data.frame(res$x[,1:15], target = y)

  return(res)
}






