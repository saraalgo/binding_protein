prob_1 = function (x) {
  cross_1 = c("A", "C", "G", "T")
  prob = summary(factor(strsplit(x, split = "")[[1]], levels = cross_1), 
                maxsum = 5)/nchar(x)
  return(prob)
}


prob_2 = function (x) {
  
  cross_1 = c("A", "C", "G", "T")
  cross_2 = as.vector((outer(cross_1, cross_1, paste, sep = "")))
  xSplitted = strsplit(x, split = "")[[1]]
  n = nchar(x)
  prob = summary(factor(paste0(xSplitted[-n], xSplitted[-1]), 
                        levels = cross_2), maxsum = 17)/(n)
  return(prob)
}


prob_3 = function (x) {
  
  cross_1 = c("A", "C", "G", "T")
  cross_2 = as.vector((outer(cross_1, cross_1, paste, sep = "")))
  cross_3 = as.vector((outer(cross_2, cross_1, paste, sep = "")))
  xSplitted = strsplit(x, split = "")[[1]]
  n = nchar(x)
  prob = summary(factor(paste0(paste(xSplitted[-c(n, n - 1)], 
                        xSplitted[-c(1, n)]), xSplitted[-c(1, 2)]), 
                        levels = cross_3), maxsum = 65)/(n)
  return(prob)
}


prob_4 = function (x) {
  
  cross_1 = c("A", "C", "G", "T")
  cross_2 = as.vector((outer(cross_1, cross_1, paste, sep = "")))
  cross_3 = as.vector((outer(cross_2, cross_1, paste, sep = "")))
  cross_4 = as.vector((outer(cross_3, cross_1, paste, sep = "")))
  xSplitted = strsplit(x, split = "")[[1]]
  n = nchar(x)
  prob = summary(factor(paste0(paste0(paste0(xSplitted[-c(n, n - 1, n - 2)], 
                                    xSplitted[-c(1, n)]), xSplitted[-c(1, 2)]),
                                    xSplitted[-c(1,2,3)]), 
                                    levels = cross_4), maxsum = 257)/(n)
  return(prob)
}

prob_5 = function (x) {
  
  cross_1 = c("A", "C", "G", "T")
  cross_2 = as.vector((outer(cross_1, cross_1, paste, sep = "")))
  cross_3 = as.vector((outer(cross_2, cross_1, paste, sep = "")))
  cross_4 = as.vector((outer(cross_3, cross_1, paste, sep = "")))
  cross_5 = as.vector((outer(cross_4, cross_1, paste, sep = "")))
  xSplitted = strsplit(x, split = "")[[1]]
  n = nchar(x)
  prob = summary(factor(paste0(paste0(paste0(paste0(xSplitted[-c(n, n - 1, n - 2)], 
                              xSplitted[-c(1, n)]), xSplitted[-c(1, 2)]),
                               xSplitted[-c(1,2,3)]), xSplitted[-c(1,2,3,4)]), 
                        levels = cross_5), maxsum = 1025)/(n)
  return(prob)
}

prob_6 = function (x) {
  
  cross_1 = c("A", "C", "G", "T")
  cross_2 = as.vector((outer(cross_1, cross_1, paste, sep = "")))
  cross_3 = as.vector((outer(cross_2, cross_1, paste, sep = "")))
  cross_4 = as.vector((outer(cross_3, cross_1, paste, sep = "")))
  cross_5 = as.vector((outer(cross_4, cross_1, paste, sep = "")))
  cross_6 = as.vector((outer(cross_5, cross_1, paste, sep = "")))
  xSplitted = strsplit(x, split = "")[[1]]
  n = nchar(x)
  prob = summary(factor(paste0(paste0(paste0(paste0(paste0(xSplitted[-c(n, n - 1, n - 2)], 
                          xSplitted[-c(1, n)]), xSplitted[-c(1, 2)]), xSplitted[-c(1,2,3)]), 
                          xSplitted[-c(1,2,3,4)]),xSplitted[-c(1,2,3,4,5)]), 
                        levels = cross_6), maxsum = 4097)/(n)
  return(prob)
}

prob_7 = function (x) {
  
  cross_1 = c("A", "C", "G", "T")
  cross_2 = as.vector((outer(cross_1, cross_1, paste, sep = "")))
  cross_3 = as.vector((outer(cross_2, cross_1, paste, sep = "")))
  cross_4 = as.vector((outer(cross_3, cross_1, paste, sep = "")))
  cross_5 = as.vector((outer(cross_4, cross_1, paste, sep = "")))
  cross_6 = as.vector((outer(cross_5, cross_1, paste, sep = "")))
  cross_7 = as.vector((outer(cross_6, cross_1, paste, sep = "")))
  xSplitted = strsplit(x, split = "")[[1]]
  n = nchar(x)
  prob = summary(factor(paste0(paste0(paste0(paste0(paste0(paste0(xSplitted[-c(n, n - 1, n - 2)], 
                                      xSplitted[-c(1, n)]), xSplitted[-c(1, 2)]), xSplitted[-c(1,2,3)]), 
                                      xSplitted[-c(1,2,3,4)]),xSplitted[-c(1,2,3,4,5)]),xSplitted[-c(1,2,3,4,5,6)]), 
                                      levels = cross_7), maxsum = 16385)/(n)
  return(prob)
}


prob_8 = function (x) {
  
  cross_1 = c("A", "C", "G", "T")
  cross_2 = as.vector((outer(cross_1, cross_1, paste, sep = "")))
  cross_3 = as.vector((outer(cross_2, cross_1, paste, sep = "")))
  cross_4 = as.vector((outer(cross_3, cross_1, paste, sep = "")))
  cross_5 = as.vector((outer(cross_4, cross_1, paste, sep = "")))
  cross_6 = as.vector((outer(cross_5, cross_1, paste, sep = "")))
  cross_7 = as.vector((outer(cross_6, cross_1, paste, sep = "")))
  cross_8 = as.vector((outer(cross_7, cross_1, paste, sep = "")))
  xSplitted = strsplit(x, split = "")[[1]]
  n = nchar(x)
  prob = summary(factor(paste0(paste0(paste0(paste0(paste0(paste0(paste0(xSplitted[-c(n, n - 1, n - 2)],
                        xSplitted[-c(1, n)]), xSplitted[-c(1, 2)]), xSplitted[-c(1,2,3)]), 
                        xSplitted[-c(1,2,3,4)]), xSplitted[-c(1,2,3,4,5)]),xSplitted[-c(1,2,3,4,5,6)]), 
                        xSplitted[-c(1,2,3,4,5,6,7)]), levels = cross_8), maxsum = 65537)/(n)
  return(prob)
}


