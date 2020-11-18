
library(e1071)

## Load data, factorize the labels for SVM.##

data_original <- read.table('C://Users//WENBO JING//Desktop//NYU//Foundations of ML//HW2//abalone.data',sep=',')
data <- data_original[,2:9]
colnames(data) <- c('Length','Diameter','Height','Wweight','Sweight',
                    'Vweight','Sweight','Rings')
data[ ,8] <- as.factor(1*(data_original[,9]<=9))
trainingset <- data[1:3133, ]
testset <- data[-(1:3133), ]

# Create weak base learners. Here I use SVMs as base learners, each trained on a random sample set of size 10
base_learners <- list()
accuarcy_base <- rep(NA, 30) # record the test accuarcy of base learners to see how weak they are
for (j in 1:30) {
  subset_index <- sample(3133, 10,replace = F)
  h <- svm(x = trainingset[subset_index ,1:7], y = trainingset[subset_index,8], 
           scale = TRUE, kernel =
             "polynomial", degree = 1)
  base_learners <- c(base_learners, list(h))
  y_predict <- predict(h, testset[ ,1:7])
  accuarcy_base[j] <- sum(y_predict == testset[,8])/ nrow(testset)
}
summary(accuarcy_base)
sd(accuarcy_base)

#Adaboost  
#Input: base learners list, training set (x,y) , iteration number (iter), test set (xtest,ytest) 
#Output: A list of three arguments : vector alpha, vector k (the base learner chosen in each iteration), and test accuarcy

Adaboost <- function(base_learners, x, y, iter, xtest, ytest){
  m = nrow(x)
  D = rep(1/m, m)
  p = length(base_learners)
  alpha <- rep(0, iter)
  k <- rep(NA, iter)
  for (t in 1:iter) {
    
    # Compute all the \epsilon_{k,t} and choose the minimum
    epsilon <- rep(NA, p)
    for (j in 1:p) {
      y_predict <- predict(base_learners[[j]], x)
      epsilon[j] <- 1*(y_predict != y) %*% D
    } 
    epsilon_t <- min(epsilon)
    k[t] <- which.min(epsilon)
    
    # Update D
    h <- predict(base_learners[[k[t]]], x) # compute h_t(x_i)
    alpha[t] <- 0.5*(log((1-epsilon_t)/(epsilon_t)))
    Z_t <- 2*(epsilon_t*(1-epsilon_t))^{0.5}
    D <- (D * exp(-alpha[t]*(2*(y==h)-1)) )/Z_t
  }
 
  # test
  n <- nrow(xtest)
  y_pred <- rep(NA, n)
  for (i in 1:n) {
    # Compute \sum_{t=1}^T \alpha_t * h_t(x) on test set
    sum <- 0
    for (t in 1:iter) {
      h <- 2*(as.numeric(predict(base_learners[[k[t]]], xtest[i, ])))-3 #the predict function returns a factor, which becomes 1 and 2 when turned into numeric form. Time 2 and minus 3, then it becomes +1 and -1.
      sum <- sum + alpha[t] * h
    }
    y_pred[i] <- 1*(sum>0)
  }
  accuarcy_test <- sum(ytest == as.factor(y_pred)) / n
  return(list(alpha, k, accuarcy_test))
}


#10-fold Cross Validation to determine iteration number T 
accuarcy_CV <- matrix(NA, nrow=10, ncol=10)

for (j in 1:10) {
   for (i in seq(10)) {
     iter <- 100 * i
     ind <- sample(3133, 3133, replace=F)
  
     ind_train <- ind[(1+313*(j-1)):(1+313*j)]
     result <- Adaboost(base_learners, trainingset[ind_train, 1:7], trainingset[ind_train, 8], 
                       iter, trainingset[-ind_train, 1:7], trainingset[-ind_train, 8])
    
     accuarcy_CV[i,j] <- result[[3]]
  }
}

iter <- which.min(lapply(accuarcy_CV, sum))

Ada_result <- Adaboost(base_learners, trainingset[ ,1:7], trainingset[ ,8], 
                       iter, testset[ ,1:7], testset[ ,8]) 


#boost with logistic loss. Input: the same as above. Output: vector alpha and  test accuarcy
LogBoost <- function(base_learners, x, y, iter, xtest, ytest){
  p <- length(base_learners)
  alpha <- rep(0, p)
  m <- nrow(x)
  D <- rep(1/m, m)
  Z <- m / (2*log(2))
  DZ <- D * Z     # useful vector in computation
  for (t in 1:iter) {
    
    # Compute all the \epsilon_{k,t} and choose the minimum
    epsilon <- rep(NA, p)
    for (j in 1:p) {
      y_predict <- predict(base_learners[[j]], x)
      epsilon[j] <- 1*(y_predict != y) %*% D
    }
    epsilon_t <- min(epsilon)
    k_t <- which.min(epsilon)
    h <- predict(base_learners[[k_t]], x)
    
    # Update alpha, D, Z, and DZ
    eta <- 0.5*(log((1-epsilon_t)/(epsilon_t)))
    alpha[k_t] <- alpha[k_t] + eta
    DZ = 1/(log(2)+exp(eta*(2*(h==y)-1))*(1/(DZ)-log(2))) # this formula can be derived by definition
    Z <- sum(DZ) 
    D <- DZ/Z
  }
  
   # test
  n <- nrow(xtest)
  y_pred <- rep(NA, n)
  for (i in 1:n) {
    sum <- 0
    for (j in 1:p) {
      h <- 2*(as.numeric(predict(base_learners[[j]], xtest[i, ])))-3
      sum <- sum + alpha[j] * h
    }
    y_pred[i] <- 1*(sum>0)
  }
  accuarcy_test <- sum(ytest == as.factor(y_pred)) / n
  return(list(alpha, accuarcy_test))
}

#10-fold Cross Validation to determine iteration number T 
accuarcy_CV_Log <- matrix(NA, nrow=10, ncol=10)

for (j in 1:10) {
  for (i in seq(10)) {
    iter <- 100 * i
    ind <- sample(3133, 3133, replace=F)
    
    ind_train <- ind[(1+313*(j-1)):(1+313*j)]
    result <- LogBoost(base_learners, trainingset[ind_train, 1:7], trainingset[ind_train, 8], 
                       iter, trainingset[-ind_train, 1:7], trainingset[-ind_train, 8])
    
    accuarcy_CV_Log [i,j] <- result[[3]]
  }
}

iter <- which.min(lapply(accuarcy_CV_Log , sum))
Log_result <- LogBoost(base_learners, trainingset[ ,1:7], trainingset[ ,8], 
                       iter, testset[ ,1:7], testset[ ,8]) 
