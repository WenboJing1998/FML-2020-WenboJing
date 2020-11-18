
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






## Do experiments and record the results. Choose d=1,2,3,4. C=2^k (k=-13,...,7)
for(d in 1:4){
  for(k in 1:20){
    cat("C =",2^(k-14),"d =",d,'\n')
    ## SVM with certain d, C and 10-fold cross validation
    model <- svm(x = trainingset[,1:7], y = trainingset[,8], 
                 scale = TRUE, kernel =
                   "polynomial", degree = d, cost = 2^(k-14),
                 cross = 10)
    ## model[["accuracies"]] : CV errors
    CVerrors[d, k, ] <- c(1-mean(model[["accuracies"]])/100, 
                          sd(model[["accuracies"]])/100)
    ## model[["tot.nSV"]] : Number of Support Vectors. Compute the number
    ## of Marginal SVs by model[["tot.nSV"]] - sum_i(w_i==C)
    SVnumbers[d, k, ] <- c(model[["tot.nSV"]], 
                           model[["tot.nSV"]]-sum(abs(model[["coefs"]])==2^(k-14)))
    ## Compute test errors
    y_predict <- predict(model, testset[,1:7])
    accuarcy_test <- sum(y_predict == testset[,8])/ nrow(testset)
    TestErrors[d, k] <- 1-accuarcy_test
  }
}

################# plot requried figures #########################
par(mfrow=c(2,2))
x1 <- -13:6
y1 <- CVerrors[1, ,1]
upper1 <- y1 + CVerrors[1, ,2]
lower1 <- y1 - CVerrors[1, ,2]
plot(x1, y1, main='d=1', ylim = range(c(0.15, 0.55)),xlab=expression(paste(Log[2],'C')),ylab="CV error", type='o',col='red')
arrows(x1, upper1, x1, lower1, length =.05, angle = 90, code = 3, col='blue')

x1 <- -13:6
y1 <- CVerrors[2, ,1]
upper1 <- y1 + CVerrors[2, ,2]
lower1 <- y1 - CVerrors[2, ,2]
plot(x1, y1 , main='d=2', ylim = range(c(0.15, 0.55)),xlab=expression(paste(Log[2],'C')),ylab="CV error", type='o',col='red')
arrows(x1, upper1, x1, lower1, length =.05, angle = 90, code = 3, col='blue')

x1 <- -13:6
y1 <- CVerrors[3, ,1]
upper1 <- y1 + CVerrors[3, ,2]
lower1 <- y1 - CVerrors[3, ,2]
plot(x1, y1 , main='d=3', ylim = range(c(0.15, 0.55)),xlab=expression(paste(Log[2],'C')),ylab="CV error", type='o',col='red')
arrows(x1, upper1, x1, lower1, length =.05, angle = 90, code = 3, col='blue')

x1 <- -13:6
y1 <- CVerrors[4, ,1]
upper1 <- y1 + CVerrors[4, ,2]
lower1 <- y1 - CVerrors[4, ,2]
plot(x1, y1 , main='d=4', ylim = range(c(0.15, 0.55)),xlab=expression(paste(Log[2],'C')),ylab="CV error", type='o',col='red')
arrows(x1, upper1, x1, lower1, length =.05, angle = 90, code = 3, col='blue')

# find the minimal CV errors
which.min(CVerrors[,,1])

x = 1:4
y1 = CVerrors[ ,18,1]
y2 = TestErrors[, 18]
plot(x, y1 , ylim = range(c(0.15, 0.55)),xlab='d',ylab='', type='o',col='red',lty=2,pch=15)
lines(x, y2 ,ylim = range(c(0.15, 0.55)), type='o',col='blue', lty=1, pch=17)
legend("topleft", inset=.05, c("CV error","Test error"), lty=c(2, 1), pch=c(15, 17), col=c("red", "blue"))

x = 1:4
y3 = SVnumbers[ ,17,1]
y4 = SVnumbers[ ,17,2]
plot(x, y3 , ylim = range(c(0, 2500)),xlab='d',ylab='', type='o',col='green',lty=2,pch=18)
lines(x, y4 , type='o',col='purple', lty=1, pch=1)
legend("right", inset=.05, c("Total","Marginal"), lty=c(2, 1), pch=c(18, 1), col=c("green", "purple"))

##############Problem 6(d), similar as above, just changing the dataframe###################
# y = 2*as.numeric(data[ ,8])-3
CVerrors_new <- array(NA, dim=c(4,20,2))
TestErrors_new <- matrix(NA, nrow = 4, ncol = 20)

for(d in 1:4){
  
    data_new_temp <- as.matrix(data[,1:7]) %*% t(as.matrix(data[,1:7]))
    ## Transform data into a new form as derived in C.6(a), under a certain polynomial kernal with dimension d
     for (i in 1:ncol(data_new_temp)) {
      data_new_temp[,i] <- (data_new_temp[,i])^d * y[i]
    }
    data_new <- cbind(data_new_temp, as.factor(y))
    trainingset_new <- data_new[1:3133,]
    testset_new <- data_new[-(1:3133),]
    for(k in 1:20){
      k=1
       cat("C =",2^(k-14),"d =",d,'\n')
       model <- svm(x = trainingset_new[ ,1:4177], y = trainingset_new[ ,4178], 
                 scale = TRUE, kernel =
                   "polynomial", degree = 1, cost = 2^(k-14),
                 cross = 10)
    CVerrors_new[d, k, ] <- c(1-mean(model[["accuracies"]])/100, 
                          sd(model[["accuracies"]])/100)
 
    y_predict <- predict(model, testset_new[,1:7])
    accuarcy_test <- sum(y_predict == testset_new[,8])/ nrow(testset_new)
    TestErrors_new[d, k] <- 1-accuarcy_test
  }
}

################# plot requried figures #########################

# find the minimal CV errors
which.min(CVerrors[ , ,1])

x = 1:4
y1 = CVerrors_new[ ,18,1]
y2 = TestErrors_new[, 18]
plot(x, y1 , ylim = range(c(0.15, 0.55)),xlab='d',ylab='', type='o',col='red',lty=2,pch=15)
lines(x, y2 ,ylim = range(c(0.15, 0.55)), type='o',col='blue', lty=1, pch=17)
legend("topleft", inset=.05, c("CV error","Test error"), lty=c(2, 1), pch=c(15, 17), col=c("red", "blue"))

