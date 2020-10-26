## R package e1071 contains libsvm
library(e1071)

## Load data, factorize the labels for SVM. The scaling process will be
## automatically done in the SVM function, by setting scale=TRUE.

data_original <- read.table('C://Users//WENBO JING//Desktop//NYU//Foundations of ML//HW2//abalone.data',sep=',')
data <- data_original[,2:9]
colnames(data) <- c('Length','Diameter','Height','Wweight','Sweight',
                    'Vweight','Sweight','Rings')
data[ ,8] <- as.factor(1*(data_original[,9]<=9))
trainingset <- data[1:3133, ]
testset <- data[-(1:3133), ]

## The result matrix. Choose d=1,2,3,4. C=2^k (k=-13,...,7)
CVerrors <- array(NA, dim=c(4,20,2)) # means and standard errors of CV errors 
TestErrors <- matrix(NA, nrow = 4, ncol = 20) # test errors
SVnumbers <- array(NA, dim=c(4,20,2)) # numbers of SV and marginal SV.

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

