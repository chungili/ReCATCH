# Introduction
ReCATCH package is applied to high-dimensional tensor classification modeling with heterogeneous covariates.
# Install the package from GitHub
1. Install the devtools package
```{R}
install.packages("devtools")
````
2. Load the devtools package
```{R}
library(devtools)
````
3. Install the package
```{R}
install_github("chungili/ReCATCH")
````
# A simulated data
```{R}
# Generate a training dataset
set.seed(123)
Sim_test<-GenCovTensor(q = 110, ng=c(500,500), pk = c(30,36,30), beta2 = 0.8, phi2 = 0.2, alpha2 = 0.5)
Sim_train<-GenCovTensor(q = 110, ng=c(50,50), pk = c(30,36,30), beta2 = 0.8, phi2 = 0.2, alpha2 = 0.5)
x <- Sim_train$x
U <- Sim_train$U
y <- Sim_train$y
testx <- Sim_test$x
testU <- Sim_test$U
testy <- Sim_test$y
# fit the model
out<-ReCatch(x=x, U=U, y=y, testx=testx, testU=testU, setseed=120, imbalanced=T, OptMethod = "kf")
# predict categorical responses
out$pred
````
# Real Dataset: Breast Cancer
```{R}
data("breast")
help(breast)
# test for homogeneity of covariance matrices 
library(heplots)
res <- boxM(trainU, trainOP)
# fit CATCH model
library(catch)
set.seed(1234)
catch.lambda <- cv.catch(x = trainX, z = trainU, y = trainOP, nfolds = 4)$lambda.min
fit.catch <- catch(x=trainX, z = trainU, y = trainOP,
             testx = testX, testz = testU, lambda = catch.lambda)
library(caret)
CM.catch <- confusionMatrix(factor(fit.catch$pred), factor(testOP), positive = "2")
round(c(CM.catch$overall["Accuracy"], CM.catch$byClass[c("Recall", "Precision", "F1")]), 2)
# fit CATCH+ model
fit.CatchPlus <- CatchPlus(x = trainX, U = trainU, y = trainOP,
                     testx = testX, testU = testU, testy = testOP,
                     setseed = 1234, imbalanced = FALSE, n.folds = 4)
CM.CatchPlus <- confusionMatrix(factor(fit.CatchPlus$pred), factor(testOP), positive = "2")
round(c(CM.CatchPlus$overall["Accuracy"], CM.CatchPlus$byClass[c("Recall", "Precision", "F1")]), 2)
```
# Real Dataset: Colorectal Cancer

```{R}
data("CRC")
help(CRC)
# Tests for Homogeneity of Covariance Matrices
library(covTestR)
Srivastava2014(list(DE[recurrence==1,], DE[recurrence==2,]))
# Calculate the performance metrics based on random splits 
library(caret)
library(MASS)
# To decrease computational expenses, 
# the number of iterations is set to 10
B <- 10
res_ReCATCH <- matrix(NA, B, 4)
res_CATCH <- res_ReCATCH
N <- length(recurrence)
for (j in 1:B){
  cat("Iteration:", j, "\n")
  set.seed(1234+j)
  ind_train <- createDataPartition(as.factor(recurrence), p = 0.7)$Resample1
  ind_test <- setdiff(1:N, ind_train)
  ###
  n <- length(ind_train)
  d <- N-n
  train_U <- DE[ind_train, ]
  test_U <- DE[-ind_train, ]
  train_X <- array(list(), n)
  for(w in 1:n) train_X[[w]] <- RF[[ ind_train[w]  ]]
  test_X <- array(list(), d)
  for(w in 1:d) test_X[[w]] <- RF[[ ind_test[w]  ]]
  train_y <- recurrence[ind_train]
  test_y <- recurrence[-ind_train]
  ### CATCH
  cat("Fit CATCH model: ", "\n")
  fit.catch <- catch(x = train_X, z = train_U, y = train_y, 
                     testx = train_X, testz = train_U )
  predva <- fit.catch$pred
  refe<-ifelse(train_y==1, "No", "Yes")
  Fsc <- rep(0, length(fit.catch$lambda))
  for (i in 1:length(fit.catch$lambda)){
    pred_vai <- ifelse(predva[,i]==1, "No", "Yes")
    CM <- confusionMatrix(factor(pred_vai), factor(refe), 
                          positive = "Yes", mode = "prec_recall")
    Fsc[i] <- CM$byClass["F1"]
  }
  lamda <- fit.catch$lambda[which.max(Fsc)]
  fit.catch <- catch(x = train_X, z = train_U, y = train_y, 
                     testx = test_X, testz = test_U, lambda = lamda )
  pred.catch <- as.vector(fit.catch$pred)
  out.catch <- confusionMatrix(factor(fit.catch$pred), factor(test_y), positive = "2")
  res_CATCH[j,] <- c(out.catch$overall["Accuracy"],
                     out.catch$byClass[c("Recall", "Precision", "F1")])
  ### ReCATCH
  cat("Fit ReCATCH model: ", "\n")
  yhat.Re <- ReCatch(x = train_X, U=train_U, y=train_y,
                     testx = test_X, testU = test_U, OptMethod = "loo",
                     imbalanced=TRUE, setseed=(1234+j) )$pred
  CM.Re <- confusionMatrix(data = as.factor(yhat.Re),
                           reference = as.factor(test_y), positive = "2")
  res_ReCATCH[j,] <- c(CM.Re$overall["Accuracy"],
                       CM.Re$byClass[c("Recall", "Precision", "F1")])
}
colnames(res_ReCATCH) <- c("Accuracy", "Recall", "Precision", "F1")
colnames(res_CATCH) <- c("Accuracy", "Recall", "Precision", "F1")
# performance metrics
round(rbind(colMeans(res_CATCH, na.rm = TRUE), 
            colMeans(res_ReCATCH, na.rm = TRUE)), 2)
```
