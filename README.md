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
# An example
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
