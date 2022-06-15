
# install.packages("'mvtnorm")


require(mvtnorm)
require(MASS)

for(global_iter in 1:30){
  
  d = 5
  K = 100
  # N <- c(rep(20,1), rep(200,1))
  N <- sample(c(190:200), size = K, replace = TRUE)
  N_test <- 1000
  
  
  A <- matrix(runif(K^2)*2-1, ncol=K) 
  Omega <- (t(A)%*%A)
  # Omega <- Omega/sum(diag(Omega))
  mu <- rep(0, K)
  Theta <- mvrnorm(d,mu,Omega)
  
  X <- list()
  Y <- list()
  X_test <- list()
  Y_test <- list()
  
  for(i in 1:K){
    
    X1 <- matrix(rnorm(d*N[i]), nrow=d,ncol=N[i])
    Y1 <- t(X1)%*%Theta[,i] + mvrnorm(1, rep(0, N[i]), diag(0.05,N[i]))
    X <- append(X, list(X1))
    Y <- append(Y, list(Y1))
    
    X2 <- matrix(rnorm(d*N_test), nrow=d,ncol=N_test)
    Y2 <- t(X2)%*%Theta[,i] # + mvrnorm(1, rep(0, N_test), diag(0.05,N_test))
    X_test <- append(X_test, list(X2))
    Y_test <- append(Y_test, list(Y2))
    
  }
  
  local_SGD <- function(theta, X, Y, LR, iter, corr, use=1){
    
    for(i in 1:iter){
      theta <- theta + 2*LR*X%*%(Y-t(X)%*%theta)
    }
    
    if(use==1) theta <- theta - 2*LR*corr
    
    return(theta)
    
  }
  
  Theta0 <- matrix(rep(0, d*K), nrow=d, ncol=K)
  Theta0_no_comm <- matrix(rep(0, d*K), nrow=d, ncol=K)
  Omega0 <- diag(1, K)
  
  LR <- 0.0001
  iter <- 20
  
  # theta01_no_corr <- matrix(c(0,0),2,1)
  # theta02_no_corr <- matrix(c(0,0),2,1)
  
  acc <- NULL
  
  for(comm in 1:100){
    
    communication_corr <- matrix( rep(0, d*K), nrow=d, ncol=K )
    Omega0_inv <- solve(Omega)
    
    for(k in 1:K){
      for(i in 1:K){
        communication_corr[,k] <- communication_corr[,k] + Theta0[,i]*Omega0_inv[i,k] 
      }
    }
    
    for(k in 1:K){
      
      Theta0[,k] <- local_SGD(Theta0[,k], X[[k]], Y[[k]], LR, iter, communication_corr[,k])
      
    }
    
    for(k in 1:K){
      
      Theta0_no_comm[,k] <- local_SGD(Theta0_no_comm[,k], X[[k]], Y[[k]], LR, iter, communication_corr[,k], 0)
      
    }
    
    
    # corr <- theta01*solve(Omega)[1,1] + theta02*solve(Omega)[2,1]
    # theta01 <- local_SGD(theta01, X1, Y1, LR, iter, corr)
    # theta01_no_corr <- local_SGD(theta01_no_corr, X1, Y1, LR, iter, corr, 0)
    # 
    # corr <- theta02*solve(Omega)[1,2] + theta02*solve(Omega)[2,2]
    # theta02 <- local_SGD(theta02, X2, Y2, LR, iter, corr)
    # theta02_no_corr <- local_SGD(theta02_no_corr, X2, Y2, LR, iter, corr, 0)
    # 
    # Theta_new <- matrix(cbind(theta01, theta02), 2, 2)
    # Theta_new_no_corr <- matrix(cbind(theta01_no_corr, theta02_no_corr), 2, 2)
    
    new_corr <- t(Theta0)%*%Theta0
    # eigen_decomposition <- eigen(new_corr)
    # e_vec <- eigen_decomposition$vectors
    # new_corr <- e_vec %*%diag(sqrt(round(eigen_decomposition$values, 5)),K)%*%solve(e_vec )
    # Omega0 <- 0.9*Omega0+0.1*( new_corr/ sum(diag(new_corr)) )
    
    Omega0 <- 0.9*Omega0+0.1/d*new_corr#*( sqrt(new_corr)/ sum(diag(sqrt(new_corr))) )
    
    acc <- c(acc, norm(Theta0-Theta)/sqrt(K))
    
  }
  
  #Omega0
  #norm(Omega0/sum(diag(Omega0)))
  #norm(Omega/sum(diag(Omega)))
  
  aa1 = t(X_test[[1]])%*%Theta0[,1]; RMSE1 <- sqrt(sum((aa1-Y_test[[1]])^2)/1000) # norm(aa1-Y_test[[1]])/sqrt(N_test)
  aa2 = t(X_test[[1]])%*%Theta0_no_comm[,1]; RMSE2 <- norm(aa2-Y_test[[1]], "2")/sqrt(N_test)
  c(RMSE1,RMSE2)
  # 
  # plot(aa1, type="l")
  # lines(Y[[1]], col="red")
  
  
  if(global_iter == 1 & acc[100]<=1) plot(acc, type="l")
  if(global_iter == 1 & acc[100]>1) global_iter = 1
  if(global_iter >1 & acc[100]<=1) lines(acc, col=global_iter)
}

#Theta0[,1:6]
#Theta[,1:6]
