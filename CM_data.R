
setwd("CMAPSSData")
FD001 <- read.table("train_FD001.txt", head=FALSE)
#View(FD001)


# Among 6-26, remove V6, 10, 11, 15, 21, 22, 23, 24
# 6 7 8 9 10 11 12 13 14 15
# 1 2 3 4 5 6 7 8 9 10
# plot(FD001[1:192,8])


# Create Training Data
x <- list()
x_test <- list()
y <- list()
y_test <- list()

sensor_index <- 7

range01 <- function(x, min_x, max_x){
  
  x_new <- (x-min_x)/(max_x-min_x)
  return(x_new)
  
  # x_new <- (x-min(x))/(max(x)-min(x))
  # return(list(x_new, min(x), max(x)))

}

RMSE_function <- function(pred, test, N_test){
  
  RMSE <- norm(pred-test, "2")/sqrt(N_test)
  return(RMSE)
  
}

local_SGD <- function(theta, X, Y, LR, iter, corr, use=1){
 
  # SGD algorithm 
  # theta <- matrix(theta, ncol = 1)
  
  for(i in 1:iter){
    theta <- theta + 2*LR*X%*%(Y-t(X)%*%theta)
  }
  
  if(use==1) theta <- theta - 2*LR*corr
  
  return(theta)
  
}

d = 4
percentage = 0.7

for(i in 1:100){
  
  Data <- FD001[(FD001[, 1]==i),]
  data_size <- dim(Data)[1]
  train_size <- round(dim(Data)[1]*0.6)
  test_size <- data_size - train_size

  train_data <- Data[1:train_size,]
  test_data <- Data[(train_size+1):data_size,]
  
  train_data[,2] <- range01(train_data[,2], 1, 300)
  temp_matrix <- matrix(0, nrow=train_size, ncol=d+1)
  for(j in 1:(d+1)) temp_matrix[,j] = train_data[,2]^(j-1)
  temp_matrix <- t(temp_matrix)
  
  x <- append(x, list(temp_matrix))
  
  temp_y <- range01(train_data[,sensor_index], min(train_data[,sensor_index]), max(train_data[,sensor_index]))
  y <- append(y, list(temp_y))
  
  
  ###  test
  test_data[,2] <- range01(test_data[,2], 1, 300)
  temp_matrix <- matrix(0, nrow=test_size, ncol=d+1)
  for(j in 1:(d+1)) temp_matrix[,j] = test_data[,2]^(j-1)
  temp_matrix <- t(temp_matrix)
  
  x_test <- append(x_test, list(temp_matrix))
  
  # x_test <- append(x_test, list(test_data[,2]))
  temp_y <- range01(test_data[,sensor_index], min(test_data[,sensor_index]), max(test_data[,sensor_index]))
  y_test <- append(y_test, list(temp_y))
  # y_test <- append(y_test, list(test_data[,sensor_index]))
  
  # x <- append(x, list(range01(Data[,2])))
  # y_local <- Data[,7]
  # mean_y = mean(y_local)
  # std_y = sqrt(var(y_local))
  # y_local <- (y_local-mean_y)/std_y
  # y <- append(y, list(y_local))
  
}

K = 100
LR <- 0.001
iter <- 20

Theta0 <- matrix(rep(0, (d+1)*K), nrow=d+1, ncol=K)
Theta0_no_comm <- matrix(rep(0, (d+1)*K), nrow=d+1, ncol=K)
Omega0 <- diag(1, K)


for(comm in 1:100){
  
  communication_corr <- matrix( rep(0, (d+1)*K), nrow=d+1, ncol=K )
  Omega0_inv <- solve(Omega0)
  
  for(k in 1:K){
    for(i in 1:K){
      communication_corr[,k] <- communication_corr[,k] + Theta0[,i]*Omega0_inv[i,k] 
    }
  }
  
  for(k in 1:K){
    
    Theta0[,k] <- local_SGD(Theta0[,k], x[[k]], y[[k]], LR, iter, communication_corr[,k])
    
  }
  
  for(k in 1:K){
    
    Theta0_no_comm[,k] <- local_SGD(Theta0_no_comm[,k], x[[k]], y[[k]], LR, iter, communication_corr[,k], 0)
    
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
  
  Omega0 <- 0.8*Omega0+0.8/d*new_corr#*( sqrt(new_corr)/ sum(diag(sqrt(new_corr))) )
  
  # acc <- c(acc, norm(Theta0-Theta)/sqrt(K))
  
}

# plot(y[[1]])
# lines(t(x[[1]])%*%Theta0[,1])
# 
# plot(y_test[[1]]~x_test[[1]][2,])
# lines(t(x_test[[1]])%*%Theta0[,1]~x_test[[1]][2,])
# lines(t(x_test[[1]])%*%Theta0_no_comm[,1]~x_test[[1]][2,])
# 
# 
# 
# xxx <- cbind(x[[1]], x_test[[1]])
# yyy <- c(y[[1]], y_test[[1]])
# plot(yyy)
# lines(t(xxx)%*%Theta0[,1])
# lines(t(xxx)%*%Theta0_no_comm[,1], col="red")
# 
# 
# 
# 
# RMSE_function(t(xxx)%*%Theta0[,1], yyy, length(yyy))
# RMSE_function(t(xxx)%*%Theta0_no_comm[,1], yyy, length(yyy))
# 
# 
# RMSE_function(t(x_test[[1]])%*%Theta0[,1], y_test[[1]], length(y_test[[1]]))
# RMSE_function(t(x_test[[1]])%*%Theta0_no_comm[,1], y_test[[1]], length(y_test[[1]]))

