sim_cure_GOR <- function(N, cov_theta_0,cov_beta_0,cov_theta_xi,cov_beta_xi,theta_0 = c(-1,1,0),theta_xi=c(1,1),
                        beta_0 = c(0,0.5),beta_xi=c(1,1),A = 5, B = 15) {
  u <- stats::runif(N)
  a <- stats::runif(N)
  C <- cbind(A,a * B)
  C <- C[,1] * (C[,1] <= C[,2]) + C[,2] * (C[,1] > C[,2])
  theta = c(theta_0,theta_xi)
  beta = c(beta_0,beta_xi)
  cov_theta= cbind(cov_theta_0,cov_theta_xi)
  colnames(cov_theta)=c("i",paste0("theta_z",seq(1,dim(cov_theta_0)[2]-1)),paste0("theta_xi",seq(1,dim(cov_theta_xi)[2])))
  
  
  cov_beta= cbind(cov_beta_0,cov_beta_xi)
  colnames(cov_beta)=c(paste0("beta_z",seq(1,dim(cov_beta_0)[2])),paste0("beta_xi",seq(1,dim(cov_beta_xi)[2])))
  eta <- exp(as.vector(theta %*% t(cov_theta)))
  
  cure = 0*eta
  beta_x <- as.vector(beta %*% t(cov_beta))
  exp_pred_beta <- exp(beta_x)
  
  tempos = C
  for(i in 1:N){
  piz = eta[i]/(eta[i]+1)
  cure[i] = rbinom(1,1,piz)
  tempos[i] <- ifelse(cure[i] != 0, rexp(1,exp_pred_beta[i]), Inf)
  }
  
  
  
  
  
  #Z <- ifelse(tempos < C, tempos, C)
  Z <- ifelse(tempos < C, tempos, C)
  delta <- ifelse(tempos < C, 1, 0)
  L <- R <- Z * NA
  for (i in 1:N) {
    if(delta[i] == 0) {
      L[i] <- Z[i]
      R[i] <- Inf
    }
    else {
      L[i] <- 0
      add <- stats::runif(1, 0.1, 0.5)
      R[i] <- add
      check <- (L[i] <= Z[i] & Z[i] < R[i])
      while(!check) {
        L[i] <- L[i] + add
        add <- stats::runif(1, 0.1, 0.5)
        R[i] <- R[i] + add
        check <- (L[i] <= Z[i] & Z[i] < R[i])
      }
    }
  }
  dados <- data.frame(Z, L, R, delta, cov_theta[,2:dim(cov_theta)[2]], cov_beta)
  return(dados)
}






sim_cure <- function(N, cov_theta_0,cov_beta_0,cov_theta_xi,cov_beta_xi,theta_0 = c(-1,1,0),theta_xi=c(1,1),
                        beta_0 = c(0,0.5),beta_xi=c(1,1),A = 5, B = 15) {
  u <- stats::runif(N)
  a <- stats::rexp(N)
  C <- cbind(A,a * B)
  C <- C[,1] * (C[,1] <= C[,2]) + C[,2] * (C[,1] > C[,2])
  theta = c(theta_0,theta_xi)
  beta = c(beta_0,beta_xi)
  cov_theta= cbind(cov_theta_0,cov_theta_xi)
  colnames(cov_theta)=c("i",paste0("theta_z",seq(1,dim(cov_theta_0)[2]-1)),paste0("theta_xi",seq(1,dim(cov_theta_xi)[2])))
  
  
  cov_beta= cbind(cov_beta_0,cov_beta_xi)
  colnames(cov_beta)=c(paste0("beta_z",seq(1,dim(cov_beta_0)[2])),paste0("beta_xi",seq(1,dim(cov_beta_xi)[2])))
  eta <- exp(as.vector(theta %*% t(cov_theta)))
  K_vector <- stats::rpois(N, eta / 2)
  U_vector <- K_vector * NA
  for(i in 1:length(K_vector)) {
    if(K_vector[i] == 0) U_vector[i] <- 0
    else{
      U_vector[i] <- 0
      for (j in 1:K_vector[i])
        U_vector[i] <- U_vector[i] + stats::rchisq(1, 2, ncp = 0)
    }
  }
  beta_x <- as.vector(beta %*% t(cov_beta))
  exp_pred_beta <- exp(beta_x)
  num <- -2 * log(1 - u)
  den <- U_vector * exp_pred_beta
  tempos <- ifelse(U_vector != 0, (num / den), Inf)
  Z <- ifelse(tempos < C, tempos, C)
  delta <- ifelse(tempos < C, 1, 0)
  L <- R <- Z * NA
  for (i in 1:N) {
    if(delta[i] == 0) {
      L[i] <- Z[i]
      R[i] <- Inf
    }
    else {
      L[i] <- 0
      add <- stats::runif(1, 0.1, 0.5)
      R[i] <- add
      check <- (L[i] <= Z[i] & Z[i] < R[i])
      while(!check) {
        L[i] <- L[i] + add
        add <- stats::runif(1, 0.1, 0.5)
        R[i] <- R[i] + add
        check <- (L[i] <= Z[i] & Z[i] < R[i])
      }
    }
  }
  dados <- data.frame(Z, L, R, delta, cov_theta[,2:dim(cov_theta)[2]], cov_beta)
  return(dados)
}
