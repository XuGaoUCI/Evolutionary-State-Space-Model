#import packages
install.packages("fields")
install.packages("Rsolnp")
require("fields")
require("Rsolnp")


#utility functions
###Estimation of the mixing matrix###
Shuffle_main <- function(Signals, Initializing){
  ###This function is to estimate the mixing matrix using resampling
  ###Input: "Signals" The observed LFPs signals
  ###       "n_channel" the number of nodes, 
  ###       "n_indpt" the number of latent sources
  ###       "ini_P" the initial covariance matrix for KF
  ###       "max_boots" the number of resampler 
  ###       "window" the size of the moving windwo for resampling,
  ###       "trials" the number of trials
  ###       "ini_state" the initial states of latent sources 
  ###       "max_iter" the number of loops for estimating mixing matrix
  ###Output: "mixing"  The estimate of mixing matrix
  ######################################################################
  ###Loglikelihood function for optimization###
  LogLike <- function(par){
    ###Kalman Filter step
    KF_procedure <- function (Phi1, Phi2, M, tau, sigma, Y, ini_state, ini_P){
      KF <- function(X, P, Phi1, Phi2, M, tau, sigma, Y){
        temp_Phi = rbind(cbind(Phi1, Phi2), cbind(diag(ncol(Phi1)), matrix(0, ncol(Phi2), ncol(Phi2))))
        temp_X = temp_Phi %*% X
        temp_P = temp_Phi %*% P %*% t(temp_Phi) + sigma^2*diag(c(rep(1,ncol(Phi1)), rep(0, ncol(Phi2))))
        con_M = cbind(M, zeros(nrow(M), ncol(M)))
        K_t = temp_P %*% t(con_M) %*% solve(con_M%*%temp_P%*%t(con_M) + tau^2*diag(rep(1,length(Y))))
        update_X = temp_X + K_t %*% (Y - con_M%*%temp_X)
        update_P = (diag(1, nrow(temp_X)) - K_t%*%con_M)%*%temp_P
        cov_matrix = con_M %*% temp_P %*% t(con_M) + tau^2*diag(1, nrow(con_M) , nrow(con_M))
        resid = Y - con_M %*% temp_X
        return(list("state" = update_X, "cov" = update_P, "cov_obs" = cov_matrix, "resid" = resid, 
                    "logdeter" = log(det(cov_matrix)), "qua" = t(resid) %*% solve(cov_matrix) %*% resid))
      }
      n <- ncol(Y)
      covv <- array(0, c(n, 2*nrow(Phi1), 2*nrow(Phi1)))
      State <- matrix(0, 2*nrow(Phi1), n)
      State[,1] <- ini_state
      covv[1,,] <- ini_P
      P <- ini_P
      kf.0 <- KF(ini_state, P, Phi1, Phi2, M, tau, sigma, Y[,1])
      Linn <- .5*(kf.0$logdeter+kf.0$qua)
      for (t in 2:n){
        kf.m <- KF(ini_state, P, Phi1, Phi2, M, tau, sigma, Y[,t])
        ini_state <- kf.m$state
        P <- kf.m$cov
        covv[t,,] <- P
        State[,t] <- ini_state
        Linn = Linn + .5*(kf.m$logdeter + kf.m$qua)
      }
      return(list("Covariance" = covv, "States" = State, "loglike" = Linn))
    }
    auto_par <- c(tran_AR2(1+par[1]/1000, 2/1000), tran_AR2(1+par[2]/1000, 8/1000), tran_AR2(1+par[3]/1000, 15/1000))
    Phi1 <- diag(auto_par[seq(1, length(auto_par), 2)]);
    Phi2 <- diag(auto_par[seq(1, length(auto_par), 2)+1]);
    tau <- par[4];
    sigma <- par[5];
    return(KF_procedure(Phi1, Phi2, M, tau, sigma, Y, ini_state, ini_P)$loglike)
  }
  ####frequency band###
  tran_AR2 <- function(mod, phase){
    phi1 <- (1/mod)*cos(phase)*2
    phi2 <- -1/(mod^2)
    return(c(phi1, phi2))
  }
  ###Estimation of mixing matrix###
  Mixing <- function(Y, State){
    temp <- sapply(1:nrow(Y),function(x){return(solve(State %*% t(State)) %*% State %*% Y[x,])})
    return(t(temp))
  }
  ###Kalman Filter Procedure####
  KF_procedure <- function (Phi1, Phi2, M, tau, sigma, Y, ini_state, ini_P){
    KF <- function(X, P, Phi1, Phi2, M, tau, sigma, Y){
      temp_Phi = rbind(cbind(Phi1, Phi2), cbind(diag(ncol(Phi1)), matrix(0, ncol(Phi2), ncol(Phi2))))
      temp_X = temp_Phi %*% X
      temp_P = temp_Phi %*% P %*% t(temp_Phi) + sigma^2*diag(c(rep(1,ncol(Phi1)), rep(0, ncol(Phi2))))
      con_M = cbind(M, zeros(nrow(M), ncol(M)))
      K_t = temp_P %*% t(con_M) %*% solve(con_M%*%temp_P%*%t(con_M) + tau^2*diag(rep(1,length(Y))))
      update_X = temp_X + K_t %*% (Y - con_M%*%temp_X)
      update_P = (diag(1, nrow(temp_X)) - K_t%*%con_M)%*%temp_P
      cov_matrix = con_M %*% temp_P %*% t(con_M) + tau^2*diag(1, nrow(con_M) , nrow(con_M))
      resid = Y - con_M %*% temp_X
      return(list("state" = update_X, "cov" = update_P, "cov_obs" = cov_matrix, "resid" = resid, 
                  "logdeter" = log(det(cov_matrix)), "qua" = t(resid) %*% solve(cov_matrix) %*% resid))
    }
    n <- ncol(Y)
    covv <- array(0, c(n, 2*nrow(Phi1), 2*nrow(Phi1)))
    State <- matrix(0, 2*nrow(Phi1), n)
    State[,1] <- ini_state
    covv[1,,] <- ini_P
    P <- ini_P
    Linn <- 0
    for (t in 2:n){
      kf.m <- KF(ini_state, P, Phi1, Phi2, M, tau, sigma, Y[,t])
      ini_state <- kf.m$state
      P <- kf.m$cov
      covv[t,,] <- P
      State[,t] <- ini_state
      Linn = Linn + .5*(kf.m$logdeter + kf.m$qua)
    }
    return(list("Covariance" = covv, "States" = State, "loglike" = Linn))
  }
  ###optimize steps###########
  max_iter <- Initializing$max_iter
  n_channel <- Initializing$n_channel
  n_indpt <- Initializing$n_indpt
  ini_P <- Initializing$ini_P
  max_boots <- Initializing$max_boots
  trials <- Initializing$trials
  window <- Initializing$window
  result_M <- array(0, c(n_channel, n_indpt, window, max_boots))
  result_S <- array(0, c(n_indpt, n_time, window, max_boots))
  for(rep in 1:max_boots){
    star <- sample(1:(trials-window+1),1)
    index_shu <- star:(star+window-1)
    if(max_boots == 1){
      realization_shuffle <-Signals
    } else {
      realization_shuffle <-Signals[,,index_shu]}
    ini_shuffle <- Initializing$ini_state[,,index_shu]
    err <-c()
    ini_M = matrix(0, n_channel, n_indpt)
    M = ini_M
    for (tri in 1:length(index_shu)) {
      if (length(index_shu) == 1){
        Y <- realization_shuffle
        ini_state <- c(ini_shuffle[,2],ini_shuffle[,1])
      } else {
        Y <- realization_shuffle[,,tri]
        ini_state <- c(ini_shuffle[,2,tri],ini_shuffle[,1,tri])
      }
      LB=c(1,1,1,0.1,0.1)
      UB=c(2,2,2,0.5,0.5)
      old_para = (LB+UB)*.5
      for(iter in 1:max_iter){
        temp = constrOptim((LB+UB)*.5, LogLike,NULL,ui = diag(c(rep(-1,3),1,1)), ci = c(rep(-2,3),0,0), metho="Nelder-Mead")
        auto_par <- c(tran_AR2(1+temp$par[1]/1000, (2/1000)*2*pi), tran_AR2(1+temp$par[2]/1000, (8/1000)*2*pi), tran_AR2(1+temp$par[3]/1000, (15/1000)*2*pi))
        Phi1 <- diag(auto_par[seq(1, length(auto_par), 2)]);
        Phi2 <- diag(auto_par[seq(1, length(auto_par), 2)+1]);
        tau <- abs(temp$par[4]);
        sigma <- abs(temp$par[5]);
        kf.mod <- KF_procedure(Phi1, Phi2, M, (tau), (sigma), Y, ini_state, ini_P)
        States <- kf.mod$States[1:3,]
        scale_States <- States/apply(States,1,sd)  ###rescale to unit variance identifiability problem
        new_M <- abs(Mixing(Y, scale_States))
        err[iter] <- sum((new_M - M)^2)
        M <- new_M
        if (sum((old_para-temp$par)^2)<1e-3) break
        old_para <- temp$par
        print(iter)
      }
      result_M[,,tri,rep] = new_M
      result_S[,,tri,rep] = scale_States
    }
    
  }
  return(list("matrix" = result_M))
}
###Reconstruct the latent signals###
Backevolve <- function(est_M, Signals, Initializing){
  ###This function is to estimate the latent sources
  ###Input: "est_M" the estimation of mixing matrix
  ###       "Signals" the observed LFPs signals
  ###       "n_channel" the number of nodes, 
  ###       "n_indpt" the number of latent sources
  ###       "ini_P" the initial covariance matrix for KF
  ###       "max_boots" the number of resampler 
  ###       "window" the size of the moving windwo for resampling,
  ###       "trials" the number of trials
  ###       "ini_state" the initial states of latent sources 
  ###       "max_iter" the number of loops for estimating mixing matrix
  ###Output: "latent"  the estimate of latent sources
  ###        "para"   the estimate of paramters
  ######################################################################
  ###Loglikelihood for optimization
  LogLike <- function(par){
    KF_procedure <- function (Phi1, Phi2, M, tau, sigma, Y, ini_state, ini_P){
      KF <- function(X, P, Phi1, Phi2, M, tau, sigma, Y){
        temp_Phi = rbind(cbind(Phi1, Phi2), cbind(diag(ncol(Phi1)), matrix(0, ncol(Phi2), ncol(Phi2))))
        temp_X = temp_Phi %*% X
        temp_P = temp_Phi %*% P %*% t(temp_Phi) + sigma^2*diag(c(rep(1,ncol(Phi1)), rep(0, ncol(Phi2))))
        con_M = cbind(M, zeros(nrow(M), ncol(M)))
        K_t = temp_P %*% t(con_M) %*% solve(con_M%*%temp_P%*%t(con_M) + tau^2*diag(rep(1,length(Y))))
        update_X = temp_X + K_t %*% (Y - con_M%*%temp_X)
        update_P = (diag(1, nrow(temp_X)) - K_t%*%con_M)%*%temp_P
        cov_matrix = con_M %*% temp_P %*% t(con_M) + tau^2*diag(1, nrow(con_M) , nrow(con_M))
        resid = Y - con_M %*% temp_X
        return(list("state" = update_X, "cov" = update_P, "cov_obs" = cov_matrix, "resid" = resid, 
                    "logdeter" = log(det(cov_matrix)), "qua" = t(resid) %*% solve(cov_matrix) %*% resid))
      }
      n <- ncol(Y)
      covv <- array(0, c(n, 2*nrow(Phi1), 2*nrow(Phi1)))
      State <- matrix(0, 2*nrow(Phi1), n)
      State[,1] <- ini_state
      covv[1,,] <- ini_P
      P <- ini_P
      Linn <- 0
      for (t in 2:n){
        kf.m <- KF(ini_state, P, Phi1, Phi2, M, tau, sigma, Y[,t])
        ini_state <- kf.m$state
        P <- kf.m$cov
        covv[t,,] <- P
        State[,t] <- ini_state
        Linn = Linn + .5*(kf.m$logdeter + kf.m$qua)
      }
      return(list("Covariance" = covv, "States" = State, "loglike" = Linn))
    }
    auto_par <- c(tran_AR2(1+par[1]/1000, 2/1000), tran_AR2(1+par[2]/1000, 8/1000), tran_AR2(1+par[3]/1000, 15/1000))
    Phi1 <- diag(auto_par[seq(1, length(auto_par), 2)]);
    Phi2 <- diag(auto_par[seq(1, length(auto_par), 2)+1]);
    tau <- par[4];
    sigma <- par[5];
    return(KF_procedure(Phi1, Phi2, M, tau, sigma, Y, ini_state, ini_P)$loglike)
  }
  ####frequency band#####
  tran_AR2 <- function(mod, phase){
    phi1 <- (1/mod)*cos(phase)*2
    phi2 <- -1/(mod^2)
    return(c(phi1, phi2))
  }
  ###function that with maxing matrix unknown###
  Mixing <- function(Y, State){
    temp <- sapply(1:nrow(Y),function(x){return(solve(State %*% t(State)) %*% State %*% Y[x,])})
    return(t(temp))
  }
  ##Kalman filter procedure
  KF_procedure <- function (Phi1, Phi2, M, tau, sigma, Y, ini_state, ini_P){
    KF <- function(X, P, Phi1, Phi2, M, tau, sigma, Y){
      temp_Phi = rbind(cbind(Phi1, Phi2), cbind(diag(ncol(Phi1)), matrix(0, ncol(Phi2), ncol(Phi2))))
      temp_X = temp_Phi %*% X
      temp_P = temp_Phi %*% P %*% t(temp_Phi) + sigma^2*diag(c(rep(1,ncol(Phi1)), rep(0, ncol(Phi2))))
      con_M = cbind(M, zeros(nrow(M), ncol(M)))
      K_t = temp_P %*% t(con_M) %*% solve(con_M%*%temp_P%*%t(con_M) + tau^2*diag(rep(1,length(Y))))
      update_X = temp_X + K_t %*% (Y - con_M%*%temp_X)
      update_P = (diag(1, nrow(temp_X)) - K_t%*%con_M)%*%temp_P
      cov_matrix = con_M %*% temp_P %*% t(con_M) + tau^2*diag(1, nrow(con_M) , nrow(con_M))
      resid = Y - con_M %*% temp_X
      return(list("state" = update_X, "cov" = update_P, "cov_obs" = cov_matrix, "resid" = resid, 
                  "logdeter" = log(det(cov_matrix)), "qua" = t(resid) %*% solve(cov_matrix) %*% resid))
    }
    n <- ncol(Y)
    covv <- array(0, c(n, 2*nrow(Phi1), 2*nrow(Phi1)))
    State <- matrix(0, 2*nrow(Phi1), n)
    State[,1] <- ini_state
    covv[1,,] <- ini_P
    P <- ini_P
    Linn <- 0
    for (t in 2:n){
      kf.m <- KF(ini_state, P, Phi1, Phi2, M, tau, sigma, Y[,t])
      ini_state <- kf.m$state
      P <- kf.m$cov
      covv[t,,] <- P
      State[,t] <- ini_state
      Linn = Linn + .5*(kf.m$logdeter + kf.m$qua)
    }
    return(list("Covariance" = covv, "States" = State, "loglike" = Linn))
  }
  
  ###estimation of latent sources###
  
  trials <- Initializing$trials
  ensem_para <- matrix(0, trials, 5)
  ensem_State <- array(0, c(3,1000, trials))
  M<- est_M
  temp <- list()
  for (tri in 1:trials){
    ini_state <- c(Initializing$ini_state[,2,tri],Initializing$ini_state[,1,tri])
    if(trials == 1){
      Y <- Signals
    }else{
      Y <- Signals[,,tri]}
    LB=c(1,1,1,0.1,0.1)
    UB=c(2,2,2,0.5,0.5)
    temp = constrOptim((LB+UB)*.5, LogLike,NULL,ui = diag(c(rep(-1,3),1,1)), ci = c(rep(-2,3),0,0))
    auto_par <- c(tran_AR2(1+temp$par[1]/1000, (2/1000)*2*pi), tran_AR2(1+temp$par[2]/1000, (8/1000)*2*pi), tran_AR2(1+temp$par[3]/1000, (15/1000)*2*pi))
    Phi1 <- diag(auto_par[seq(1, length(auto_par), 2)]);
    Phi2 <- diag(auto_par[seq(1, length(auto_par), 2)+1]);
    tau <- abs(temp$par[4]);
    sigma <- abs(temp$par[5]);
    ensem_para[tri,] <- temp$par
    kf.mod <- KF_procedure(Phi1, Phi2, M, (tau), (sigma), Y, ini_state, ini_P)
    States <- kf.mod$States[1:3,]
    scale_States <- States/apply(States,1,sd)  ###rescale to unit variance identifiability problem
    ensem_State[,,tri] <- scale_States
    print(tri)
  }
  return(list("latent" = ensem_State, "para" = ensem_para))
}
###Function transferring module to AR2 parameters
tran_AR2 <- function(mod, phase){
  phi1 <- (1/mod)*cos(phase)*2
  phi2 <- -1/(mod^2)
  return(c(phi1, phi2))
}
###Estimation of mixing matrix###
Mixing <- function(Y, State){
  temp <- sapply(1:nrow(Y),function(x){return(solve(State %*% t(State)) %*% State %*% Y[x,])})
  return(t(temp))
}
###Theoretical spectrum####
spec <- function(phi, sigma){
  temp = 2*phi[1]*(phi[2] - 1)*cos(2*pi*seq(-500:500)/1000) - 2*phi[2]*cos(4*pi*seq(-500:500)/1000)+1+sum(phi^2)
  return(sigma^2/temp)
}
zeros <- function(a,b){
  return(matrix(0, a, b))
}
###Time series generation
generate_time <- function(module, sigma, tau, n_channel, n_indpt, trials, n_time, M){
  con_Y <- array(0, dim = c(n_channel, n_time, trials))
  con_State <- array(0,  dim = c(n_indpt, n_time, trials))
  for (tri in 1:trials){
    phi1 = tran_AR2(module[1, tri], (2/1000)*2*pi)
    phi2 = tran_AR2(module[2, tri], (8/1000)*2*pi)
    phi3 = tran_AR2(module[3, tri], (15/1000)*2*pi)
    s1 <- arima.sim(model = list(ar = phi1), n = n_time,
                    innov = rnorm(n_time) * sigma[tri])
    s2 <- arima.sim(model = list(ar = phi2), n = n_time,
                    innov = rnorm(n_time) * sigma[tri])
    s3 <- arima.sim(model = list(ar = phi3), n = n_time,
                    innov = rnorm(n_time) * sigma[tri])
    
    S <- rbind(s1/sd(s1), s2/sd(s2), s3/sd(s3))
    con_State[,,tri] <- S
    con_Y[,,tri] <- M %*% S + tau[tri]*matrix(rnorm(n_time*nrow(M)), nrow(M), n_time) 
  }
  return(list("realization" = con_Y, "latent" = con_State))
}
###Modulus generation function
evol <- function(ini_mod, incre, trials, n_indpt){
  module <- matrix(0, n_indpt, trials)
  module[,1] <- ini_mod
  for (tri in 2:trials){
    module[,tri] <- module[,tri-1] + incre
  }
  return(module)
}

#generate data and initializing parameters
hp_para = read.csv("~/.../para_sim-3.csv")
true_M = read.csv("~/../mixingmatrix_sim-3.csv")

set.seed(2016)
phi1 = tran_AR2(1.001, (2/1000)*2*pi)
phi2 = tran_AR2(1.001, (15/1000)*2*pi)
phi3 = tran_AR2(1.001, (25/1000)*2*pi)
auto_par = c(phi1, phi2, phi3)
Phi1 = diag(auto_par[seq(1, length(auto_par), 2)]);
Phi2 = diag(auto_par[seq(1, length(auto_par), 2)+1]);
ini_P = rbind(cbind(Phi1, Phi2), cbind(diag(ncol(Phi1)), matrix(0, ncol(Phi2), ncol(Phi2)))) %*% rbind(cbind(Phi1, Phi2), cbind(diag(ncol(Phi1)), matrix(0, ncol(Phi2), ncol(Phi2))))
trials = 247
max_iter = 4
window = 2
max_boots = 4
n_time = 1000
n_channel = 12
n_indpt = 3
tau = hp_para[,4]
sigma = hp_para[,5]
module = t(hp_para[,1:3])
sim = generate_time(module, sigma, tau, n_channel, n_indpt, trials, n_time, as.matrix(true_M))
Signals = sim$realization
ini_state = sim$latent[,c(1,2),]
Initializing = list("n_channel" = n_channel, "n_indpt" = n_indpt, "ini_P" = ini_P, 
                      "max_boots" = max_boots, "window" = window,
                     "trials" = trials, "ini_state" = ini_state,"max_iter" = max_iter)

#Estimation steps
temp_real = Shuffle_main(Signals, Initializing)
avg_chan_est = apply(temp_real$matrix,c(1,2),mean)
sum((avg_chan_est - true_M)^2) #0.2237568
reconst = Backevolve(avg_chan_est,Signals, Initializing)

#temp = readRDS("test.rds")
