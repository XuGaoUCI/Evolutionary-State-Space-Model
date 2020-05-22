#import packages
install.packages("fields")
install.packages("Rsolnp")
require("fields")
require("Rsolnp")

#utility functions
###Estimation of the mixing matrix###
Shuffle_main <- function(Signals, Initializing){
  ###This function is to estimate the mixing matrix using resampling
  ###Input: "Signals" the observed LFPs signals
  ###       "n_channel" the number of nodes, 
  ###       "n_indpt" the number of latent sources
  ###       "ini_P" the initial covariance matrix for KF
  ###       "max_boots" the number of resampler 
  ###       "window" the size of the moving windwo for resampling,
  ###       "trials" the number of trials
  ###       "ini_state" the initial states of latent sources 
  ###       "max_iter" the number of loops for estimating mixing matrix
  ###Output: "mixing"  the estimate of mixing matrix
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
    auto_par <- c(tran_AR2(1+par[1]/1000, 2/1000), tran_AR2(1+par[2]/1000, 6/1000), tran_AR2(1+par[3]/1000, 8/1000),  tran_AR2(1+par[4]/1000, 15/1000),  tran_AR2(1+par[5]/1000, 32/1000))
    Phi1 <- diag(auto_par[seq(1, length(auto_par), 2)]);
    Phi2 <- diag(auto_par[seq(1, length(auto_par), 2)+1]);
    tau <- par[6];
    sigma <- par[7];
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
      LB=c(1,1,1,1,1,0.1,0.1)
      UB=c(2,2,2,2,2,0.5,0.5)
      old_para = (LB+UB)*.5
      for(iter in 1:max_iter){
        temp = constrOptim((LB+UB)*.5, LogLike,NULL,ui = diag(c(rep(-1,5),1,1)), ci = c(rep(-2,5),0,0), metho="Nelder-Mead")
        auto_par <- c(tran_AR2(1+temp$par[1]/1000, (2/1000)*2*pi), tran_AR2(1+temp$par[2]/1000, (6/1000)*2*pi), tran_AR2(1+temp$par[3]/1000, (8/1000)*2*pi), tran_AR2(1+temp$par[4]/1000, (15/1000)*2*pi), tran_AR2(1+temp$par[5]/1000, (32/1000)*2*pi))
        Phi1 <- diag(auto_par[seq(1, length(auto_par), 2)]);
        Phi2 <- diag(auto_par[seq(1, length(auto_par), 2)+1]);
        tau <- abs(temp$par[6]);
        sigma <- abs(temp$par[7]);
        kf.mod <- KF_procedure(Phi1, Phi2, M, (tau), (sigma), Y, ini_state, ini_P)
        States <- kf.mod$States[1:5,]
        scale_States <- States/apply(States,1,sd)  ###rescale to unit variance identifiability problem
        new_M <- abs(Mixing(Y, scale_States))
        err[iter] <- sum((new_M - M)^2)
        M <- new_M
        if (sum((old_para-temp$par)^2)<1e-3) break
        old_para <- temp$par
        print(iter)
        print(err)
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
    auto_par <- c(tran_AR2(1+par[1]/1000, 2/1000), tran_AR2(1+par[2]/1000, 6/1000), tran_AR2(1+par[3]/1000, 8/1000), tran_AR2(1+par[4]/1000, 15/1000), tran_AR2(1+par[5]/1000, 32/1000))
    Phi1 <- diag(auto_par[seq(1, length(auto_par), 2)]);
    Phi2 <- diag(auto_par[seq(1, length(auto_par), 2)+1]);
    tau <- par[6];
    sigma <- par[7];
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
  ensem_para <- matrix(0, trials, 7)
  ensem_State <- array(0, c(5,1000, trials))
  M<- est_M
  temp <- list()
  for (tri in 1:trials){
    ini_state <- c(Initializing$ini_state[,2,tri],Initializing$ini_state[,1,tri])
    if(trials == 1){
      Y <- Signals
    }else{
    Y <- Signals[,,tri]}
    LB=c(1,1,1,1,1,0.1,0.1)
    UB=c(2,2,2,2,2,0.5,0.5)
    temp = constrOptim((LB+UB)*.5, LogLike,NULL,ui = diag(c(rep(-1,5),1,1)), ci = c(rep(-2,5),0,0))
    auto_par <- c(tran_AR2(1+temp$par[1]/1000, (2/1000)*2*pi), tran_AR2(1+temp$par[2]/1000, (6/1000)*2*pi), tran_AR2(1+temp$par[3]/1000, (8/1000)*2*pi), tran_AR2(1+temp$par[4]/1000, (15/1000)*2*pi), tran_AR2(1+temp$par[5]/1000, (32/1000)*2*pi))
    Phi1 <- diag(auto_par[seq(1, length(auto_par), 2)]);
    Phi2 <- diag(auto_par[seq(1, length(auto_par), 2)+1]);
    tau <- abs(temp$par[6]);
    sigma <- abs(temp$par[7]);
    ensem_para[tri,] <- temp$par
    kf.mod <- KF_procedure(Phi1, Phi2, M, (tau), (sigma), Y, ini_state, ini_P)
    States <- kf.mod$States[1:5,]
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
    phi1 = tran_AR2(module[1, tri], (2/1000)*2*pi) #delta
    phi2 = tran_AR2(module[2, tri], (6/1000)*2*pi) #theta
    phi3 = tran_AR2(module[3, tri], (8/1000)*2*pi) #alpha
    phi4 = tran_AR2(module[4, tri], (15/1000)*2*pi) #beta
    phi5 = tran_AR2(module[5, tri], (32/1000)*2*pi) #gamma 
    
    s1 <- arima.sim(model = list(ar = phi1), n = n_time,
                    innov = rnorm(n_time) * sigma[tri])
    s2 <- arima.sim(model = list(ar = phi2), n = n_time,
                    innov = rnorm(n_time) * sigma[tri])
    s3 <- arima.sim(model = list(ar = phi3), n = n_time,
                    innov = rnorm(n_time) * sigma[tri])
    s4 <- arima.sim(model = list(ar = phi4), n = n_time,
                    innov = rnorm(n_time) * sigma[tri])
    s5 <- arima.sim(model = list(ar = phi5), n = n_time,
                    innov = rnorm(n_time) * sigma[tri])
    
    
    
    S <- rbind(s1/sd(s1), s2/sd(s2), s3/sd(s3), s4/sd(s4), s5/sd(s5))
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




###generate the time series
###initializing########
set.seed(2016)
n_time <- 1000
trials <- 20
n_channel <- 20
n_indpt <- 5
tau <- rep(1, trials)
sigma <- rep(.5, trials)
module <- evol(rep(1.001, 5), 0.00005, trials, n_indpt)
true_M <- matrix(runif(n_channel*n_indpt, min = 0, max = 1), n_channel, n_indpt)
true_M[,2] <- true_M[,4] <- 0
sim <- generate_time(module, sigma, tau, n_channel, n_indpt, trials, n_time, true_M)
###initializing parameters###

###############################################
###############################################

#simultion 1###################################
phi1 = tran_AR2(module[1,1], (2/1000)*2*pi)  #delta
phi2 = tran_AR2(module[2,1], (6/1000)*2*pi)  #theta
phi3 = tran_AR2(module[3,1], (8/1000)*2*pi) #alpha
phi4 = tran_AR2(module[4,1], (15/1000)*2*pi) #beta
phi5 = tran_AR2(module[5,1], (32/1000)*2*pi) #gamma 
auto_par = c(phi1, phi2, phi3, phi4, phi5)
Phi1 <- diag(auto_par[seq(1, length(auto_par), 2)]);
Phi2 <- diag(auto_par[seq(1, length(auto_par), 2)+1]);
ini_P <- rbind(cbind(Phi1, Phi2), cbind(diag(ncol(Phi1)), matrix(0, ncol(Phi2), ncol(Phi2)))) %*% rbind(cbind(Phi1, Phi2), cbind(diag(ncol(Phi1)), matrix(0, ncol(Phi2), ncol(Phi2))))

max_iter = 5
max_boots = 1
window = 1
Signals <- sim$realization[,,5]
n_indpt_updated <- 5
Initializing <- list("n_channel" = n_channel, "n_indpt" = n_indpt_updated, "ini_P" = ini_P, 
                     "max_boots" = max_boots, "window" = window,
                     "trials" = 1, "ini_state" = sim$latent[,1:2,], "max_iter" = max_iter)
#estimation of mixing matrix
temp_resu <- Shuffle_main(Signals, Initializing)

avg_chan = apply(temp_resu$matrix,c(1,2),mean)
sum((avg_chan - true_M)^2)
#estimation of latent sources
reconst <- Backevolve(avg_chan,Signals, Initializing)
#plots in the paper
peri_true = matrix(0, 5, 1000)
for(i in 1:5){
  peri_true[i,]= scale(abs(fft(sim$latent[i,,5])/sqrt(1000))^2)
}
peri_est = matrix(0, 5, 1000)
for(i in 1:5){
  peri_est[i,]= scale(abs(fft(reconst$latent[i,,1])/sqrt(1000))^2)

}
par(mfrow=c(3,1))
plot(peri_true[1,1:50], type="l", xlab = "Frequency in Hz", ylab="Power", cex.lab=1.5,cex.axis=1.5, lwd=2)
lines(peri_est[1,1:50],col=c(2), lwd=.5)
plot(peri_true[2,1:50], type="l", xlab = "Frequency in Hz", ylab="Power", cex.lab=1.5,cex.axis=1.5, lwd=2)
lines(peri_est[2,1:50],col=c(2), lwd=.5)
plot(peri_true[3,1:50], type="l", xlab = "Frequency in Hz", ylab="Power", cex.lab=1.5,cex.axis=1.5, lwd=2)
lines(peri_est[3,1:50],col=c(2), lwd=.5)
plot(peri_true[4,1:50], type="l", xlab = "Frequency in Hz", ylab="Power", cex.lab=1.5,cex.axis=1.5, lwd=2)
lines(peri_est[4,1:50],col=c(2), lwd=.5)
plot(peri_true[5,1:50], type="l", xlab = "Frequency in Hz", ylab="Power", cex.lab=1.5,cex.axis=1.5, lwd=2)
lines(peri_est[5,1:50],col=c(2), lwd=.5)



#simulation 2##################################
phi1 = tran_AR2(module[1,1], (2/1000)*2*pi)  #delta
phi2 = tran_AR2(module[2,1], (6/1000)*2*pi)  #theta
phi3 = tran_AR2(module[3,1], (8/1000)*2*pi) #alpha
phi4 = tran_AR2(module[4,1], (15/1000)*2*pi) #beta
phi5 = tran_AR2(module[5,1], (32/1000)*2*pi) #gamma 
auto_par = c(phi1, phi2, phi3, phi4, phi5)
Phi1 <- diag(auto_par[seq(1, length(auto_par), 2)]);
Phi2 <- diag(auto_par[seq(1, length(auto_par), 2)+1]);
ini_P <- rbind(cbind(Phi1, Phi2), cbind(diag(ncol(Phi1)), matrix(0, ncol(Phi2), ncol(Phi2)))) %*% rbind(cbind(Phi1, Phi2), cbind(diag(ncol(Phi1)), matrix(0, ncol(Phi2), ncol(Phi2))))

max_iter = 5
max_boots = 3
window = 5
Signals <- sim$realization
Initializing <- list("n_channel" = n_channel, "n_indpt" = n_indpt, "ini_P" = ini_P, 
                     "max_boots" = max_boots, "window" = window,
                     "trials" = trials, "ini_state" = sim$latent[,1:2,], "max_iter" = max_iter)

#estimation procedure
temp_resu <- Shuffle_main(Signals, Initializing)
avg_chan = apply(temp_resu$matrix,c(1,2),mean)
sum((avg_chan - true_M)^2) #0.02813096
reconst = Backevolve(avg_chan,Signals, Initializing)


#heatmap plot
#heatmap geneareted for channel 1
par(mfrow=c(1,1))
par(mar = c(5,4,4,11))  
hmcols = colorRampPalette(c("blue","red"))(256)
peri_multi = matrix(0, dim(sim$realization)[3], 1000)
for(tri in 1:dim(sim$realization)[3]){
  peri_multi[tri,]= abs(fft(sim$realization[1,,tri])/sqrt(1000))^2
}
x = 1:nrow(peri_multi[,1:40])
y = (1:ncol(peri_multi[,1:40]))
temp = (peri_multi[,1:40])
image(x, y, temp,xaxt = "n",xlab = "Epochs",ylab = "Frequency in Hertz", col = hmcols)
axis(1,c(1,20,40, 60, 80, 100),c(1,20,40, 60, 80, 100))
par(mfrow=c(1,1))
image.plot(((temp)), legend.only = T, col = hmcols)

#periodogram of true and est 3 bands
par(mar = c(5,4.5,4,2))  
#simulation heatmap true
par(mfrow=c(3,1))
for(i in c(1,3,5)){
  hmcols = colorRampPalette(c("blue","red"))(256)
  peri_multi = matrix(0, dim(sim$latent)[3], 1000)
  for(tri in 1:dim(sim$latent)[3]){
    peri_multi[tri,]= abs(fft(sim$latent[i,,tri])/sqrt(1000))^2
  }
  x = 1:nrow(peri_multi[,1:40])
  y = (1:ncol(peri_multi[,1:40]))
  temp = (peri_multi[,1:40])
  image(x, y, temp,  xaxt="n", xlab="Epochs", ylab = "Frequency in Hertz", col = hmcols, cex.lab=2,cex.axis=1.8)
  
}
par(mfrow=c(1,1))
#simulation heatmap est
par(mar = c(5,4.5,4,12))  
par(mfrow=c(3,1))
for(i in c(1,3,5)){
  hmcols = colorRampPalette(c("blue","red"))(256)
  peri_multi = matrix(0, dim(reconst$latent)[3], 1000)
  for(tri in 1:dim(reconst$latent)[3]){
    peri_multi[tri,]= abs(fft(reconst$latent[i,,tri])/sqrt(1000))^2
  }
  x = 1:nrow(peri_multi[,1:40])
  y = (1:ncol(peri_multi[,1:40]))
  temp = (peri_multi[,1:40])
  image(x, y, temp,  xaxt="n", xlab="Epochs", ylab = "Frequency in Hertz", col = hmcols,cex.lab=2,cex.axis=2.2)
}
par(mfrow=c(1,1))
image.plot(((temp)), legend.only = T, col = hmcols)



