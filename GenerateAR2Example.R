###Function transferring module to AR2 parameters
tran_AR2 <- function(mod, phase){
  phi1 <- (1/mod)*cos(phase)*2
  phi2 <- -1/(mod^2)
  return(c(phi1, phi2))
}


###############################################
###############################################

#simultion 1###################################
phi1 = tran_AR2(1.001, (2/1000)*2*pi)  #delta
phi2 = tran_AR2(1.001, (6/1000)*2*pi)  #theta
phi3 = tran_AR2(1.001, (8/1000)*2*pi) #alpha
phi4 = tran_AR2(1.001, (15/1000)*2*pi) #beta
phi5 = tran_AR2(1.001, (32/1000)*2*pi) #gamma 
auto_par = c(phi1, phi2, phi3, phi4, phi5)
