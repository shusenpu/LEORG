############### MLE ############
### fix alpha
MLE_fun_param2 <- function(param){
  alpha = 1.5
  beta = param[1]
  gamma = param[2]
  k = param[3]
  LnL <- n*log(alpha) + n*log(beta) + n*log(k) + 
    sum(log(gamma*exp(-gamma*x))) + 
    (beta-1)*sum(log(1-exp(-gamma*x))) - (beta+1)*sum(-gamma*x) -
    (k+1)*sum(log(1+alpha*(exp(gamma*x)-1)^beta))
  LnL <- -LnL
  return(LnL)
}




LS_fun_param2 <- function(param) {
  alpha = 1.5
  beta = param[1]
  gamma = param[2]
  k = param[3]
  LS <- sum((1-(1+alpha*(1/(1-exp(-gamma*x))-1)**(-beta))**(-k)-ls_ins)**2)
  LS
}




CVM_fun_param2 <- function(param) { 
  alpha = 1.5
  beta = param[1]
  gamma = param[2]
  k = param[3]
  CVM <- 1/n*sum((1-(1+alpha*(1/(exp(gamma*x)-1))^beta)^(-k)-ins)**2)
  CVM
} 



ADE_fun_param2 <- function(param) {
  alpha = 1.5
  beta = param[1]
  gamma = param[2]
  k = param[3]
  ADE <- -sum(i21*(log(1-(1+alpha*(exp(gamma*x)-1)^(beta))^(-k))
                   +log((1+alpha*(exp(gamma*x_rev)-1)^(beta))^(-k))))
  ADE
}
