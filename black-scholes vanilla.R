BSM = function(tau, s0, r, sigma, k, isCall){
  d1 = (log(s0/k, base = exp(1)) + ((r +(sigma)^2)/2)*tau)/(sigma*sqrt(tau))
  d2 = d1 - (sigma*sqrt(tau))  
  
  if(isCall == TRUE){
    optionPrice = s0 * pnorm(d1) - k * exp(-r * tau) * pnorm(d2)
  }else{
    optionPrice = k * exp(-r * tau) * pnorm(-d2) - s0 * pnorm(-d1)  
  }
  
  return(optionPrice)
}