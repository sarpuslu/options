#implied volatility is the sigma st someFunc ~= 0
someFunc = function(sigma, k, bid, ask){
  return(BSM(tau, s0, r, sigma, k, isCall) - (bid + ask)/2)
}

#implementation of bisection method
bisection = function(lowerBound, upperBound, tolerance, k, bid, ask){
  count = 0 #counting the number of steps it takes to prevent infinite loops
  while(abs(upperBound - lowerBound) > tolerance && count < 10000){
    mid = (lowerBound + upperBound) / 2
    
    if((someFunc(lowerBound, k, bid, ask)*someFunc(mid, k, bid, ask)) < 0){
      upperBound = mid 
    }
    else if((someFunc(mid, k, bid, ask)*someFunc(upperBound, k, bid, ask)) < 0){
      lowerBound = mid  
    }
    
    count = count + 1
  }
  
  impliedVolAndCount = c(lowerBound, count)
  #if no solution is obtained by 10000 steps return error code -1
  if(count >= 10000){
    return(-1)  
  }else{
    return(impliedVolAndCount)
    
  }
}