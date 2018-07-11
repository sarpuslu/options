explicitGreeks = function(N, r, sigma, maturity, s0, k, isCall, divYield){
  delta.t = maturity / N    
  nu = r - (sigma^2)/2 - divYield
  delta.x = sigma * sqrt(3*delta.t)
  discount.factor = exp(-r*delta.t)
  
  up.factor = exp(delta.x)
  down.factor = exp(-delta.x)
  
  probUp = (sigma^2 * delta.t)/(2*delta.x^2) + (nu*delta.t)/(2*delta.x)
  probMid = 1 - (sigma^2*delta.t)/delta.x^2
  probDown = (sigma^2 * delta.t)/(2*delta.x^2) - (nu*delta.t)/(2*delta.x)
  
  stock.tree = matrix(0, nrow = 2*N+1, ncol = N+1)
  rownames(stock.tree) = seq(from = N, to = -N, by = -1)
  stock.tree["0",] = s0
  
  #populate upper triangle of stock tree
  for(j in 2:(N+1)){
    for(i in 1:(j-1)){
      row = as.character(i)
      stock.tree[row, j] = s0 * up.factor^i
    }  
  }
  #populate lower triange of stock tree
  for(j in 2:(N+1)){
    for(i in -1:-(j-1)){
      row = as.character(i)
      stock.tree[row, j] = s0 * down.factor^abs(i)
    }  
  }
  
  #evaluate option tree at maturity
  option.tree = matrix(0, nrow = 2*N+1, ncol = N+1)
  rownames(option.tree) = seq(from = N, to = -N, by = -1)
  if(isCall == TRUE){#call
    option.tree[, N+1] = stock.tree[, N+1] - k
    option.tree[which(option.tree[,N+1] < 0), N+1] = 0
  }else{#put
    option.tree[, N+1] = k - stock.tree[, N+1]  
    option.tree[which(option.tree[,N+1] < 0), N+1] = 0
  }
  
  for(j in N:1){
    for(i in seq(from = N-1, to = -(N-1), by = -1)){
      row = as.character(i)
      rowUp = as.character(i+1)
      rowDown = as.character(i-1)
      option.tree[row,j] = discount.factor * (probUp * as.numeric(option.tree[rowUp, j+1]) +
                                                probMid * as.numeric(option.tree[row, j+1]) +
                                                probDown * as.numeric(option.tree[rowDown, j+1]))
    }
    #upper boundary
    if(isCall == TRUE){
      option.tree[1, j] = option.tree[2, j] + (stock.tree[1,N+1] - stock.tree[2,N+1])  
    }else{
      option.tree[1, j] = option.tree[2, j]    
    }
    #lower boundary
    if(isCall == TRUE){
      option.tree[2*N+1, j] = option.tree[2*N, j]  
    }else{
      option.tree[2*N+1, j] = option.tree[2*N,j] - (stock.tree[2*N+1,N] - stock.tree[2*N,N])
    }
  }
  
  
  delta = (option.tree["1",1] - option.tree["-1",1])/(stock.tree["1",N] - stock.tree["-1", N])
  gamma = ((option.tree["1", 1] - option.tree["0", 1])/(stock.tree["1",N] - stock.tree["0", N])
           - (option.tree["0", 1] - option.tree["-1", 1])/(stock.tree["0",N] - stock.tree["-1", N]))/ (0.5 * (stock.tree["1",N] - stock.tree["-1", N]))
  theta = (option.tree["0", 2] - option.tree["0",1])/delta.t
  vega = (explicitFiniteDiff(N, r, (sigma + 0.001*sigma), maturity, s0, k, isCall, divYield)
          - explicitFiniteDiff(N, r, (sigma - 0.001*sigma), maturity, s0, k, isCall, divYield))/ (2*0.001*sigma)
  rho = (explicitFiniteDiff(N, r + 0.001*r, sigma, maturity, s0, k, isCall, divYield)
         - explicitFiniteDiff(N, r + 0.001*r, sigma, maturity, s0, k, isCall, divYield))/ (2*0.001*r)
  greeks = c(delta, gamma, theta, vega, rho)
  names(greeks) = c("delta", "gamma", "theta", "vega", "rho")
  return(greeks)
}