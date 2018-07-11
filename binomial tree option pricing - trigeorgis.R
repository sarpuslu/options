#trigeorgis parameters
binomTree = function(N, r, sigma, maturity, s0, k, isCall, isAmerican){
  delta.t = maturity / N    
  nu = r - (sigma^2)/2
  delta.x = sqrt(nu^2 * delta.t^2 + sigma^2 * delta.t)
  discount.factor = exp(-r*delta.t)
  
  up.factor = exp(delta.x)
  down.factor = exp(-delta.x)
  
  probUp = 0.5 + (nu*delta.t)/(2*delta.x)
  probDown = 1 - probUp
  
  stock.tree = matrix(0, nrow = N+1, ncol = N+1)
  
  #create stock price tree
  for(i in (N+1):1){
    for(j in 1:i){
      stock.tree[j, i] = s0 * up.factor^(i-j) * down.factor^(j-1)
    }
  }
  
  #evaluate option tree at maturity
  option.tree = matrix(0, nrow = N+1, ncol = N+1)
  if(isCall == TRUE){#call
    option.tree[, N+1] = stock.tree[, N+1] - k
    option.tree[which(option.tree[,N+1] < 0), N+1] = 0
  }else{#put
    option.tree[, N+1] = k - stock.tree[, N+1]  
    option.tree[which(option.tree[,N+1] < 0), N+1] = 0
  }
  
  for(i in N:1){
    for(j in 1:i){
      if(isAmerican == FALSE){#european options
        option.tree[j, i] = discount.factor * (probUp * option.tree[j,i+1] + probDown *                option.tree[j+1,i+1])
      }else{#american options
        if(isCall == TRUE){
          option.tree[j, i] = max(discount.factor * (probUp * option.tree[j,i+1] + probDown *                option.tree[j+1,i+1]), (stock.tree[j,i]-k))
        }else{
          option.tree[j, i] = max(discount.factor * (probUp * option.tree[j,i+1] + probDown *                option.tree[j+1,i+1]), (k-stock.tree[j,i]))
        }
      }
    }
  }
  
  return(option.tree[1,1])
}