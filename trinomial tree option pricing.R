trinomTree = function(N, r, sigma, maturity, s0, k, isCall, isAmerican){
  delta.t = maturity / N    
  nu = r - (sigma^2)/2
  delta.x = sigma * sqrt(3*delta.t)
  discount.factor = exp(-r*delta.t)
  
  up.factor = exp(delta.x)
  down.factor = exp(-delta.x)
  
  probUp = 0.5*((sigma^2*delta.t + nu^2 * delta.t^2)/delta.x^2 + (nu*delta.t)/delta.x)
  probMid = 1 - ((sigma^2*delta.t + nu^2 * delta.t^2)/delta.x^2)
  probDown = 0.5 * ((sigma^2 * delta.t + nu^2 * delta.t^2)/delta.x^2 - (nu*delta.t)/delta.x)
  
  stock.tree = matrix(0, nrow = 2*N+1, ncol = N+1)
  
  #create stock price tree
  for(i in 1:(2*N+1)){
    for(j in 1:(N+1)){
      stock.tree[i,j] = s0 * down.factor^max((i-j),0) * up.factor^max((j-i),0)
    }
  }
  
  for(j in 1:(N+1)){
    if(j*2 < (2*N+1)){
      stock.tree[(j*2):(2*N+1),j] = NA
    }
  }
  
  
  #evaluate option tree at maturity
  option.tree = matrix(0, nrow = 2*N+1, ncol = N+1)
  if(isCall == TRUE){#call
    option.tree[, N+1] = stock.tree[, N+1] - k
    option.tree[which(option.tree[,N+1] < 0), N+1] = 0
  }else{#put
    option.tree[, N+1] = k - stock.tree[, N+1]  
    option.tree[which(option.tree[,N+1] < 0), N+1] = 0
  }
  
  for(j in N:1){
    i = 1
    while(!is.na(stock.tree[i,j])){
      if(isAmerican == FALSE){#european
        option.tree[i,j] = discount.factor * (probDown * option.tree[i+2,j+1] + 
                                                probMid * option.tree[i+1,j+1] +
                                                probUp * option.tree[i,j+1])
      }else{#american
        if(isCall == TRUE){#american call
          option.tree[i,j] = max(discount.factor * (probDown * option.tree[i+2,j+1] + 
                                                      probMid * option.tree[i+1,j+1] +
                                                      probUp * option.tree[i,j+1]), (stock.tree[i,j]-k))
        }else{#american put
          option.tree[i,j] = max(discount.factor * (probDown * option.tree[i+2,j+1] + 
                                                      probMid * option.tree[i+1,j+1] +
                                                      probUp * option.tree[i,j+1]), (k-stock.tree[i,j]))
        }
      }
      i = i + 1
    }
  }
  
  return(option.tree[1,1])
}

idx = 1
for(i in 1:6){
  s0 = 1003
  tau = timeToMat[i]
  trinomPriceEuro = vector(mode = "numeric", length = nrow(intersectData[[i]]))
  trinomPriceUs = vector(mode = "numeric", length = nrow(intersectData[[i]]))
  
  if((i %% 2) == 1){
    isCall = TRUE  
  }else{
    isCall = FALSE  
  }
  for(j in 1:(nrow(intersectData[[i]]))){
    trinomPriceEuro[j] = trinomTree(200, r, intersectData[[i]][j,"impliedVol"], tau, s0, intersectData[[i]][j, "strike"], isCall, FALSE)
    trinomPriceUs[j] = trinomTree(200, r, intersectData[[i]][j,"impliedVol"], tau, s0, intersectData[[i]][j, "strike"], isCall, TRUE)
  } 
  
  AMZNoptions[[idx]] = cbind(AMZNoptions[[idx]], trinomPriceEuro, trinomPriceUs)
  
  idx = idx + 1
}