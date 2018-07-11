crankNicholson = function(N, r, sigma, maturity, s0, k, isCall, divYield){
  startTime = Sys.time()
  delta.t = maturity / N    
  nu = r - (sigma^2)/2 - divYield
  delta.x = sigma * sqrt(3*delta.t)
  discount.factor = exp(-r*delta.t)
  
  up.factor = exp(delta.x)
  down.factor = exp(-delta.x)
  
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
  
  alpha = -0.25 * delta.t * (sigma^2/delta.x^2 + nu/delta.x)
  beta = 1 + delta.t * sigma^2 / (2*delta.x^2) + r*delta.t/2
  gamma = -0.25 * delta.t * (sigma^2/delta.x^2 - nu/delta.x)
  
  
  A = matrix(data = 0, nrow = 2*N+1, ncol = 2*N+1)
  A[1,1] = 1
  A[1,2] = -1
  A[2*N+1, 2*N] = 1
  A[2*N+1, 2*N+1] = -1
  
  for(i in 2:(2*N)){
    A[i,i-1] = alpha
    A[i,i] = beta
    A[i,i+1] = gamma
  }
  
  # optionColAtMaturity = as.numeric(option.tree[,N+1])
  # for(z in N:1){
  #   b = option.tree[,z+1]
  # 
  #   #initial guess is a vector of ones
  #   solution = rep.int(1, times = 2*N+1)
  #   prevSolution = rep.int(20, times = 2*N+1)
  #   #tolerance is used for solving linear system of equations
  #   #it is the stopping condition for jacobi iterations
  #   tolerance = 0.5
  #   #w is the relaxation parameter. typically it is between 1 & 2 for over relaxation
  #   w = 1.5
  # 
  #   while(normOfVector(prevSolution - solution) > tolerance){
  #     prevSolution = solution
  #     for(i in 1:(2*N+1)){
  #       summation1 = 0
  #       summation2 = 0
  #       for(j in 1:(2*N+1)){
  #         if((i != j) && (A[i,j] != 0) && (j<=i-1)){
  #           summation1 = summation1 + A[i,j] * solution[j]
  #         }
  #         if((i != j) && (A[i,j] != 0) && (j>=i+1)){
  #           summation2 = summation2 + A[i,j] * solution[j]
  #         }
  #       }
  #       if((i != 1) && (i != 2*N+1)){
  #         solution[i] = (1-w)*solution[i] + 
  #         (w/A[i,i]) * ((-alpha*b[i-1] + (-beta+2)*b[i] + -gamma*b[i+1]) - summation1 - summation2)
  #       }else{
  #         solution[i] = (1-w)*solution[i] + (w/A[i,i]) * (b[i] - summation1 - summation2)
  #       }
  #     }
  #   }
  #   option.tree[,z] = solution
  # }
  
  b = as.numeric(option.tree[,N+1])
  for(z in N:1){
    b = option.tree[,z+1]
    
    #initial guess is a vector of ones
    xCurr = rep.int(1, times = 2*N+1)
    xNext = vector(mode = "numeric", length = (2*N+1))
    #tolerance is used for solving linear system of equations
    #it is the stopping condition for jacobi iterations
    tolerance = 0.5
    
    while(normOfVector(xNext - xCurr) > tolerance){
      xCurr = xNext
      for(i in 1:(2*N+1)){
        summation = (A[i,] %*% xCurr) - A[i,i]*xCurr[i]
        if((i != 1) && (i != 2*N+1)){
          xNext[i] = (1/A[i,i]) * ((-alpha*b[i-1] + (-beta+2)*b[i] + -gamma*b[i+1]) - summation)
        }else{
          xNext[i] = (1/A[i,i]) * (b[i] - summation)
        }
      }
    }
    #after the loop xNext or xCurr holds option prices for a given time
    option.tree[,z] = xNext
  }
  
  endTime = Sys.time()
  timeTaken = endTime - startTime 
  return(c(as.numeric(option.tree["0", 1]), timeTaken))
}