simpson = function(fun, a, b, n) {
  # numerical integral using Simpson's rule
  # assume a < b and n is an even positive integer
  h = (b-a)/n
  x = seq(a, b, by=h)
  sum = 0
  for(i in 2:length(x)){
    sum = sum + fun(x[i-1]) + 4*fun((x[i] + x[i-1])/2) + fun(x[i])
  }
  return(h/6 * sum)
}

simpsonsError = function(a, N){
  return(pi - simpson(f, -a, a, N))
}