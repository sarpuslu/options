trapezoid = function(fun, a, b, n) {
  # numerical integral of fun from a to b
  # using the trapezoid rule with n subdivisions
  # assume a < b and n is a positive integer
  h = (b-a)/n
  x = seq(a, b, by=h)
  y = vector(mode = "numeric", length = length(x))
  for(i in 1:length(x)){
    y[i] = fun(x[i])
  }
  s = h * (y[1]/2 + sum(y[2:n]) + y[n+1]/2)
  return(s)
}