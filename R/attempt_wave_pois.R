
library(smashr)


N=100
Tp <- 128
x <-matrix( rpois(N*Tp, lambda = 2), ncol=Tp)




# Change matrix x to vector.
x = as.vector(x[1,-1])





J = log2(length(x))

# If ncol(x) is not a power of 2, reflect x.
if ((J %% 1) != 0)
{
  reflect = TRUE
}

# Reflect signal.
if (reflect ) {
  reflect.res     = reflect(x)
  reflect.indices = reflect.res$idx
  x               = reflect.res$x
}


reflect.indices

n = length(x)
J = log2(n)
