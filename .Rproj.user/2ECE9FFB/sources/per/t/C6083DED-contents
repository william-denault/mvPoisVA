
library(smashr)


N=100
Tp <- 128
x <-matrix( rpois(N*Tp, lambda = 2), ncol=Tp)




    # Change matrix x to vector.
    x = as.vector(x[1,-1])
x_orig <-x




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

lol<- ParentTItable(x)$parent
n = length(x)
J = log2(n)

# Create the parent TI table for each signal, and put into rows of
# matrix y.
ls = sum(x)
plot(lol[1,reflect.res$idx],x_orig)


y = as.vector(t(ParentTItable(x)$parent))

zdat = withCallingHandlers(do.call(glm.approx, c(list(x = y, g = NULL),
                                                 glm.approx.param)))

res = list()

# Loop through resolutions, smoothing each resolution separately
for (j in 1:(J - lev)) {
  res = getlist.res(res, j, n, zdat, log, TRUE, ashparam)
}


recons = recons.mv(ls, res, log, n, J)
if (reflect==TRUE | floor(log2(length(x)))!=ceiling(log2(length(x)))) {
  recons$est.mean = recons$est.mean[reflect.indices]
  recons$est.var = recons$est.var[reflect.indices]
}

