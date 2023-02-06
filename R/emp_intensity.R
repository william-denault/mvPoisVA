
n = length(sig)
J = log2(n)
dmat = matrix(NA, nrow = J + 1, ncol = n)
dmat[1, ] = sig



ind = (l * nD + 1):((l + 1) * nD)
ind2 = (l * twonD + 1):((l + 1) * twonD)
x = dmat[D + 1, ind]
lsumx = x[seq(from = 1, to = nD - 1, by = 2)] +
  x[seq(from = 2, to = nD, by = 2)]
rx = rshift(x)
rsumx = rx[seq(from = 1, to = nD - 1, by = 2)] +
  rx[seq(from = 2, to = nD, by = 2)]
dmat[D + 2, ind] = c(lsumx, rsumx)
dmat2[D + 1, ind2] = c(x, rx)



