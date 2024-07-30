tt <-  Y[1,]
plot(tt)

emp_prop <- Y_min[1,]/Y_tot[1,]
emp_log_prop <- log(emp_prop[- length(emp_prop)])
plot(emp_prop)
emp_lambda_tot <- Y_tot[1,ncol(Y_tot)]
emp_log_lambda_tot <- log(emp_lambda_tot )

tt <- reverse_intensity_transform(vec_int = c(log( sigmoid(emp_prop[-length(emp_prop)])) ,
                                              emp_lambda_tot),
                                  indx_lst = indx_lst,
                                  is.logprob=TRUE,
                                  is.log_int =TRUE)





plot( tt)
lines(Y[1,])

reverse_intensity_transform2 <- function(vec_int,indx_lst ){


  vec_int = c(emp_prop[-length(emp_prop)] ,
              emp_lambda_tot)
  J=log2(length(vec_int))

  lp <- log(vec_int[- length(vec_int)])
  lq = log(1-pmin(exp(lp),1-1e-10))



  print( "assuming input is total intensity ")
  log_lambda_tot <- log(vec_int[  length(vec_int)])




  out <- rep( log_lambda_tot , 2^J)

  for(s in (J ):1){

    nD = 2^(J-s+1)
    nDo2 = nD/2
    tt <-1
    for(l in 0:(2^(s-1)-1)){
      ind = (l*nD+1):((l+1)*nD) # all "sub index for coef s,l (here s=D)

      # print(ind)
      ind_l <-  ind[1:nDo2] #all "sub index in the left for coef s,l (here s=D)
      ind_r <-  ind[(nDo2+1):nD] # all "sub index in the right for coef s,l (here s=D)
      out[ind_l] <- out[ind_l]+ lp[indx_lst[[(s )]][tt]]
      out[ind_r] <- out[ind_r]+ lq[indx_lst[[(s )]][tt]]
      tt <- tt+1
    }
  }
  return(out)
}


plot(  reverse_intensity_transform2(vec_int = c(emp_prop[-length(emp_prop)] ,
                                                   emp_lambda_tot),
                                       indx_lst = indx_lst))
