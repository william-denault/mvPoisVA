library(susiF.alpha)
library(sim1000G)
examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")

vcf = readVCF( vcf_file , maxNumberOfVariants = 200000 , min_maf = 0.02 , max_maf = NA )

# downloadGeneticMap( 4 )
readGeneticMap( chromosome = 4)

startSimulation( vcf )



'%!in%' <- function(x,y)!('%in%'(x,y))
id = c()
for(i in 1:100) id[i] = SIM$addUnrelatedIndividual()

# Show haplotype 1  of first 5 individuals
#print(SIM$gt1[1:5,1:6])

# Show haplotype 2
#print(SIM$gt1[1:5,1:6])



genotypes = SIM$gt1[1:100,] + SIM$gt2[1:100,]


caca <- genotypes
print(dim(genotypes))

str(genotypes)

if(file.exists("D:/Document/Serieux/Travail/Package/mvPoisVA/inst/check_pois_fsusie.RData")){
  load("D:/Document/Serieux/Travail/Package/mvPoisVA/inst/check_pois_fsusie.RData")

}else{
  res <-list()
}
pb_input <-list()

N <- 50
cal_purity <- function(l_cs,X){
  tt <- list()
  for (k in 1:length(l_cs)){
    if(length(unlist(l_cs[[k]]))==1 ){
      tt[[k]] <- 1
    }else{
      x <-abs( cor(X[,unlist(l_cs[[k]]   ) ]))


      tt[[k]] <-  min( x[col(x) != row(x)])
    }
  }
  return( mean( unlist(tt )))
}



input_pb <- list()
for( o in 1:2000){
  L <- sample(1:20, size=1)#actual number of effect
  lf <-  list()
  for(l in 1:L){
    lf[[l]] <- abs(simu_IBSS_per_level(lev_res=7)$sim_func) #functional effect for effect l
  }


  tt <- sample(0:4,1)
  G <- genotypes[sample( 1:nrow(genotype), size=N, replace=FALSE),]

  if( length(which(apply(G,2,var)==0))>0){
    G <- G[,-which(apply(G,2,var)==0)]
  }
  # G <- matrix( rnorm(100*300), nrow = 100)
  true_pos <- sample( 1:ncol(G), L)

  Y <- matrix(0 , ncol=  2^7 , nrow = N)
  for ( i in 1:N){
    for ( l in 1:L){
      predictor <- Y[i,]+ lf[[l]]*G[i,true_pos[[l]]]


    }
    Y[i,] <-   rpois(n= ncol(Y) ,
                     lambda =predictor  )
  }


  input <- list(Y=Y, X=G, true_pos=true_pos)
  m2 <- HF_susiF(Y=Y, X=G,L=20,L_start=11 ,nullweight=10  )
  m2$cs
  m1 <- susiF   (Y=Y, X=G,L=20,L_start=11 ,nullweight=10 )
  m1$cs
  m3 <- mv_Poisproc_reg(Y=Y, X=G,L=20,L_start=11 ,nullweight=10 , maxit=3 ,  verbose=FALSE)
  m3$cs


  out <- c( length(m1$cs), #number of CS
            length(which(true_pos%in% do.call(c, m1$cs))), #number of effect found
            Reduce("+",sapply(1:length(m1$cs), function(k)
              ifelse( length(which(true_pos%in%m1$cs[[k]] ))==0, 1,0)
            )
            ),#number of CS without any effect
            cal_purity(m1$cs, X=as.matrix(G)),#mean purity
            mean(sapply( m1$cs, length)), #CS size


            length(m2$cs),
            length(which(true_pos%in% do.call(c, m2$cs))),
            Reduce("+",sapply(1:length(m2$cs), function(k)
              ifelse( length(which(true_pos%in%m2$cs[[k]] ))==0, 1,0)
            )
            ),#number of CS without any effect
            cal_purity(m2$cs, X=as.matrix(G)),#mean purity
            mean(sapply( m2$cs, length)), #CS size



            length(m3$cs),
            length(which(true_pos%in% do.call(c, m3$cs))),
            Reduce("+",sapply(1:length(m3$cs), function(k)
              ifelse( length(which(true_pos%in%m3$cs[[k]] ))==0, 1,0)
            )
            ),#number of CS without any effect
            cal_purity(m3$cs, X=as.matrix(G)),#mean purity
            mean(sapply( m3$cs, length)), #CS size

            L,tt)
  res[[o]] <-unlist(out)
  if( out[3]>=1 ||out[8]>=1  || out[13]>=1   ){

    input <- list(Y=Y,
                  X=G,
                  true_pos=true_pos,
                  out=out)

    pb_input [[length(pb_input)+1]] <- input
    save(pb_input, file="D:/Document/Serieux/Travail/Package/mvPoisVA/inst/check_pb_pois_fsusie.R.RData")
  }

  save(res, file="D:/Document/Serieux/Travail/Package/mvPoisVA/inst/check_pois_fsusie.RData")


}

