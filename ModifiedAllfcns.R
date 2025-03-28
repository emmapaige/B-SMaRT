# Install packages
library(coda)
library(parallel)
library(bmixture)
library(fields)
library(matrixStats) 
library(gtools)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(dplyr)
library(reshape2)

#############################
#Random parameters (J=5) #
#############################
# Generate newP from dirichlet(p+main)
NucleoProb = function( J, J.main = 4 ){ 
  #J: number of possible reads; J.main: number of possible invariant reads 
  p = runif( J )
  main = rep( 0, J )
  r = runif(1)
  if( r > .5 ){ main[ sample.int( J.main, size = 1 ) ] = 100  
  }else if( r > .25 ){ main[ sample.int( J.main, size = 2 ) ] = c(25, 75)  
  }else{ main[ sample.int( J.main, size = 2 ) ] = c(50, 50) } 
  newp = rdirichlet( 1, p + main ) 
  return( newp )  
}
#Let the second largest element of newp be 10, and normalize newp
NewNucleoProb = function( oldp, read.main = 4 ){
  oldp.main = oldp[ 1:read.main ]
  oldp.2max = which( oldp.main == sort( oldp.main, decreasing = T )[2] )
  newp = oldp; newp[ oldp.2max ] = 10; newp = newp / sum( newp )
  return( newp )
}


##################################
#Consolidate homogenous positions#
##################################
#g0[1:J] records the invariant read sites. g0[J+1] records the noninvariant read sites.
consoData<-function( Data )
{
  J = nrow( Data )
  N0 = ncol( Data ) 
  C0 = rep( 0, N0 )
  for ( j in 1:J ){
    for ( i in 1:N0 ){
      if ( sum( Data[ -j, i ] == 0 ) == ( J - 1 ) ){ C0[ i ] = j } #record invariant read position
    }
  }  
  g0 = list( NULL ) #Original indices for the consolidated sites
  na.list = list( NULL )
  newD = matrix( NA, nrow = J, ncol = J )
  for ( j in 1:J ){
    g0[[ j ]] = which( C0 == j )
    l = length( which( C0 == j ) )
    if ( l == 0 ){ na.list
    }else if ( l == 1 ){ newD[ , j ] = Data[ , g0[[ j ]] ]
    }else{ newD[ , j ] = rowSums( Data[ , g0[[ j ]] ] ) #consolidate invariant read sites
    }
  }
  g0[[ J+1 ]] = which( C0==0 )
  g0.length = lapply( g0, length )
  emptyg0 = which( g0.length == 0 )
  if ( length( emptyg0 ) > 0 ){ newD = newD[ , -emptyg0 ] }  
  newData = cbind( newD, Data[ , g0[[ J + 1 ]] ] )
  return( list( newData, g0 ) )
}

#########################################
#Create test data                       #          
#K=15; Substitutions=5; n=300; N=1200;  #
#########################################
#What are these numbers? What are the outputs?
#K0 is group of component; newK is number of mutation
#d: number of possible reads
#n: sample size
#tests: n*tests=total sample size
#sitect: multinomial(sitect,p)
#postD: the test with mutations
#test1, 2, 3 have no mutation; test 4 has mutation on 5 positions
createtestdata<-function(copy, K0 = 15, newK = 5, d = 5, n = 300, 
                         sitect = 5000, tests = 4, postD = 4, folder )
{
  PP = sapply( 1:K0, function(x){ return( NucleoProb( d ) ) } )
  newPP = apply( PP[ , 1:newK ], 2, NewNucleoProb )
  b = n / K0
  rawdata = matrix( NA, nrow = d, ncol = tests * n ) # tests*n?
  for( i in 1:K0 ){
    for( t in 1:tests ){
      rawdata[ , ( b * ( i - 1 ) + 1 ) : ( b * i ) + n * ( t - 1 ) ] = 
        rmultinom( b, sitect, PP[ , i ] )  
    }
  }
  for( i in 1:newK ){
    for( t in postD ){
      rawdata[ , ( b * ( i - 1 ) + 1 ) + n * ( t - 1 ) ] = 
        rmultinom( 1, sitect, newPP[ , i ] )  
    }
  }
  save( rawdata, PP, newPP, tests, sitect, postD, 
        file = paste(folder,"/Data/TestData_",copy,".Rdata",sep=""))
  return(paste("Test", copy," created!",sep=""))
}

#########################################################
#Thresholding the reads by percentage of the total count#
#########################################################
filterdata<-function(Data,threshold)
{
  threshold.list = lapply( 1 : ncol( Data ), function( x ){
    Col = Data[ , x ]; kills = which( Col < (sum( Col ) * threshold) ); Col[ kills ] = 0; return( Col ) } )
  newData = do.call( cbind, threshold.list)
  return( newData )
}


##################################
#Log posterior without normalizer#
##################################
# computes Equation 2 from paper

# supporting function: N is the number of locations,
# L is lambda the expected number of clusters determines a priori, 
# u is the number of unique clusters, and 
# K is some large number which truncates the infinite sum
inf_sum_approx <- function(u, K, N, L) {
  r = seq(u, max(u, K))
  results = r * log(L) - L - lgamma(r - u + 1) - N * log(r)
  return(logSumExp(results, na.rm = TRUE))
}

# main log posterior function: Y is genomic data, 
# C are cluster assignment labels, 
# K is some large number which truncates the infinite sum of the posterior
# a is Dirichlet weight
# L is lambda expected number of clusters
LP_Emma <-
  function(Y, C, K, a, N, L)
    #a is the even weight in the Dirichlet prior, whose concentration=5a.
  {
    J = nrow(Y)
    u = length(unique(C))
    ### Case 1: Y is a vector
    
    if (is.null(dim(Y)) == TRUE) {
      lp = sum(lgamma(Y + (1 / J ^ 2))) - lgamma(sum(Y) + (1 / J)) + inf_sum_approx(u, K, N, L) + u * lgamma(1 / J) - J * u * lgamma(1 / J ^ 2)
      warning("Y IS A VECTOR in LP_Emma")
    }
    
    ### Case 2: Y is a matrix
    M = colSums(Y)
    Vec = rep(0, u)
    for (k in unique(C)) {
      #k in unique C
      I = which(C == k)
      #Case 2.1: I is empty
      if (length(I) == 0) {
        #Vec[k] = J * lgamma(1 / J ^ 2) - lgamma(1 / J)
        stop("ERROR I EMPTY in LP_Emma")
      }
      
      #Case 2.2: I contains one number
      if (length(I) == 1) {
        Vec[which(unique(C)==k)] = sum(lgamma(Y[, I] + a)) - lgamma((M[I]) + J*a)
        #print(Vec)
      }
      #Case 2.3: I is a set of numbers
      else{
        Vec[which(unique(C)==k)] = sum(lgamma(rowSums(Y[, I]) + a)) - lgamma(sum(M[I]) + J*a)
        #print(Vec[k])
      }
    }
    lp = sum(Vec) + inf_sum_approx(u, K, N, L) + u * lgamma(a*J) - J * u * lgamma(a)
    
    return(lp)
  }

###################
#Geweke diagnostic#
###################
geweke.diag2 <- function (x, frac1 = 0.1, frac2 = 0.5) 
{
  if (frac1 < 0 || frac1 > 1) {
    stop("frac1 invalid")
  }
  if (frac2 < 0 || frac2 > 1) {
    stop("frac2 invalid")
  }
  if (frac1 + frac2 > 1) {
    stop("start and end sequences are overlapping")
  }
  if (is.mcmc.list(x)) {
    return(lapply(x, geweke.diag2, frac1, frac2))
  }
  x <- as.mcmc(x)
  xstart <- c(start(x), floor(end(x) - frac2 * (end(x) - start(x))))
  xend <- c(ceiling(start(x) + frac1 * (end(x) - start(x))), 
            end(x))
  y.variance <- y.mean <- vector("list", 2)
  for (i in 1:2) {
    y <- window(x, start = xstart[i], end = xend[i])
    y.matrix <- as.matrix(y)
    y.mean[[i]] <- apply(y.matrix, 2, mean)
    if (identical(sum(y.matrix==y.matrix[1])==length(y.matrix), TRUE)){
      y.variance[[i]] <- 0 
    }else{
      y.variance[[i]] <- spectrum0.ar(y)$spec/niter(y)
    }
  }
  z <- (y.mean[[1]] - y.mean[[2]])/sqrt(y.variance[[1]] + y.variance[[2]])
  out <- list(z = z, frac = c(frac1, frac2))
  class(out) <- "geweke.diag"
  return(out)
}


################################
#2-component Metropolis-Hasting#
################################
split<-function(data,G,s,a,L)
{
  Y=data[,G]
  N=length(G)
  #print(paste0("N = ", N))
  if (N==1) {MC=list(left=G,right=NULL)}
  else if (sum(0==(Y-Y[,1]))==N*(nrow(Y))) {{MC=list(left=G,right=NULL)}}
  else{
    if (N==2){Co=c(1,2)
    }else{Co=kmeans(t(Y),2)[[1]]}
    lp=rep(NA,s)
    CC=matrix(NA,nrow=s,ncol=N)
    #print(Co)
    po=LP_Emma(Y,Co,2,a,N,L) #first occurence
    #print(po)
    t=0
    absZ=10
    while (absZ>2) #Metroplis checked with Geweke diagnostic
    {
      C=Co
      i=ceiling(N*runif(1))
      C[i]=C[i]%%2+1
      p=LP_Emma(Y,C,2,a,N,L)
      d=p-po
      #print(d)
      if (runif(1)<exp(d)){ #Accept the proposal
        Co=C
        po=p
      }
      lp[t%%s+1]=po
      CC[t%%s+1,]=Co			
      if (t%%s+1==s) 
      {
        MClp=mcmc(lp)
        geweke=geweke.diag2(MClp)
        z=geweke$z
        #print(z)
        if (is.na(z)){absZ=0 
        }else if (z=="Inf"){z=10
        }else{absZ=abs(z)}
        c=colSums(CC==1)/s
        MC=list(RunningC=CC, LP=lp, left=G[which(c>=.5)],right=G[which(c<.5)])
      }
      t=t+1
    }
    c=colSums(CC==1)/s
    MC=list(RunningC=CC, LP=lp, left=G[which(c>=.5)],right=G[which(c<.5)])
  }
  return(MC)
}

#############################
#Log ratio of two posteriors#
#############################
LR<-function(Y,C1,C2,i,a,N,L,K)
{
  U = length(unique(C1))                                
  W = length(unique(C2))
  A=which(C1==C1[i])
  B=which(C2==C2[i])
  M=colSums(Y)
  mi=M[i]
  J=nrow(Y)
  lr=rep(0,J)
  for (j in 1:J) 
  {
    y=Y[j,i]
    if (y>0) {lr[j]=-lbeta(sum(Y[j,B])-y+a,y)+lbeta(sum(Y[j,A])-y+a,y)}
  }
  LR=sum(lr)+lbeta(sum(M[B])-mi+a*J,mi)-lbeta(sum(M[A])-mi+a*J,mi) + 
    inf_sum_approx(W, K, N, L) - inf_sum_approx(U, K, N, L) 
  + W*(lgamma(a*J)) - J*W*lgamma(a) - U*(lgamma(a*J)) + U*J*lgamma(a)
  return(LR)
}

#########################################
#Bundle data according to tree end nodes#
#########################################
#add together the number of reads in each group
grpData<-function(Y,Gp)
{
  K=length(Gp)
  Z=matrix(NA,nrow=nrow(Y),ncol=K)
  for (k in 1:K) 
  {
    if (length(Gp[[k]])==1) {Z[,k]=Y[,Gp[[k]]]}
    else {Z[,k]=rowSums(Y[,Gp[[k]]])}	
  }
  return(Z)
}

###################
#Combine groups#
###################
# RunningBC(BCC) is the labels generated in the process of MCMC
# LP is the log posterior generated in the process of MCMC
blockcomb<-function(Z,Gp,s,a,bk,N,L) #bk is the burnin number in terms of multiples of K.
{
  K=length(Gp)
  if (K==1) {print("one group")
  }else {
    lp=rep(NA,s)
    BCC=matrix(NA,nrow=s,ncol=K)
    BC0=seq(1,K,1)
    p0=LP_Emma(Z,BC0,K,a,N,L)
    t=0
    absZ=10
    burn=bk*K
    maxtotal=burn+100000*K #Set an upper bound for number of iterations.
    while (absZ>=2) #Check stationary using Geweke diagnostic.
    {
      BC=BC0
      i=ceiling(K*runif(1))  # randomly choose a coordinate
      Bci=(ceiling((K-1)*runif(1))+BC[i]) %% K # transition kernel f(x,y)=1/(K-1) if x!=y; f(x,x)=0
      if (Bci==0) {BC[i]=K}
      else{BC[i]=Bci}
      d=exp(LR(Z,BC0,BC,i,a,N,L,K=length(Gp)))
      u=runif(1)
      if (t<burn){
        if (u<d) {BC0=BC}
        else {BC0=BC0}	
      }else{
        tb=t-burn
        if (u<d)	
        {
          BC0=BC
          p0=LP_Emma(Z,BC,K,a,N,L)
        }
        else 
        {
          BC0=BC0
          p0=p0
        }
        lp[tb%%s+1]=p0  # store the results of latest s iterations 
        BCC[tb%%s+1,]=BC0			
        
        if (tb%%s+1==s) # diagnose every s steps
        {
          MClp=mcmc(lp) 
          geweke=geweke.diag2(MClp)
          z=geweke$z
          if (is.na(z)){absZ=0 
          }else if (z=="Inf"){absZ=10
          }else{absZ=abs(z)}
        }
      }
      t=t+1
      if (t==maxtotal)
      {
        absZ=0
        print(paste("Stopped after ",maxtotal," iterations",sep=""))
      }
    }   
    MC=list(RunningBC=BCC, LP=lp,t=t,zscore=z)
    return(MC)
  }
}

##############################
#treeGp&BC->DatafinalGp&DataC#
##############################
newGp<-function(treeGp,BC,consoData)
{
  n=dim(consoData)[2]
  C=rep(NA,n)
  K=length(BC)
  nK=length(unique(BC))
  if (nK==K) 
  {
    for (k in 1:K){C[ treeGp[[k]] ]=k} 
    return(list(group=treeGp, index=C))
  }
  else
  {
    G=list(rep(NA,nK))
    for (k in 1:K) {C[treeGp[[k]]]=BC[k]}
    for (j in 1:nK) {G[[j]]=which(C==unique(BC)[j])}
    for (i in 1:nK){C[G[[i]]]=i} 
    return(list(group=G, index=C))
  }	
}

###############################
#runningBC=>GibbsinitialDataCs#
###############################
GibbsC0<-function(runningBC,treeGp,consoData)
{
  if (is.null(dim(runningBC))){return(newGp(treeGp,runningBC,consoData)$index)
  }else{
    GibbC<-function(BC){return(newGp(treeGp,BC,consoData)$index)}
    C0s=apply(runningBC,1,GibbC)
    return(t(C0s))
  }
}

#############################
#rename C
#############################
renameC<-function(C)
{
  newC = rep( NA,length( C ) )
  uniqueC = unique( C )
  for ( i in 1:length(uniqueC) ){ 
    newC[ which( C == uniqueC[i] ) ] = i }
  return( newC )	
}

###############################
#sums of y & m with indicators#
###############################
myfcn<-function(data,C,K)
{
  J=nrow(data)
  csY=colSums(data)
  yInd=matrix(NA,nrow=J,ncol=K)	
  for (k in 1:K)
  {
    Ind=which(C==k)
    l=length(Ind)
    if (l==0){yInd[,k]=rep(0,J)
    }else{ for (j in 1:J) {yInd[j,k]=sum(data[j,Ind])} }
  }
  mInd=colSums(yInd)
  my=list(sMind=mInd,sYind=yInd)	
  return(my)
}

##########################
#Choose a C from runningC#
##########################
convC<-function(CC,lp)
{
  aclp=abs(lp-mean(lp))
  state<-which(aclp==min(aclp))
  if (length(state)==1) {return(CC[state,])}
  else 
  {
    s=state[1]
    return(CC[s,])
  }
}

convC.median<-function(CC,lp)
{
  state<-which(lp==median(lp))
  if (length(state)==1) {return(CC[state,])}
  else 
  {
    s=state[1]
    return(CC[s,])
  }
}

###############
#Gibbs
##############

# Supporting function
# The goal is to calculate the relative probability for 
# assigning a data point i to cluster ki.
# Y is joint genomic data
# C is current cluster assignments
# i is a specific site on genome
# ki is the new group label for position i
# a is dirichlet weight
# L is lambda from prior
# K is some upper bound which truncates infinite sum

LR.Gibbs<-function(Y,C,i,ki,a,N,L,K) 
{
  
  U=length(unique(C))
  A=which(C==C[i])
  
  tempC=C
  tempC[i]=ki
  
  W= length(unique(tempC))
  B=which(tempC==ki)
  
  M=colSums(Y)
  mi=M[i]
  J=nrow(Y)
  lr=rep(0,J)
  
  for (j in 1:J) 
  {
    y=Y[j,i]
    if (y>0) {lr[j]=-lbeta(sum(Y[j,B])-y+a,y)+lbeta(sum(Y[j,A])-y+a,y)}
  }
  LR=sum(lr)+lbeta(sum(M[B])-mi+a*J,mi)-lbeta(sum(M[A])-mi+a*J,mi) + 
    inf_sum_approx(W, K, N, L) - inf_sum_approx(U, K, N, L) 
  + W*(lgamma(a*J)) - J*W*lgamma(a) - U*(lgamma(a*J)) + U*J*lgamma(a)
  return(LR) 
}

GibbsFixScan_EmmaContinue <- function(Data, leftoff, initialC, Scan, a, Record = F, CORE = 2, TempFile, L)
{
  N = length(initialC)
  K = 2*max(initialC) # maximum number of clusters we are willing to consider 
  #change added 07/12/2023 to K
  #numClust = max(unique(initialC))
  
  # Initialize the cluster assignment vector C0 with initialC
  C0 = initialC
  
  # Create placeholders for likelihood probabilities and cluster assignments, if they are to be recorded
  LPs = Cs = NULL
  if( Record == T ){
    LPs=rep(NA,Scan)
    Cs=matrix(NA,nrow=Scan,ncol=N)  
  }
  
  # Initialize the number of completed scans
  Scan0 = 0
  
  # Start a loop until the desired number of scans is reached
  while(Scan0 < Scan){
    
    # Loop through each data point
    for( ii in leftoff:N){
      
      # Update the number of cluster options LR.Gibbs cycles through
      Kseq = c(unique(C0),max(C0)+1)
      #Kseq = seq(1,numClust+1,1) #change here 2/10/24
      #print(paste0("Max Kseq is",max(Kseq)))
      
      # Calculate the log likelihood ratio for each cluster assignment possibility
      #start = proc.time()
      lrs=simplify2array(mclapply(Kseq, LR.Gibbs, Y=Data, C=C0, i=ii, a = a, N=N, L=L, K=K, mc.cores = CORE))
      #end = proc.time() - start
      #print(paste0("time for location ",ii," is ",end[1]+end[5]))
      # Convert log likelihood ratios to probability form and normalize
      rs=exp(lrs-max(lrs)) 
      pp=rs/sum(rs)
      C0[ii]=sample(Kseq, size = 1, prob = pp)
      if( ii %% 10 == 0 ){
        save( ii, Scan0, K, C0, file = TempFile)
      }
      
    }
    # Increment the number of completed scans
    Scan0 = Scan0 + 1
    
    # If recording is enabled, save the likelihood and cluster assignment after each scan
    if( Record == T ){
      LPs[Scan0] = LP_Emma(Data,C0,K,a,N,L)
      Cs[Scan0, ] = C0  
    }
    #save(Scan0, C0, K, file = TempFile)
  }
  
  # Return the final cluster assignment and record of likelihoods and cluster assignments (if recorded)
  MC=list(ScanC=Cs, ScanLP=LPs, Scan=Scan, C = C0)
  return(MC)
}
# Direct Gibbs sampler for our method
# One scan for last step of our prcoedure described in methodology section
# Data is genomic joint data
# initialC is clustering from block MH
# Scan is how many runs of the Gibbs sampler (one for ours)
# a is dirichlet weight
# Record = T saves likelihood and cluster assignment after each scan
GibbsFixScan_Emma <- function(Data, initialC, Scan, a, Record = F, CORE = 2, TempFile, L)
{
  N = length(initialC)
  K = 2*max(initialC) # maximum number of clusters we are willing to consider 
  #change added 07/12/2023 to K
  #numClust = max(unique(initialC))
  
  # Initialize the cluster assignment vector C0 with initialC
  C0 = initialC
  
  # Create placeholders for likelihood probabilities and cluster assignments, if they are to be recorded
  LPs = Cs = NULL
  if( Record == T ){
    LPs=rep(NA,Scan)
    Cs=matrix(NA,nrow=Scan,ncol=N)  
  }
  
  # Initialize the number of completed scans
  Scan0 = 0
  
  # Start a loop until the desired number of scans is reached
  while(Scan0 < Scan){
    
    # Loop through each data point
    for( ii in 1:N){
      
      # Update the number of cluster options LR.Gibbs cycles through
      Kseq = c(unique(C0),max(C0)+1)
      #Kseq = seq(1,numClust+1,1) #change here 2/10/24
      #print(paste0("Max Kseq is",max(Kseq)))
      
      # Calculate the log likelihood ratio for each cluster assignment possibility
      #start = proc.time()
      lrs=simplify2array(mclapply(Kseq, LR.Gibbs, Y=Data, C=C0, i=ii, a = a, N=N, L=L, K=K, mc.cores = CORE))
      #end = proc.time() - start
      #print(paste0("time for location ",ii," is ",end[1]+end[5]))
      # Convert log likelihood ratios to probability form and normalize
      rs=exp(lrs-max(lrs)) 
      pp=rs/sum(rs)
      C0[ii]=sample(Kseq, size = 1, prob = pp)
      if( ii %% 10 == 0 ){
        save( ii, Scan0, K, C0, file = TempFile)
      }
      
    }
    # Increment the number of completed scans
    Scan0 = Scan0 + 1
    
    # If recording is enabled, save the likelihood and cluster assignment after each scan
    if( Record == T ){
      LPs[Scan0] = LP_Emma(Data,C0,K,a,N,L)
      Cs[Scan0, ] = C0  
    }
    #save(Scan0, C0, K, file = TempFile)
  }
  
  # Return the final cluster assignment and record of likelihoods and cluster assignments (if recorded)
  MC=list(ScanC=Cs, ScanLP=LPs, Scan=Scan, C = C0)
  return(MC)
}

# Direct Gibbs comparison. initialC is kmeans initialization. same function as GibbsFixScan_Emma but includes burn and s and multiple Scans
Direct_Gibbs_Sampler <- function(Data, initialC, Scan, a, CORE, TempFile, L, burn, s){
  N = length(initialC)
  K = 2*max(initialC) # maximum number of clusters we are willing to consider (infinite sum truncation)
  # Initialize the cluster assignment vector C0 with initialC
  C0 = initialC
  # Create placeholders for likelihood probabilities and cluster assignments
  #if( Record == T ){
  LPs=rep(NA,Scan)
  Cs=matrix(NA,nrow=Scan,ncol=N)
  #}
  absZ_storage = list()
  absZ = 10
  # Initialize the number of completed scans
  Scan0 = 0
  # Start a loop until the desired number of scans is reached or convergence
  while(Scan0 < Scan && absZ > 2){
    #START = proc.time()
    # Loop through each data point
    for( ii in 1:N){
      # Update the number of cluster options LR.Gibbs cycles through
      Kseq = c(unique(C0),max(C0)+1)
      # Calculate the log likelihood ratio for each cluster assignment possibility
      lrs=simplify2array(mclapply(Kseq, LR.Gibbs, Y=Data, C=C0, i=ii, a = a, N = N,L= L,K= K,mc.cores = CORE))
      # Convert log likelihood ratios to probability form and normalize
      rs=exp(lrs-max(lrs)) 
      pp=rs/sum(rs)
      C0[ii]=sample(Kseq, size = 1, prob = pp)
    }
    # Increment the number of completed scans
    Scan0 = Scan0 + 1
    #END = proc.time() - START
    #print(END)
    # If recording is enabled, save the likelihood and cluster assignment after each scan
    #if( Record == T ){
    LPs[Scan0] = LP_Emma(Data,C0,K,a,N,L)
    Cs[Scan0, ] = C0
    #}
    save(Scan0, C0, K, file = TempFile)
    # Test for convergence 
    if ((Scan0-burn) %% s == 0 && Scan0 > burn) # diagnose every s steps
    {
      MC=mcmc(LPs[1:Scan0])
      gew = geweke.diag(MC)
      z=gew$z
      absZ = abs(z)
      absZ = sapply(absZ, function(x) {
        if (x == "Inf") {
          return(10)
        } else {
          return(x)
        }
      })
      absZ_storage[[ceiling((Scan0 - burn) / s)]] = absZ
    }
  }
  # Return the final cluster assignment and record of likelihoods and cluster assignments (if recorded)
  MC=list(ScanC=Cs, ScanLP=LPs, Scan=Scan, C = C0, Z = absZ_storage)
  return(MC)
}



###########################
#consoDataC=>originalDataC#
###########################
CtoC<-function(consoC,g0,originalData)
{
  J=nrow(originalData)
  C=rep(NA,ncol(originalData))
  gl = lapply( g0, length )
  gl = do.call( rbind, gl )
  glJ = which( gl[ 1:J ] > 0 )
  jj = length( glJ )
  if( jj > 0 ){ 
    for (i in 1 : jj ) {
      j = glJ[ i ]
      C[g0[[j]]]=consoC[i]    
    }
  }
  s=seq(1,length(g0[[J+1]]))
  C[g0[[J+1]][s]]=consoC[s+jj]  
  return(C)
}

#######################
#Gather Gibbs results#
#######################
# function for direct gibbs only
Gibbsresult3<-function(C0s,name,rawdata,GibbsMC=T,G0)
{
  if( GibbsMC == T ){
    cputs = maxt = NULL
    #GibbsFiles = Sys.glob(paste(name,"/FixGibbsDone/FixGibbsMC*.Rdata",sep="")) 
    GibbsFiles = Sys.glob(paste(name,"/GibbsOnly*.Rdata",sep=""))
    # Preallocate a matrix for 100 iterations per chain
    totalRows <- length(GibbsFiles) * 100
    finalCs <- matrix(NA, nrow = totalRows, ncol = ncol(rawdata))
    
    rowIndex <- 1
    chainIndex <- 1
    
    for (filename in GibbsFiles) {
      load(filename)  
      nonNA_rows <- which(apply(GibbsMC$ScanC, 1, function(x) all(!is.na(x))))
      nonNA_matrix <- GibbsMC$ScanC[nonNA_rows, , drop = FALSE]
      niter <- nrow(nonNA_matrix)
      
      # Select 100 evenly spaced indices (thinned iterations)
      indices <- round(seq(1, niter, length.out = 100))
      
      for (i in seq_along(indices)) {
        this_c <- nonNA_matrix[indices[i], ]
        finalCs[rowIndex, ] <- CtoC(renameC(this_c), g0 = G0, originalData = rawdata)
        rowIndex <- rowIndex + 1
      }
    
      cputs[chainIndex] <- cput
      chainIndex <- chainIndex + 1
    }
        maxt <- max(cputs)
  }else{
    finalCs = t( apply( C0s, 1, CtoC, g0=G0, originalData=rawdata ) )
    cputs = maxt = NULL
  }  
  Gibbs=list(FinalCs=finalCs, Runningtimes=cputs, maxCPUtime=maxt)
  return(Gibbs)
}

Gibbsresult2<-function(C0s,name,rawdata,GibbsMC=T,G0) 
{
  if( GibbsMC == T ){
    cputs = maxt = NULL
    GibbsFiles = Sys.glob(paste(name,"/FixGibbsDone/FixGibbsMC*.Rdata",sep="")) 
    GibbsFiles = Sys.glob(paste(name,"/FixGibbsMC*.Rdata",sep=""))
    
    finalCs = matrix(NA,nrow=length(GibbsFiles),ncol=ncol(rawdata))
    i = 1
    for (filename in GibbsFiles)
    {
      load(filename)
      finalCs[i,]=CtoC(renameC(GibbsMC$C), g0=G0, originalData=rawdata)
      cputs[i]=cput
      i = i + 1
    }
    maxt=max(cputs)
  }else{
    finalCs = t( apply( C0s, 1, CtoC, g0=G0, originalData=rawdata ) )
    cputs = maxt = NULL
  }  
  Gibbs=list(FinalCs=finalCs, Runningtimes=cputs, maxCPUtime=maxt)
  return(Gibbs)
}


################################
#Transformed Hellinger distance#
#ln(1-ln(1-H^2))               #
################################
HtMatrix<-function(data,C,a)
{
  K=length(unique(C))
  if (K==1)
  {
    return(paste("Single group, Ht matrix is zero!",sep=""))
  }else{
    Ht=diag(K,x=0)
    indsum=myfcn(data,C,K)
    J=nrow(data)
    sMind=indsum[[1]]
    sYind=indsum[[2]]
    for (i in 1:(K-1))
    {
      for (j in (i+1):K)
      {
        L1j=sYind[,i]+a
        L2j=sYind[,j]+a
        L1=sMind[i]+a*J
        L2=sMind[j]+a*J
        logfrac=sum(lgamma((L1j+L2j)/2))-lgamma((L1+L2)/2)-
          .5*(sum(lgamma(L1j)+lgamma(L2j))-lgamma(L1)-lgamma(L2))
        fH=log(1-logfrac)
        Ht[i,j]=fH
        Ht[j,i]=fH
      }
    }
  } 
  return(Ht)
}

HtDiv<-function(data,runningC,a) #Three time points pairwise Ht distances
{ 
  L=nrow(runningC)
  n=ncol(data)/3
  Ht12=matrix(NA,nrow=L,ncol=n)
  Ht13=matrix(NA,nrow=L,ncol=n)
  Ht23=matrix(NA,nrow=L,ncol=n)
  for (l in 1:L)
  {
    Cl=runningC[l,]
    Htl=HtMatrix(data,Cl,a)
    for (i in 1:n)
    {
      Ht12[l,i]=Htl[Cl[i],Cl[i+n]]
      Ht13[l,i]=Htl[Cl[i],Cl[i+n*2]]
      Ht23[l,i]=Htl[Cl[i+n],Cl[i+n*2]]
    }
  }
  return(list(Ht12,Ht13,Ht23))
}

HtList<-function(data,runningC,a){ 
  L=nrow(runningC)
  Htls=lapply(1:L, function(x){return(HtMatrix(data,runningC[x,],a))} )
  return(Htls)
}

pairHt<-function(Htls,t1,t2,runningC,n)
{
  L=nrow(runningC)
  if (t1==t2){Ht=matrix(0,nrow=L,ncol=n)}
  else{
    Ht=matrix(NA,nrow=L,ncol=n)
    subCs=runningC[,c(((t1-1)*n+1):(t1*n),((t2-1)*n+1):(t2*n))]  
    for (l in 1:L)
    {
      Htl=Htls[[l]]
      subC=subCs[l,]
      Ht[l,]=sapply(1:n,function(x){return(Htl[subC[x],subC[x+n]])})
    }
  }
  return(Ht)
}



#################################
#Extract inference result
##################################

TestResult <- function(TestFiles, mut.true = c(1, 21, 41, 61, 81)){
  # Initialize a list to hold the filenames for each outcome
  filenames <- list(TP = vector(), FN = vector(), FP = vector(), FPN = vector())
  
  for (file in TestFiles){
    load(file)
    freaks = as.numeric(result$Substitution)
    
    if(identical(freaks, mut.true)){
      filenames$TP <- c(filenames$TP, file)
    } else {
      diff.length = length(setdiff(freaks, mut.true))
      diff.length2 = length(setdiff(mut.true, freaks))
      
      if(diff.length == 0){
        filenames$FN <- c(filenames$FN, file)
      } else if(diff.length2 == 0){
        filenames$FP <- c(filenames$FP, file)
      } else {
        filenames$FPN <- c(filenames$FPN, file)
      }
    }
  }
  
  return(filenames)
}

TestTime <- function( TimeFiles ){
  testtime = sapply( TimeFiles, function(x){
    load(x); return(CPUtime)
  })
  return(testtime)
}

############################################
#functions for Dirichlet Process comparison#
############################################

# Helper functions
K <- function(y, theta) {
  # Calculate the likelihood of data given cluster parameter theta
  lh <- dmultinom(x = y, prob = theta) 
  return(lh)
}

P0 <- function(J) {
  # Generate initial parameter values from the prior distribution
  prior <- rdirichlet(1, alpha = rep(1/J^2, J))
  return(prior)
}

relabel_clusters <- function(assignments) {
  unique_clusters <- unique(assignments)
  new_labels <- 1:length(unique_clusters)
  relabeled_assignments <- as.integer(factor(assignments, levels = unique_clusters, labels = new_labels))
  return(relabeled_assignments)
}

# Main Gibbs sampling function
gibbs_sampler <- function(y, alpha, max_clusters, s, burn, name, initialC, n_iterations) {
  # Initialize variables
  #n_iterations = burn+30000
  n <- ncol(y) # number of data points
  J <- nrow(y) # number of categories within data points
  theta_star <- matrix(NA, nrow = J, ncol = max_clusters*n_iterations)
  
  cluster_assignments = initialC # kmeans initialization
  absZ = 10
  iter = 1
  
  theta_star_storage = list()
  cluster_assignments_storage <- matrix(0, nrow = n_iterations, ncol = n)
  num_clusters_storage = rep(0,n_iterations)
  absZ_storage = list()
  
  for (c in unique(initialC)){
    current_cluster_data <- y[, c == cluster_assignments ]
    if (is.null(nrow(current_cluster_data))){
      new_theta <- rdirichlet(1, alpha = (current_cluster_data + rep(1/J^2, J)))
    } else {
      new_theta <- rdirichlet(1, alpha = (rowSums(current_cluster_data) + rep(1/J^2, J)))
    }
    theta_star[, c] <- new_theta
    
  }
  
  # Gibbs sampling loop
  while (any(absZ > 2) && iter <= n_iterations){
    for (subject in 1:n) {  
      
      # Calculate the probabilities of subject 'i' belonging to each existing cluster
      p_cluster <- sapply(unique(cluster_assignments), function(c) {
        K(y[, subject], theta_star[, c])
      })
      
      # Probability of forming a new cluster
      log_base_prior <- sum(lgamma(rep(1/J^2, J))) - lgamma(sum(rep(1/J^2, J)))
      log_data_prob <- sum(lgamma(rep(1/J^2, J) + y[, subject])) - lgamma(sum(rep(1/J^2, J) + y[, subject]))
      log_multinom_coeff <- sum(lfactorial(y[, subject]))
      # Combine all together in log space and return the exponential of the log probability
      log_p_new_cluster <- log(alpha) + log_data_prob + lfactorial(sum(y[, subject])) - log_base_prior - log_multinom_coeff
      p_new_cluster <- exp(log_p_new_cluster)
      
      # Normalize and sample new cluster assignment for the subject
      p_all <- c(p_cluster, p_new_cluster)
      p_all <- p_all / sum(p_all)
      new_assignment <- sample(c(unique(cluster_assignments),max(unique(cluster_assignments))+1), 1, prob = p_all)
      
      
      # Check if a new cluster was formed and update if necessary
      if (new_assignment > (max(cluster_assignments))) {
        cluster_assignments[subject] = new_assignment
        current_cluster_data <- y[, new_assignment == cluster_assignments ]
        if (is.null(nrow(current_cluster_data))){
          new_theta <- rdirichlet(1, alpha = (current_cluster_data + rep(1/J^2, J)))
        } else {
          new_theta <- rdirichlet(1, alpha = (rowSums(current_cluster_data) + rep(1/J^2, J)))
        }
        theta_star[, new_assignment] <- new_theta
      }
      
      # Update subject's cluster assignment if not done already
      cluster_assignments[subject] = new_assignment
      
    }
    
    # Update cluster parameters
    for (c in unique(cluster_assignments)) {
      # Select data points in the current cluster
      current_cluster_data <- y[ , cluster_assignments == c]
      # case 1: current cluster data is vector
      if (is.null(nrow(current_cluster_data))){
        new_theta <- rdirichlet(1, alpha = (current_cluster_data + rep(1/J^2,J)))
      }
      # case 2: current cluster data is matrix
      else{
        new_theta <- rdirichlet(1, alpha = (rowSums(current_cluster_data) + rep(1/J^2,J)))
      }
      theta_star[,c] <- new_theta
    }
    
    # Store 
    theta_star_select <- theta_star[,unique(cluster_assignments)]
    
    theta_star_storage[[iter]] <- list(theta_star_select)
    
    num_clusters_storage[iter] <- length(unique(cluster_assignments))
    
    cluster_assignments_storage[iter, ] <- relabel_clusters(cluster_assignments)
    
    
    # Test for convergence 
    if ((iter-burn) %% s == 0 && iter > burn) # diagnose every s steps
    {
      
      theta_star_clean = theta_star[,unique(cluster_assignments)]
      MC=mcmc(theta_star_clean) 
      gew = geweke.diag(MC)
      z=gew$z
      absZ = abs(z)
      absZ = sapply(absZ, function(x) {
        if (x == "Inf") {
          return(10) 
        } else {
          return(x)  
        }
      })
      absZ_storage[[ceiling((iter - burn) / s)]] = absZ
    }
    
    
    iter = iter + 1
    
    # Optional: print progress every s steps
    if ((iter-burn)%%s==0 && iter > burn) {
      cat("Iteration", iter, "- number of clusters:", length(unique(cluster_assignments)), "\n")
    }
    
  }
  
  if (iter == n_iterations)
  {
    print(paste("Stopped after ",n_iterations," iterations",sep=""))
  }
  
  # Return a list with all the stored data
  return(list(
    theta_star = theta_star_storage[[length(theta_star_storage)]], 
    clusters = cluster_assignments_storage[(nrow(cluster_assignments_storage)-(s-1)):nrow(cluster_assignments_storage),],
    num_clusters = num_clusters_storage[(length(num_clusters_storage)-(s-1)):length(num_clusters_storage)], 
    Z = absZ_storage
    
    
  ))
}


NonParamresult2 <- function(name, rawdata, G0) {
  # Find all IMM cluster files
  NonParamFiles <- Sys.glob(paste(name, "/IMMclusters*.Rdata", sep = ""))
  
  # Preallocate a matrix for 100 iterations per chain
  totalRows <- length(NonParamFiles) * 100
  finalCs <- matrix(NA, nrow = totalRows, ncol = ncol(rawdata))
  
  # Initialize a vector to store CPU times from each chain
  cputs <- numeric(length(NonParamFiles))
  
  rowIndex <- 1
  chainIndex <- 1
  
  for (filename in NonParamFiles) {
    load(filename)  # This loads variables including 'results', 'cput', and 'G0'
    
    # Determine the number of iterations in this chain
    niter <- nrow(results$clusters)
    # Select 100 evenly spaced indices (thinned iterations)
    indices <- round(seq(1, niter, length.out = 100))
    
    for (i in seq_along(indices)) {
      this_c <- results$clusters[indices[i], ]
      finalCs[rowIndex, ] <- CtoC(renameC(this_c), g0 = G0, originalData = rawdata)
      rowIndex <- rowIndex + 1
    }
    
    # Record the CPU time for this chain
    cputs[chainIndex] <- cput
    chainIndex <- chainIndex + 1
  }
  
  maxt <- max(cputs)
  
  NonParam <- list(FinalCs = finalCs, Runningtimes = cputs, maxCPUtime = maxt)
  return(NonParam)
}

#############################################
######### New Post Processing code ##########
##############################################

# Takes in assignment labels from all s chains of the Gibbs sampler and outputs the hellinger 
# distances for every site and time point (called Hts) for each chain. Hts is a list
# where Hts[[1]] decribes the t1t2 comparisons, the Hts[[2]] describes the remaining control
# comparisons and Hts[[3]] describes the treatment comparisons. Hts[[2]][[1]] is all left over control comparisons involving t1 while 
# Hts[[2]][[2]] is all control compairons involving t2, etc. Hts[[2]][[1]][[1]] contains all chains from the t1t3 comparison. if there are more then 3
# timepoints there will be more sublists i.e. Hts[[2]][[1]][[2]] would be t1t4 comparison.

Ht.Emma <- function(GibbsCs, Rawdata, a, PostD, n, FolderName, Control = T, mHtsFile = "mHts.Rdata") {
  Finalcs <- t(apply(GibbsCs, 1, renameC))  
  Htlist <- HtList(Rawdata, Finalcs, a)    
  Hts <- list(t1t2 = NULL, N = NULL, D = NULL)  
  
  # Compute t1t2 values for all chains
  ht1 <- pairHt(Htlist, 1, 2, Finalcs, n)  # pairwise distances for t1 vs t2
  Hts[[1]] <- ht1  # Store all chain-specific values
  
  D.length <- length(PostD)
  
  # Add comparisons for control and treatment
  if (Control == T) {
    for (i in 2:3) {
      for (j in 1:2) {
        Hts[[i]][[j]] <- list()
        for (t in 1:D.length) {
          k <- (i - 2) * D.length + 2 + t
          htijt <- pairHt(Htlist, j, k, Finalcs, n)
          Hts[[i]][[j]][[t]] <- htijt  # Save chain-specific values
        }
      }
    }
  } else {
    Hts[[2]] <- NA
    for (j in 1:2) {
      Hts[[3]][[j]] <- list()
      for (t in 1:D.length) {
        k <- 2 + t
        htijt <- pairHt(Htlist, j, k, Finalcs, n)
        Hts[[3]][[j]][[t]] <- htijt  # Save chain-specific values
      }
    }
  }
  
  # Construct mHts.pool as a matrix
  # Flatten all chain-specific matrices into vectors
  flat_matrices <- c(
    list(Hts[[1]]),  # t1t2 comparisons
    unlist(Hts[[2]], recursive = FALSE),  # Control comparisons
    unlist(Hts[[3]], recursive = FALSE)   # Treatment comparisons
  )
  
  # Combine into a single matrix
  # Each row = a chain; each column = flattened comparison values
  num_chains <- nrow(Finalcs)  # Number of chains
  num_sites <- n  # Number of sites
  num_comparisons <- length(flat_matrices) * num_sites
  mHts.pool.mat <- matrix(NA, nrow = num_chains, ncol = num_comparisons)
  
  for (chain in 1:num_chains) {
    mHts.pool.mat[chain, ] <- unlist(lapply(flat_matrices, function(mat) {
      mat[chain, ]
    }))
  }
  
  # Save results
  save(Finalcs, Hts, mHts.pool.mat, file = paste(FolderName, mHtsFile, sep = ""))
}


# Computes # Ht_{s,i}^N for each chain and site and # Ht_{s,i}^D. Saved as max_N_t1t2, min_D.
# Requires InfFile which is /mHtsFixOneEmma.Rdata, the results from function HtEmma.
compute_site_chain_extremes <- function(InfFile) {
  load(InfFile)
  
  # Ensure Hts is available
  if (is.null(Hts)) stop("Hts is not available in the loaded file.")
  
  # Extract number of chains and sites
  num_chains <- nrow(Hts[['t1t2']]) 
  num_sites <- ncol(Hts[['t1t2']])
  
  # Initialize matrices to store results
  max_N_t1t2 <- matrix(NA, nrow = num_chains, ncol = num_sites)  # Ht_{s,i}^N
  min_D <- matrix(NA, nrow = num_chains, ncol = num_sites)       # Ht_{s,i}^D
  
  # combine all time points, locations, 
  for (chain in 1:num_chains) {
    for (site in 1:num_sites) {
      # Collect all values for N
      N_values <- c(
        Hts[['t1t2']][chain, site],  
        unlist(lapply(Hts[['N']], function(subset) {  
          unlist(lapply(subset, function(time_point_matrix) {  
            time_point_matrix[chain, site]
          }))
        }))
      )
      
      # Compute max for t1t2 and N
      max_N_t1t2[chain, site] <- max(N_values, na.rm = TRUE)
      
      # Collect all values for D
      D_values <- unlist(lapply(Hts[['D']], function(subset) {  
        unlist(lapply(subset, function(time_point_matrix) {         
          time_point_matrix[chain, site]
        }))
      }))
      
      # Compute min for D
      min_D[chain, site] <- min(D_values, na.rm = TRUE)
    }
  }
  
  # Return the computed max and min matrices
  return(list(max_N_t1t2 = max_N_t1t2, min_D = min_D))
}


sub.sites.Emma.new <- function(quantiles_control, quantiles_tmt, Cuts) {
  
  # Helper function to identify sites for given quantile vectors from control and treatment
  identify_sites <- function(q_control, q_tmt) {
    sapply(Cuts, function(Cut) {
      D.Site <- which(q_tmt > Cut)
      DN.Site <- intersect(which(q_tmt < q_control), which(q_control > Cut))
      S <- setdiff(D.Site, DN.Site)
      list(Substitution = S, potential = D.Site, noise = DN.Site)
    }, simplify = FALSE)
  }
  
  # For q = 0.25: treatment uses "25%" and control uses "75%"
  q25 <- identify_sites(quantiles_control[, "75%"], quantiles_tmt[, "25%"])
  
  # For q = 0.5: both treatment and control use "50%"
  qmed <- identify_sites(quantiles_control[, "50%"], quantiles_tmt[, "50%"])
  
  # For q = 0.75: treatment uses "75%" and control uses "25%"
  q75 <- identify_sites(quantiles_control[, "25%"], quantiles_tmt[, "75%"])
  
  # Return the results with improved naming
  list(
    q25 = q25,
    qmed = qmed,
    q75 = q75,
    cutoff = Cuts
  )
}

# Calculate lengths of potential, noise, and signal sets for each quantile
calculate_set_lengths <- function(sites) {
  sapply(sites, function(site) {
    c(length(site$Substitution), length(site$potential), length(site$noise))
  })
}

# very similar to orignal Jenny function (freak.result1) which computes the elbow of the noise set. We have 
# Ht for each chain now so we must pick a summary statistic of the Ht to use to compute the elbow of the
# noise set. Must pick from Q1, Q2, or Q3 of summary statistics as those are the ones computed above. 
noise.elbow.Emma.new <- function( InfFile, ResultFile, SaveFile, delta = 3, Control = T, alpha = 0){
  load( InfFile ); load( ResultFile ) #Provide mHts.pool, Pairs, PostD, noise.test (lengths of pot, sig, noise sets), quantiles_control, quantiles_tmt
  noise.test = noise.test.mat[['Q1']]
  noise = noise.test[ 3, 1:ncol(noise.test) ]
  
  noise.unique.length = length(unique(noise))
  if( noise.unique.length < delta * 2 ){ 
    delta = 1 
  }
  if( noise.unique.length == 1 ){
    cutoff.s = cutoff.s0 = 1; noi.l = length(noise); noi.ldiff = 0
  }else{
    noi.l = sapply( 0:max(noise), function(x){return(sum(noise == x))} )
    delta2 = delta*2-1
    noi.ldiff = sapply( 1:( length(noi.l) - delta2 ), function(x){
      del = sum(noi.l[x + delta:delta2]) - sum(noi.l[x + 0:(delta - 1)])
      return(del)
    })
    cutoff.s0 = max(which(noi.ldiff == min(noi.ldiff)))
    cutoff.s = cutoff.s0 + delta - 1
  }    
  if(cutoff.s0 == 1){
    delta0 = delta
    buffer = mean( noi.l[1:delta0] ) * alpha
    while( noi.l[delta0] < buffer){
      cutoff.s = cutoff.s - 1 
      delta0 = delta0 - 1
      buffer = mean( noi.l[1:cutoff.s] ) * alpha
    }     
  }else{ buffer = mean( noi.l[(cutoff.s - delta + 1):cutoff.s] ) * alpha }
  
  # Start with cutoff.s - 3 and iterate upward, staying within the delta group
  found = FALSE
  for (adjustment in 1:delta) {
    candidate = cutoff.s - adjustment  
    if (candidate < cutoff.s - delta) {
      warning("Reached the delta group lower bound without finding a match.")
      break
    }
    if (any(noise.test[3,] == candidate)) {
      d = Cuts[min(which(noise.test[3,] == candidate)) + buffer]
      found = TRUE
      break
    }
  }

  # Handle case where no match is found within the delta group
  if (!found) {
    stop("No valid match found within the delta group for cutoff.s adjustments in noise.test[3,].")
  }
  
  result_all = sub.sites.Emma.new(quantiles_control, quantiles_tmt, d)
  result <- c(result_all[['q25']][[1]], list(cutoff = d))
  save(delta, buffer, noi.l, noi.ldiff, cutoff.s, d, result, alpha, file = SaveFile)
  return(result)
}

################################################
############### New Plot Code ##################
###############################################

plot_hts <- function(savefile, resultfile, savefile2 = NULL, savefile3 = NULL) {
  # Load the data
  load(savefile)
  load(resultfile)
  
  # Extract threshold and substitution data
  threshold_d <- d
  max_N_t1t2 <- result_stats$max_N_t1t2
  min_D <- result_stats$min_D
  sub_sites <- result$Substitution
  
  # If additional files are provided, load them and extract threshold data
  if (!is.null(savefile2)) {
    load(savefile2)
    threshold_d2 <- d
  }
  if (!is.null(savefile3)) {
    load(savefile3)
    threshold_d3 <- d
  }
  
  # Compute Q1, Q2 (median), and Q3 for each site for both matrices
  q1_N <- apply(max_N_t1t2, 2, quantile, probs = 0.25, na.rm = TRUE)
  q2_N <- apply(max_N_t1t2, 2, median, na.rm = TRUE)
  q3_N <- apply(max_N_t1t2, 2, quantile, probs = 0.75, na.rm = TRUE)
  
  q1_D <- apply(min_D, 2, quantile, probs = 0.25, na.rm = TRUE)
  q2_D <- apply(min_D, 2, median, na.rm = TRUE)
  q3_D <- apply(min_D, 2, quantile, probs = 0.75, na.rm = TRUE)
  
  # Define the site indices and determine the y-axis limit
  sites <- 1:ncol(max_N_t1t2)
  ylim_max <- max(c(q3_N, q3_D), na.rm = TRUE)
  extended_ylim <- c(0, ylim_max * 1.3)  # Extend the y-axis by 30%
  
  # Create x-axis ticks starting at 0 and incrementing by 50
  x_ticks <- seq(0, max(sites), by = 50)
  
  # Define custom colors for plotting
  col_N <- "#1E88E5"      # For A_i^T(q)
  col_D <- "#D81B60"      # For A_i^C(q)
  col_signal <- "#FFC107" # For Signal
  
  # Create the main plot with an extended y-axis
  plot(
    sites, q2_N,
    type = "n",
    ylim = extended_ylim,
    xlim = c(0, max(sites)),
    xlab = "Locations on the Genome",
    ylab = "Transformed Hellinger Distance (A)",
    main = "Summary Statistics Plot",
    xaxt = "n",
    font.lab = 2,
    cex.lab = 1.5,
    cex.main = 1.5
  )
  
  # Customize the x-axis (horizontal labels)
  axis(1, at = x_ticks, labels = x_ticks, las = 1, cex.axis = 0.8)
  
  # Add slight horizontal jitter to avoid overplotting
  jitter_offset_N <- -0.09 + runif(length(sites), -0.03, 0.03)
  jitter_offset_D <-  0.09 + runif(length(sites), -0.03, 0.03)
  sites_N <- sites + jitter_offset_N
  sites_D <- sites + jitter_offset_D
  
  # Plot for A_i^C(q) using jittered x coordinates 
  segments(sites_N, q1_N, sites_N, q3_N, col = col_N, lwd = 0.5)
  points(sites_N, q2_N, pch = 19, col = col_N, cex = 0.4)
  
  # Plot for A_i^T(q) using jittered x coordinates
  segments(sites_D, q1_D, sites_D, q3_D, col = col_D, lwd = 0.5)
  points(sites_D, q2_D, pch = 19, col = col_D, cex = 0.4)
  
  # Highlight substitution sites with open circles for signal.
  points(sites_D[sub_sites], q2_D[sub_sites], col = col_signal, pch = 1, cex = 1.5, lwd = 2)
  
  # Draw horizontal threshold lines and label them with d
  if (is.null(savefile2) && is.null(savefile3)) {
    abline(h = threshold_d, col = "black", lty = 2, lwd = 1.5)
    text(0, threshold_d, expression(d), pos = 3, offset = 0.1, cex = 1.5)
  } else if (!is.null(savefile2) && is.null(savefile3)) {
    abline(h = threshold_d, col = "black", lty = 2, lwd = 1.5)
    abline(h = threshold_d2, col = "black", lty = 3, lwd = 1.5)
    text(0, threshold_d2, expression(d), pos = 3, offset = 0.1, cex = 1.5)
  } else if (!is.null(savefile2) && !is.null(savefile3)) {
    abline(h = threshold_d, col = "black", lty = 2, lwd = 1.5)
    abline(h = threshold_d2, col = "black", lty = 3, lwd = 1.5)
    abline(h = threshold_d3, col = "black", lty = 4, lwd = 1.5)
    text(0, threshold_d3, expression(d), pos = 3, offset = 0.1, cex = 1.5)
  }
  

  # Place a horizontal legend at the top middle inside the plot.
  # Set y_top as 95% of the extended y-axis maximum so it is inside the plot.
  x_center <- mean(par("usr")[1:2])
  y_top <- extended_ylim[2] * 0.95
  
  legend(
    x = x_center,
    y = y_top,
    legend = c(expression(A[i]^C * (q)), expression(A[i]^T * (q)), "Signal"),
    col = c(col_N, col_D, col_signal),
    pch = c(19, 19, 1),
    lty = c(1, 1, NA),
    horiz = TRUE,
    xjust = 0.5,
    yjust = .53,
    cex = 1.3,
    bty = "o"
  )
}

PlotThresholdEmma <- function(InfFile, ResultFile, SaveFile, SaveFile2 = NULL, SaveFile3 = NULL, 
                              Main, quantile_label1 = "Q2", quantile_label2 = NULL, quantile_label3 = NULL, 
                              save_path = NULL, x_axis = NULL) {
  # Load data files
  load(SaveFile)
  load(InfFile)
  load(ResultFile)
  d1 <- d
  print(d1)
  
  # Load thresholds from SaveFile2 and SaveFile3 if provided
  if (!is.null(SaveFile2)) {
    load(SaveFile2)
    d2 <- d
    print(d2)
  }
  if (!is.null(SaveFile3)) {
    load(SaveFile3)
    d3 <- d
  }
  
  # Subsample every threshold (by = 1, so no skipping)
  subsample_idx <- seq(1, length(Cuts), by = 1)
  subsampled_Cuts <- Cuts[subsample_idx]
  reversed_Cuts <- rev(subsampled_Cuts)  # Reverse for largest to smallest
  
  # Subsample Q2 (median) values and Q1 values for each category
  noise_signal    <- rev(noise.test.mat[["Qmed"]][1, ][subsample_idx])  # Q2 for Signal
  noise_potential <- rev(noise.test.mat[["Qmed"]][2, ][subsample_idx])  # Q2 for Potential
  noise_noise     <- rev(noise.test.mat[["Qmed"]][3, ][subsample_idx])  # Q2 for Noise
  
  signal_q1       <- rev(noise.test.mat[["Q1"]][1, ][subsample_idx])    # Q1 for Signal
  potential_q1    <- rev(noise.test.mat[["Q1"]][2, ][subsample_idx])    # Q1 for Potential
  noise_q1        <- rev(noise.test.mat[["Q1"]][3, ][subsample_idx])    # Q1 for Noise
  
  # Define color palettes
  # Q2 (main colors)
  col_signal_new     <- "#FFC107"  # Signal (Q2)
  col_potential_new  <- "#D81B60"  # Potential (Q2)
  col_noise_new      <- "#1E88E5"  # Noise (Q2)
  
  # Q1 (lighter hues)
  col_signal_q1     <- "#FFECB3"   # Signal (Q1)
  col_potential_q1  <- "#F8BBD0"   # Potential (Q1)
  col_noise_q1      <- "#90CAF9"   # Noise (Q1)
  
  # Set up plotting parameters
  par(family = "sans", mar = c(5, 5, 4, 2) + 0.1)
  
  # Open a PDF device if saving
  if (!is.null(save_path)) {
    pdf(save_path, width = 8, height = 6)
  }
  
  # Determine xlim based on x_axis input; if NULL, use reversed_Cuts
  if (is.null(x_axis)) {
    xlim_range <- c(max(reversed_Cuts), min(reversed_Cuts))
    # Use all indices if no x_axis is specified
    in_range <- rep(TRUE, length(reversed_Cuts))
  } else {
    xlim_range <- c(max(x_axis), min(x_axis))
    # Identify indices where reversed_Cuts falls within the chosen x range
    in_range <- reversed_Cuts <= max(x_axis) & reversed_Cuts >= min(x_axis)
  }
  
  # Set y-limit based on maximum of Q2 values for points within the chosen x range
  ylim_range <- c(0, max(c(noise_signal[in_range],
                           noise_potential[in_range],
                           noise_noise[in_range]), na.rm = TRUE))
  
  # Create an empty plot
  plot(
    reversed_Cuts,
    noise_signal,
    type = "n",
    xlim = xlim_range,
    ylim = ylim_range,
    main = Main,
    xlab = "Threshold d",
    ylab = "Number of Locations",
    font.lab = 2,
    cex.main = 1.5,
    cex.lab = 1.5
  )
  
  # Plot Q2 (median) as filled circles and Q1 as open circles for each threshold point
  n <- length(reversed_Cuts)
  for (i in seq_len(n)) {
    # Q2 data (filled circles)
    points(reversed_Cuts[i], noise_signal[i],    col = col_signal_new,    pch = 16, cex = 0.5)
    points(reversed_Cuts[i], noise_potential[i], col = col_potential_new, pch = 16, cex = 0.5)
    points(reversed_Cuts[i], noise_noise[i],     col = col_noise_new,     pch = 16, cex = 0.5)
    
    # Q1 data (open circles)
    points(reversed_Cuts[i], signal_q1[i],    col = col_signal_q1,    pch = 16, cex = 0.5)
    points(reversed_Cuts[i], potential_q1[i], col = col_potential_q1, pch = 16, cex = 0.5)
    points(reversed_Cuts[i], noise_q1[i],     col = col_noise_q1,     pch = 16, cex = 0.5)
  }
  
  # Plot vertical dashed lines for thresholds
  abline(v = d1, col = "black", lty = 2, lwd = 1.5)
  if (!is.null(SaveFile2)) abline(v = d2, col = "black", lty = 3, lwd = 1.5)
  if (!is.null(SaveFile3)) abline(v = d3, col = "black", lty = 4, lwd = 1.5)
  
  # Create a single legend with 6 entries.
  # For Q2 entries, use filled circles; for Q1 entries, use open circles.
  legend_labels <- c(
    expression(S[1]^{(d*","*q==.5)}),
    expression(S[2]^{(d*","*q==.5)}),
    expression(S[3]^{(d*","*q==.5)}),
    expression(S[1]^{(d*","*q==.25)}),
    expression(S[2]^{(d*","*q==.25)}),
    expression(S[3]^{(d*","*q==.25)})
  )
  
  
  legend_colors <- c(
    col_potential_new, col_noise_new, col_signal_new,
    col_potential_q1, col_noise_q1, col_signal_q1
  )
  
  legend_pch <- c(16, 16, 16, 16, 16, 16)
  
  legend(
    "topleft",  
    legend = legend_labels,      # your 6 expressions
    col    = legend_colors,      # matching colors
    pch    = legend_pch,         # same shape or different
    bty    = "o",                # box around legend
    cex    = 1.5,                # text size
    box.lwd = 2,                 
    box.col = "black",
    x.intersp = 1.5,             # horizontal spacing among legend items
    y.intersp = 1.2,             # vertical spacing among legend items
    ncol = 2                     # 2-column layout
  )
  
  
  if (!is.null(save_path)) dev.off()
}
