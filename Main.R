# Algorithm (preprocess, processing, and post-process)
library(stringr)

############ Step1: Preprocess ##################
conso=consoData(rawdata)  # consoData is in AllfcnsFinal.R
Data=conso[[1]]
G0=conso[[2]]
N = dim(Data)[2]
save(Data,G0,L,J,a,PostD,Pairs,GibbsRun,file=paste(name,"/Data.Rdata",sep=""))


########### Step2: Processsing: Depends on algorithm chosen ###############
dir.create(paste(name,"/Outputs",sep=""))
dir.create(paste(name,"/SLURMouts",sep=""))

## Our algorithm
if(algorithm == "ours"){
  # Hierarchical SCMH
  
  starttree = paste("sbatch --constraint=rhel8 -t 1- --mem=3000 -N 1 -n 1 -J ",Jname," -o ", name, "/SLURMouts/", Jname,
                    ".out.1  --wrap=\"R  CMD BATCH --vanilla --args --copy=", copy,
                    " --name=",name, " --Jname=", Jname, " --current=1 Tree.R ", name, "/Outputs/tree1.Rout \"",sep="" )
  
  print (paste("sbatch command:",starttree))
  print (" ")
  system(starttree)
  
  
  state <- system(paste("squeue -u", user,"-O  jobid,name:40"),intern=TRUE)
  # collapse character vector into on long string
  state <- paste(state,collapse=" ")
  print("printing state")
  print(state)
  isInQueue <- str_detect(state, Jname)
  print (paste("isInQueue=",isInQueue))
  while (isInQueue) {
    Sys.sleep(90)
    state <- system(paste("squeue -u", user, "-O  jobid,name:40"),intern=TRUE)
    state <- paste(state,collapse=" ")
    isInQueue <- str_detect(state, Jname)
    print (paste("whiling away isInQueue=",isInQueue))
  }
  
  
  SplitList=list() #Summarize tree result.
  print (paste("1 spl=",SplitList))
  run=T
  i=0
  endInd=3
  print ("doing while")
  while(run==T){
    i=i+1
    print (paste("i=",i))
    print(paste("Temporary directory:", tempdir()))
    if(1==i){
      load(paste(name,"/Tree_",i,".Rdata", sep=""))
      test=old.result
      SplitList[[i]]<-SplitNode(left=2*i,right=2*i+1,parent=NA,curent=i,LP=test$LP,
                                ind=sort(c(test$left,test$right)))    #SplitNode is in source.R
      if(0!=length(test$left) & 0!=length(test$right)){
        SplitList[[2*i]]<-SplitNode(left=4*i,right=4*i+1,curent=2*i,parent=i,LP=NULL,ind=test$left)
        SplitList[[2*i+1]]<-SplitNode(left=2*(2*i+1),right=2*(2*i+1)+1,curent=2*i+1,parent=i,
                                      LP=NULL,ind=test$right)
      }
      if(0==length(test$left) | 0==length(test$right)){
        noLeft(SplitList[[i]])    #noLeft is in source.R
        noRight(SplitList[[i]])  #noRight is in source.R
      }
    }
    else if(is.Node(SplitList[[i]])){
      load(paste(name,"/Tree_",i,".Rdata", sep=""))
      test=old.result
      updateLP(SplitList[[i]],test$LP)  #updateLP is in source.R
      if(0!=length(test$left) & 0!=length(test$right)){
        SplitList[[2*i]]<-SplitNode(left=4*i,right=4*i+1,curent=2*i,parent=i,LP=NULL,ind=test$left) #Function is in source.R
        SplitList[[2*i+1]]<-SplitNode(left=2*(2*i+1),right=2*(2*i+1)+1,curent=2*i+1,parent=i,
                                      LP=NULL,ind=test$right)
      }else{
        noLeft(SplitList[[i]])
        noRight(SplitList[[i]])
        SplitList[[2*i]]=NULL
        SplitList[[2*i+1]]=NULL
      }
    }else{
      SplitList[[2*i]]=NULL
      SplitList[[2*i+1]]=NULL
    }
    if (i==endInd){
      st=stopSplit(i)
      if(0==st){run=F}
      endInd=2*endInd+1
      SplitList[[2*endInd+1]]=NA
    }
  }
  print ("end while T")
  Gp=groupIndex()
  save(Gp, file=paste(name,"/TreeResult.Rdata",sep=""))
  
  
  print ("Block MH section")
  # Block MH
  Z0=grpData(Data,Gp)
  T0 = proc.time()
  combgrp=blockcomb(Z0,Gp,50000,a,1000,N,L)  
  TN = proc.time() - T0
  cput_block = TN[1]
  save(Z0,combgrp, cput_block,file=paste(name,"/BlockResult.Rdata",sep=""))
  
  print ("Gibbs mod section")
  
  BCs=combgrp[[1]][seq(1,10000,100),] #Thin-out the chain by every 100 iterations.
  C0s=GibbsC0(BCs,Gp,Data) #Assign group label to each sequence position. 100*1175 matrix
  newC0s=t(apply(C0s,1,renameC)) #Rename the labels.
  save(newC0s, file = paste( name, "/FixGibbsResult.Rdata", sep = "" ) )
  
  if (GibbsRun=="T"){ #ours with last direct gibbs step
    numC0=nrow(newC0s)
    
    
    for(C0 in 1:numC0){
      
      subFix=paste("sbatch --constraint=rhel8 -t 1- -N 1 -n ", CORE, " -J FixGibbs.", Jname, "[",C0,"] -o  ", name, "/SLURMouts/FixGibbs%j.out  --wrap=\"R CMD BATCH --vanilla --args --name=", name,
                   " --copy=",copy, " --core=", CORE, " --Ci=",C0," FixGibbsMC.R ", name,"/Outputs/FixGibbs",C0,".out\"", sep="")
      system(subFix)
      Sys.sleep(1)
    }
    
    state <- system(paste("squeue -u", user ,"-O  jobid,name:40"),intern=TRUE)
    # collapse character vector into on long string
    state <- paste(state,collapse=" ")
    print("printing state")
    print(state)
    isInQueue <- str_detect(state, Jname)
    print (paste("isInQueue=",isInQueue))
    while (isInQueue) {
      Sys.sleep(60)
      system(paste("squeue -u", user, "-O  name:40"))
      state <- system(paste("squeue -u", user, "-O  name:40"),intern=TRUE)
      state <- paste(state,collapse=" ")
      isInQueue <- str_detect(state, Jname)
      print (paste("whiling away isInQueue=",isInQueue))
    }

    system(paste("rm ", name, "/FixGibbsMC*temp.Rdata", sep="" ))
    Gibbs=Gibbsresult2(newC0s,name,rawdata, GibbsMC = T, G0)
    save( newC0s, Gibbs, file = paste( name, "/FixGibbsResult.Rdata", sep = "" ) )
  } else{ # ours without direct gibbs step
    Gibbs=Gibbsresult2(newC0s,name,rawdata, GibbsMC = F, G0) 
    save( newC0s, Gibbs, file = paste( name, "/FixGibbsResult.Rdata", sep = "" ) )
  }
} else if(algorithm == "DP"){ ## Dirichlet Process algortihm 
  chains = 3
  
  initialC <- list()
  
  # Perform k-means clustering 3 times
  for (i in 1:chains) {
    for (j in c(5,15,25))
      km.res <- kmeans(t(Data), j, nstart = 1)
    initialC[[i]] <- km.res$cluster
  }
  
  save(initialC, file = paste( name, "/NonParamResult.Rdata", sep = "" ) )
  numC0 = length(initialC)
  
  for(C0 in 1:numC0){
    subFix=paste("sbatch --constraint=rhel8 --mem=3000 -t 1- -N 1 -n ", CORE, " -J NonParamMC.", Jname, "[",C0,"] -o  ", name, "/SLURMouts/NonParamMC%a.out  --wrap=\"R CMD BATCH --vanilla --args --name=", name,
                 " --copy=",copy, " --core=", CORE, " --Ci=",C0," NonParamMC.R ", name,"/Outputs/NonParamMC",C0,".out\"", sep="")
    system(subFix)
    Sys.sleep(1)
  }
  
  state <- system(paste("squeue -u", user, " -O  jobid,name:40"),intern=TRUE)
  # collapse character vector into on long string
  state <- paste(state,collapse=" ")
  print("printing state")
  print(state)
  isInQueue <- str_detect(state, Jname)
  print (paste("isInQueue=",isInQueue))
  while (isInQueue) {
    Sys.sleep(60)
    system(paste("squeue -u", user, "-O  name:40"))
    state <- system(paste("squeue -u", user, "-O  name:40"),intern=TRUE)
    state <- paste(state,collapse=" ")
    isInQueue <- str_detect(state, Jname)
    print (paste("whiling away isInQueue=",isInQueue))
  }


  Gibbs=NonParamresult2(name, rawdata, G0) #Gather Gibbs results #without 2?
  save(initialC, Gibbs, file = paste( name, "/NonParamResult.Rdata", sep = "" ) )
}else if(algorithm == "direct"){ ## Direct Gibbs
  chain = 3
  initialC <- list() 
  
  # Perform k-means clustering 3 times
  for (i in 1:chain) {
    for (j in c(5,15,25))
      km.res <- kmeans(t(Data), j, nstart = 1)
    initialC[[i]] <- km.res$cluster
  }
  
  ClusterData <- t(sapply(initialC, function(x) x))
  newC0s = t(apply(ClusterData, 1, renameC))
  save(newC0s, file=paste(name, "/FixGibbsResult.Rdata", sep=""))
  numC0 = length(initialC)
  
  
  for(C0 in 1:numC0){
    subFix=paste("sbatch --constraint=rhel8 -t 3- -N 1 -n ", CORE, " -J GibbsMC.", Jname, "[",C0,"] -o  ", name, "/SLURMouts/GibbsMC%j.out  --wrap=\"R CMD BATCH --vanilla --args --name=", name,
                 " --copy=",copy, " --core=", CORE, " --Ci=",C0," DirectGibbsMC.R ", name,"/Outputs/GibbsMC",C0,".out\"", sep="")
    system(subFix)
    Sys.sleep(1)
  }
  
  state <- system(paste("squeue -u", user, "-O  jobid,name:40"),intern=TRUE)
  # collapse character vector into on long string
  state <- paste(state,collapse=" ")
  print("printing state")
  print(state)
  isInQueue <- str_detect(state, Jname)
  print (paste("isInQueue=",isInQueue))
  while (isInQueue) {
    Sys.sleep(60)
    system(paste("squeue -u", user, "-O  name:40"))
    state <- system(paste("squeue -u", user, "-O  name:40"),intern=TRUE)
    state <- paste(state,collapse=" ")
    isInQueue <- str_detect(state, Jname)
    print (paste("whiling away isInQueue=",isInQueue))
  }


  Gibbs = Gibbsresult3(newC0s, name, rawdata, GibbsMC=T, G0)
  save(newC0s, Gibbs, file=paste(name, "/FixGibbsResult.Rdata", sep=""))
  
}

########### Step3: post processsing:###############
Ht.Emma(Gibbs$FinalCs, rawdata, a, PostD, n, name, Control = T, mHtsFile = "/mHtsFixOne.Rdata" )
inffile = paste( name, "/mHtsFixOne.Rdata", sep = "" )
load(inffile)


resultfile =  paste( name, "/ResultLocalFixOne.Rdata", sep = "" )
result_stats = compute_site_chain_extremes(inffile)  #Computes Ht_{s,i}^N for each chain and site and Ht_{s,i}^D


# Compute quantiles of Ht(Ni) and Ht(Di)
compute_quantiles <- function(data) {
  t(apply(data, 2, function(x) round(quantile(x, probs = c(0.25, 0.5, 0.75)), 2)))
}

quantiles_control <- compute_quantiles(result_stats$max_N_t1t2)
quantiles_tmt <- compute_quantiles(result_stats$min_D)

Cuts <- seq(0, max(c(result_stats$max_N_t1t2, result_stats$min_D)), by = 0.001) 
# Compute sites for all Cuts at q=.25,.5, and .75
result_sites <- sub.sites.Emma.new(quantiles_control, quantiles_tmt, Cuts)

# Calculate lengths of potential, noise, and signal sets for median and chosen q and 1-q
noise.test.mat <- list(
  Qmed = calculate_set_lengths(result_sites$qmed),
  Q1 = calculate_set_lengths(result_sites$q25),
  Q3 = calculate_set_lengths(result_sites$q75)
)

save( Pairs, PostD, Cuts, noise.test.mat, quantiles_control, quantiles_tmt, result_stats, file = resultfile )
savefile =  paste( name, "/FinalFixOne1.Rdata", sep = "" )

result_final = noise.elbow.Emma.new( inffile, resultfile, savefile, delta = 3, Control = T, alpha = 0) 
save(Gibbs, result_final, file = paste( name, "/InferenceFix1.Rdata",sep=""))


