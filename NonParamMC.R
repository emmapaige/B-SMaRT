source("ModifiedAllfcns.R")
source("LibraryCmdline.R")

Ci<-cmdline.numeric("Ci")
name<-cmdline.strings("name")
copy<-cmdline.numeric("copy")
CORE<-cmdline.numeric("core") # use core in the SBATCH command!!!

load(paste(name,"/Data.Rdata",sep=""))
load(paste( name, "/NonParamResult.Rdata", sep = "" ))

a=1/(nrow(Data))^2
s = 5000
burn = 1000
alpha = 1
num_iter = 35000
max_clusters = ncol(Data)

set.seed(copy*20000+Ci*10+5+1/a)

startC=initialC[[Ci]]

t0=proc.time()
results <- gibbs_sampler(y = Data, alpha = alpha, max_clusters = max_clusters, s = s, burn = burn, name = name, initialC = startC, n_iterations = num_iter)
t1=proc.time()-t0
cput=t1[1]
save(results,cput, file = paste( name, "/IMMclusters", Ci, ".Rdata", sep=""))

