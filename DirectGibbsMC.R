source("ModifiedAllfcns.R")
source("LibraryCmdline.R")

Ci<-cmdline.numeric("Ci")
name<-cmdline.strings("name")
copy<-cmdline.numeric("copy")
CORE<-cmdline.numeric("core") # use core in the SBATCH command!!!

load(paste(name,"/Data.Rdata",sep=""))
load(paste( name, "/FixGibbsResult.Rdata", sep = "" ))

s = 5000
burn = 1000
num_iter = 35000

TempFile = paste0(name, "/GibbsTemp", Ci, ".Rdata")

a=1/(nrow(Data))^2

set.seed(copy*20000+Ci*10+5+1/a)

startC=newC0s[Ci,]

t0=proc.time()

print(CORE)
GibbsMC = Direct_Gibbs_Sampler(Data, startC, Scan = num_iter, a = 1/25, CORE = CORE, TempFile, L, burn, s)
t1=proc.time()-t0
cput=t1[1]
save(GibbsMC, cput, file=paste(name,"/GibbsOnly",Ci,".Rdata",sep=""))
