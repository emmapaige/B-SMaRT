source("ModifiedAllfcns.R")
source("LibraryCmdline.R")

Ci<-cmdline.numeric("Ci")
name<-cmdline.strings("name")
copy<-cmdline.numeric("copy")
CORE<-cmdline.numeric("core") # use core in the SBATCH command!!!

print(Ci)
print(name)
print(copy)
print(CORE)

load(paste(name,"/Data.Rdata",sep=""))
a=1/(nrow(Data))^2

set.seed(copy*20000+Ci*10+5+1/a)
load(paste(name,"/FixGibbsResult.Rdata",sep=""))
tempfile = paste(name, "/FixGibbsMC", Ci, "temp.Rdata", sep="")

t0=proc.time()
initialC=newC0s[Ci,]
GibbsMC=GibbsFixScan_Emma(Data=Data, initialC=initialC, Scan=1, a=a, Record=T, CORE=CORE, TempFile=tempfile,L=L) # change here 07/18 back to one iteration
t1=proc.time()-t0
print(t1)
cput=t1[1]
save(initialC,GibbsMC,cput,file=paste(name,"/FixGibbsMC",Ci,".Rdata",sep=""))

Sys.sleep(60)


