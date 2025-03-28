source("LibraryCmdline.R")
source("ModifiedAllfcns.R")
current<-cmdline.numeric("current")
name<-cmdline.strings("name")
Jname<-cmdline.strings("Jname")
parent<-floor(current/2)
copy<-cmdline.numeric("copy")
load(paste(name,"/Data.Rdata",sep=""))

a=1/(nrow(Data)^2)
print(paste("First lambda =",L))
set.seed(copy*1200+current*10+3)
print(paste('parent is',parent))
print(paste('current is',current))
s = 10000 

if(0==parent)
{	
  t0=proc.time()
  print(paste('t0=',t0))
  G=seq(1,ncol(Data),1)
  out=split(Data,G,s,a,L)	#change here with number of splits
  time=proc.time()-t0
  old.result=out
  # sub.L=paste("bsub -q week -J tree.", Jname," -o ", "~/LSFouts/", Jname, ".out
  #       R CMD BATCH --vanilla --args --copy=", copy," --current=", current*2," --name=",
  #       name, " --Jname=", Jname," Tree.R ", name,"/Outputs/Tree_",current*2,".out",sep="")
  # sub.R=paste("bsub -q week -J tree.", Jname," -o ", "~/LSFouts/", Jname, ".out
  #       R CMD BATCH --vanilla --args --copy=", copy," --current=", current*2+1," --name=",
  #       name, " --Jname=", Jname, " Tree.R ", name,"/Outputs/Tree_",current*2+1,".out",sep="")
  #sub.L=paste("sbatch -p general -J ", Jname," -o ", name, "SLURMouts/", Jname, 
  #".out.%j  --wrap=\"R CMD BATCH --vanilla --args --copy=", copy," --current=", current*2," --name=",
  #            name, " --Jname=", Jname," Tree.R ", name,"/Outputs/Tree_",current*2,".out\"",sep="")
  #sub.R=paste("sbatch -p general -J tree.", Jname," -o SLURMouts/", Jname, 
  #".out.%j  --wrap=\"R CMD BATCH --vanilla --args --copy=", copy," --current=", current*2+1," --name=",
  #            name, " --Jname=", Jname, " Tree.R ", name,"/Outputs/Tree_",current*2+1,".out\"",sep="")
  
  #change here
  sub.L = paste("sbatch -t 00:05:00 -N 1 -n 1 -J",Jname, "_%j -o ", name, "/SLURMouts/", Jname,
                "_%j.out.", current*2,"  --wrap=\"R  CMD BATCH --vanilla --args --copy=", copy," --current=", current*2,
                " --name=",name, " --Jname=", Jname, " Tree.R ", name, "/Outputs/tree_",current*2,".out\"",sep="" )
  sub.R = paste("sbatch -t 00:05:00 -N 1 -n 1 -J",Jname, "_%j -o ", name, "/SLURMouts/", Jname,
                "_%j.out.", current*2+1,"  --wrap=\"R  CMD BATCH --vanilla --args --copy=", copy," --current=", current*2+1,
                " --name=",name, " --Jname=", Jname, " Tree.R ", name, "/Outputs/tree_",current*2+1,".out\"",sep="" )	
  
  save(old.result, time, file=paste(name,"/Tree_",current,".Rdata",sep=""))
  #   save(old.result, file=paste(name,"/Tree_",current,".Rdata",sep=""))
  print ("sub.L=")
  print (sub.L)
  print ("sub.R=")
  print (sub.R)
  system(sub.L)	
  system(sub.R)
  #Sys.sleep(30)
  q(save="no")		
}	

load(paste(name,"/Tree_",parent,".Rdata",sep=""))
t0=proc.time()
print(paste('t0=',t0))
if (1==current%%2) 
{
  out=split(Data,old.result$right,s,a,L)	
}
if (0==current%%2) 
{
  out=split(Data,old.result$left,s,a,L)	
}
time=proc.time()+time #This seems questionable
old.result=out

if (0!=length(out$left) & 0!=length(out$right))				
{	
  # sub.L=paste("bsub -q week -J tree.", Jname," -o ", "~/LSFouts/", Jname, ".out
  #       R CMD BATCH --vanilla --args --copy=", copy," --name=", name, " --Jname=", Jname,
  #       " --current=", current*2," Tree.R ", name,"/Outputs/Tree_",current*2,".out",sep="")
  # sub.R=paste("bsub -q week -J tree.", Jname," -o ", "~/LSFouts/", Jname, ".out
  #       R CMD BATCH --vanilla --args --copy=", copy," --name=",name, " --Jname=", Jname,
  #       " --current=", current*2+1," Tree.R ", name,"/Outputs/Tree_",current*2+1,".out",sep="")
  #sub.L=paste("sbatch -p general -J tree.", Jname," -o SLURMouts/", Jname, 
  #            ".out.%j  --wrap=\"R CMD BATCH --vanilla --args --copy=", copy," --name=", name, " --Jname=", Jname,
  #            " --current=", current*2," Tree.R ", name,"/Outputs/Tree_",current*2,".out\"",sep="")
  #sub.R=paste("sbatch -p general -J tree.", Jname," -o SLURMouts/", Jname, 
  #".out.%j  --wrap=\"R CMD BATCH --vanilla --args --copy=", copy,
  #" --name=",name, " --Jname=", Jname, " --current=", current*2+1,
  #" Tree.R ", name,"/Outputs/Tree_",current*2+1,".out\"",sep="")
  
  sub.L = paste("sbatch -t 00:05:00  -N 1 -n 1 -J ",Jname, "_%j -o ", name, "/SLURMouts/", Jname,
                "_%j.out.", current*2,"  --wrap=\"R  CMD BATCH --vanilla --args --copy=", copy," --current=", current*2,
                " --name=",name, " --Jname=", Jname, " Tree.R ", name, "/Outputs/tree_",current*2,".out\"",sep="" )
  sub.R = paste("sbatch -t 00:05:00 -N 1 -n 1 -J ",Jname, "_%j -o ", name, "/SLURMouts/", Jname,
                "_%j.out.", current*2+1,"  --wrap=\"R  CMD BATCH --vanilla --args --copy=", copy," --current=", current*2+1,
                " --name=",name, " --Jname=", Jname, " Tree.R ", name, "/Outputs/tree_",current*2+1,".out\"",sep="" )	
  
  save(old.result, time, file=paste(name,"/Tree_",current,".Rdata",sep=""))
  #   save(old.result, file=paste(name,"/Tree_",current,".Rdata",sep=""))
  system(sub.L)  
  system(sub.R)	
  #Sys.sleep(30)
}

if (0==length(out$left) | 0==length(out$right))  
{
  
  save(old.result, time, file=paste(name,"/Tree_",current,".Rdata",sep=""))
  CPUtime=time[1]
  write.csv(CPUtime, file=paste(name,"/BranchCPUtime_",current,".csv",sep=""), row.names=F)
}






