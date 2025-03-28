source("LibraryCmdline.R")
source("Source.R")
source("ModifiedAllfcns.R")

# Global Params
shortname <- cmdline.strings("shortname") # For simulated data, shortname = Test
copy <- cmdline.numeric("copy") # For possible multiple tests
CORE <- cmdline.numeric("core") # Number of cores needed for mclapply
GibbsRun <- cmdline.strings("GibbsRun") # true or false for whether to do third step of our algortihm
algorithm <- cmdline.strings("algorithm") # choose ours, DP, or direct
PostD <- cmdline.numeric("PostD") # sample(s) with treatment
L <- cmdline.numeric("L")
folder <- cmdline.strings("folder")
n <- cmdline.numeric("n")
user <- cmdline.strings("user")

set.seed(copy*10)

PostD = PostD; K.control = 10; K.mut = 5; J = 5; n = n; m = 5000
a = 1/(J^2) # Penalty term from Dirichlet

# Create a folder "Tests" and a subfolder "Tests/Data" for the simulated data sets
Jname = paste(shortname, copy, "_a",1/a,sep="") #e.g. Jname = Test1_25
name = paste( folder, "/", Jname, sep = "" )
dir.create(folder)
dir.create(paste(folder,"/Data", sep=""))

# Create simulated data
createtestdata(copy, K0 = K.control, newK = K.mut, d = J, n = n, sitect = m, tests = max(PostD),
               postD = PostD, folder)

load(paste(folder,"/Data/TestData_",copy,".Rdata",sep=""))

Pairs = c( "t1t2", "t1t3", "t2t3", "t1t3_D", "t2t3_D" )
n = ncol(rawdata)/max(PostD)

# Create directory 
dir.create(name)

# Run Main Script with pre-processing, processing, and post processing
source("Main.R")