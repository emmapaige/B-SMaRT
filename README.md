# B-SMaRT

## Getting Started

### Prerequisites

Ensure R (version 4.1.3 or later) is installed, along with the following R packages:

```r
install.packages(c("R.oo", "bmixture", "coda", "parallel", "stringr", "fields", "matrixStats", "gtools", "patchwork", "dplyr"))
```

### Installation

Clone this repository: 
git clone [https://github.com/emmapaige/B-SMaRT.git]

### The Workflow:

1. **Setup (`TestRun.R`)**:  
   - Generates simulated genomic data.  
   - Sets up global parameters.

2. **Main Script (`Main.R`)**:  
   Performs the main analysis, which includes:
   - **Preprocessing**: Combines homogenous sites and outputs the joint data matrix with all time samples.  
   - **Processing**: Generates assignment labels using any of the three methods outlined in the manuscript.
   - **Postprocessing**: Computes Hellinger distances given the assignment labels. Outputs substitution sites, noise sites, potential sites, and the computed threshold in the output file `FinalFixOne1.Rdata`.

### Supporting Scripts:

1. `Tree.R`: The workflow for creating the 2cSCMH.
2. `source.R`: Enables automatic hierarchical divisive tree.
3. `ModifiedAllfcns.R`: Provides all the necessary functions for the Main Script.
4. `LibraryCmdline.R`: Provides all the source code needed for SLURM command line input.
5. `DirectGibbsMC.R` : Workflow for obtaining the assignment labels using direct Gibbs sampler only.
6. `NonParamMC.R` : Workflow for obtaining the assignment labels using a Dirichlet Process.
7. `FixGibbsMC.R` : Workflow for the third step of our algorithm, the single Gibbs sampler scan.
---
## How to Run
`Main.R` is designed to be run on a SLURM cluster. To execute, submit something similar to the below SLURM script:

```#!/bin/bash

#SBATCH -t 1- 
#SBATCH -N 1 
#SBATCH -n 1

mkdir -p ./Rlogs

module add r/4.1.3

R CMD BATCH --vanilla --args --shortname=Test --core=1 --copy=1 --GibbsRun=T --algorithm=ours --PostD=4  --L=1 --folder=Tests --n=300 --user=<username> TestRun.R ./Rlogs/TestRun.out
```

### Inputs for Submission
The execution script above requires six input parameters, listed after --args. The inputs are as follows:

- **`shortname`**: A short identifier for the current run (e.g., `"Test"`).
- **`core`**: Number of cores needed.
- **`copy`**: Copy number in case multiple runs are performed.
- **`GibbsRun`**: Boolean indicating whether to use the third step of our algorithm, the single-scan direct Gibbs (T or F).
- **`algorithm`**: The algorithm to use for generating the assignment labels:
  - `"ours"`: Single-coordinate MH, block MH, and one-scan Gibbs sampler.
  - `"direct"`: Direct Gibbs sampler.
  - `"DP"`: Dirichlet Process.
- **`PostD`**: Number indicating which time sample is the treated one.
- **`L`**: $\lambda$ prior parameter for poisson prior
- **`folder`**: Name of the folder the results save to
- **`n`**: Number of locations on the genome
- **`user`**: username for job submissions

