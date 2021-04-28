# temporal_Ne_linked

This repository contains the datasets and computing codes of the accepted manuscript "Contemporary Ne estimation using temporally-spaced data with linked loci. " Molecular Ecology Resources (2021)

maintained by Tin-Yu Hui <tin-yu.hui11@imperial.ac.uk>

The two datasets ```col.RData``` and ```gam.RData``` for *An. coluzzii* and *An. gambiae* locate under the ```/data/``` folder. ```Script.R``` is the main R code to generate the results in Table 2 and Supplementary Information (Part IV). There is also a C source file ```/cpp/r_matrix.c``` which needs to be compiled before running. The compiled shared object will then be called from R. 

To compile the C source file, open ```cmd``` (Windows) or the equivalent terminal in Linux/MacOS, go to the target directory, and run the following command: 
```
R CMD SHLIB r_matrix.c -fopenmp
```
The flag ```-fopenmp``` is required to enable parallelisation via OpenMP. This creates a shared object with extension ```.dll``` (Windows) or ```.so``` (Linux/MacOS). The shared object will be linked to the R session via ```dyn.load()```, and the C functions within it will be called via ```.Call()```. For details please consult the document "Writing R extensions" on CRAN. 
