# temporal_Ne_linked

This repository contains the datasets and computing codes of the now published paper "Contemporary Ne estimation using temporally-spaced data with linked loci. "

Please use the following citation: 
Hui, T. Y. J., Brenas, J. H., & Burt, A. (2021). Contemporary N e estimation using temporally spaced data with linked loci. Molecular Ecology Resources, 21(7), 2221-2230. 

https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13412

Uploaded by Tin-Yu Hui <tin-yu.hui11@imperial.ac.uk>, April 2021. 

The two datasets ```col.RData``` and ```gam.RData``` for *An. coluzzii* and *An. gambiae* locate under the ```/data/``` folder. ```Script.R``` is the main R code to generate the results in Table 2 and Supplementary Information (Part IV). There is also a C source file ```/cpp/r_matrix.c``` which needs to be compiled before running. The compiled shared object will then be called from R. 

To compile the C source file, open ```cmd``` (Windows) or the equivalent terminal in Linux/MacOS, go to the target directory, and run the following command: 
```
R CMD SHLIB r_matrix.c -fopenmp
```
The flag ```-fopenmp``` is required to enable parallelisation via OpenMP. This creates a shared object with extension ```.dll``` (Windows) or ```.so``` (Linux/MacOS). The shared object will be linked to the R session via ```dyn.load()```, and the C functions within it will be called via ```.Call()```. For details please consult the document "Writing R extensions" on CRAN. 
