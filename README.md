# temporal_Ne_linked

This repository contains the datasets and computing codes of the accepted manuscript "Contemporary Ne estimation using temporally-spaced data with linked loci. " Molecular Ecology Resources (2021)

maintained by Tin-Yu Hui <tin-yu.hui11@imperial.ac.uk>

The two datasets ```col.RData``` and ```gam.RData``` for *An. coluzzii* and *An. gambiae* locate under the ```/data/``` folder. ```Script.R``` contains the main R code to generate the results in Table 2 and Supplementary Information (Part IV). There is also a C source file ```/cpp/r_matrix.c``` which needs to be compiled before running. The compiled shared object will be callable from R. 

To compiler the C source file, open ```cmd``` (Windows) or its equivalent in other OS, go to the targeted directory, and run the following command: 
```
R CMD SHLIB r_matrix.c -fopenmp
```
The flag ```-fopenmp``` is required for OpenMP parallelisation. This creates a shared object with extension ```.dll``` (Windows) or ```.so``` (Linux). 
