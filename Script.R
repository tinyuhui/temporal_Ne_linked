# Contemporary Ne estimation using temporally-spaced data with linked loci
# Ne AND C.I. ESTIMATION FOR AN. COLUZZII AND AN. GAMBIAE
# 27/04/2021 TIN-YU J HUI

# REQUIRE PACKAGE(S)
require(compiler)
fsum<-function(x) {.Primitive('sum')(x)}

# IF USING MICROSOFT R OPEN
setMKLthreads(20)

#################################################
# CHROMOSOME 3R OF COLUZZII
#################################################
# LOAD COLUZIII DATA
load('data/col.RData')
ls()
# FIND SAMPLE SIZES
s0<-nrow(subset_3R_col_2012)
st<-nrow(subset_3R_col_2014)
c(s0, st)

# NUMBER OF LOCI
K<-ncol(subset_3R_col_2012)
K
x0<-apply(subset_3R_col_2012, 2, fsum)/(2*s0)
xt<-apply(subset_3R_col_2014, 2, fsum)/(2*st)
# Fa AND CORRESPONDING Ne
Fa_3R_col<-fsum((x0-xt)^2/(x0*(1-x0)))/K
Na_3R_col<-20/(2*(Fa_3R_col-0.5/s0-0.5/st))
c(Fa_3R_col, Na_3R_col)
# Fb AND CORRESPONDING Ne
Fb_3R_col<-fsum((x0-xt)^2)/fsum(x0*(1-x0))
Nb_3R_col<-20/(2*(Fb_3R_col-0.5/s0-0.5/st))
c(Fb_3R_col, Nb_3R_col)

# CALCULATE R MATRIX AND EIGENVALUES FOR Fa
dyn.load('cpp/r_matrix.dll')
R<-.Call('cal_corr_matrix_unphased', subset_3R_col_2012, subset_3R_col_POS, 0.014, Na_3R_col, s0, st, 20, 20)
dyn.unload('cpp/r_matrix.dll')
e<-eigen(R, only.values=TRUE)
gc()

# C.I. WITH Fa
gen_Q2<-function(eigenvalues)
{
working_eigenvalues<-eigenvalues[eigenvalues>0]
working_eigenvalues<-working_eigenvalues*sum(eigenvalues)/sum(working_eigenvalues)
dyn.load('cpp/r_matrix.dll')
y<-.Call('generate_Q2', working_eigenvalues, 50000)
dyn.unload('cpp/r_matrix.dll')
return(y)
}
# EMPIRICAL Q2 DISTRIBUTION. YOU WILL GET A SLIGHTLY DIFFERENT Q2 BECAUSE OF RANDOMNESS
Q2<-gen_Q2(e$values)
# C.I. FOR F
Fa_3R_col_lower<-K*Fa_3R_col/quantile(Q2, 0.975)
Fa_3R_col_upper<-K*Fa_3R_col/quantile(Q2, 0.025)
c(Fa_3R_col_lower, Fa_3R_col, Fa_3R_col_upper)
# C.I. FOR Ne
Na_3R_col_lower<-20/(2*(Fa_3R_col_upper-0.5/s0-0.5/st))
Na_3R_col_upper<-20/(2*(Fa_3R_col_lower-0.5/s0-0.5/st))
c(Na_3R_col_lower, Na_3R_col, Na_3R_col_upper)

# C.I. WITH Fb
rm(R); rm(Q2); gc(); 
w<-x0*(1-x0)/fsum(x0*(1-x0))
W_half<-diag(sqrt(w))
dyn.load('cpp/r_matrix.dll')
R<-.Call('cal_corr_matrix_unphased', subset_3R_col_2012, subset_3R_col_POS, 0.014, Nb_3R_col, s0, st, 20, 20)
dyn.unload('cpp/r_matrix.dll')
e<-eigen(K*W_half%*%R%*%t(W_half), only.values=TRUE)
gc()
Q2<-gen_Q2(e$values)
# C.I. FOR F
Fb_3R_col_lower<-K*Fb_3R_col/quantile(Q2, 0.975)
Fb_3R_col_upper<-K*Fb_3R_col/quantile(Q2, 0.025)
c(Fb_3R_col_lower, Fb_3R_col, Fb_3R_col_upper)
# C.I. FOR Ne
Nb_3R_col_lower<-20/(2*(Fb_3R_col_upper-0.5/s0-0.5/st))
Nb_3R_col_upper<-20/(2*(Fb_3R_col_lower-0.5/s0-0.5/st))
c(Nb_3R_col_lower, Nb_3R_col, Nb_3R_col_upper)

#################################################
# CHROMOSOME 3L OF COLUZZII
#################################################
# CLEAN UP VARIABLES
rm(K); rm(R); rm(Q2); rm(x0); rm(xt); rm(w); rm(W_half); gc(); 
# THEN REPEAT THE SAME ANALYSIS WITH 3L DATA
K<-ncol(subset_3L_col_2012)
K
x0<-apply(subset_3L_col_2012, 2, fsum)/(2*s0)
xt<-apply(subset_3L_col_2014, 2, fsum)/(2*st)
# Fa AND CORRESPONDING Ne
Fa_3L_col<-fsum((x0-xt)^2/(x0*(1-x0)))/K
Na_3L_col<-20/(2*(Fa_3L_col-0.5/s0-0.5/st))
c(Fa_3L_col, Na_3L_col)
# Fb AND CORRESPONDING Ne
Fb_3L_col<-fsum((x0-xt)^2)/fsum(x0*(1-x0))
Nb_3L_col<-20/(2*(Fb_3L_col-0.5/s0-0.5/st))
c(Fb_3L_col, Nb_3L_col)

# CALCULATE R MATRIX AND EIGENVALUES FOR Fa
dyn.load('cpp/r_matrix.dll')
R<-.Call('cal_corr_matrix_unphased', subset_3L_col_2012, subset_3L_col_POS, 0.014, Na_3L_col, s0, st, 20, 20)
dyn.unload('cpp/r_matrix.dll')
e<-eigen(R, only.values=TRUE)
gc()

# C.I. WITH Fa
# EMPIRICAL Q2 DISTRIBUTION. YOU WILL GET A SLIGHTLY DIFFERENT Q2 BECAUSE OF RANDOMNESS
Q2<-gen_Q2(e$values)
# C.I. FOR F
Fa_3L_col_lower<-K*Fa_3L_col/quantile(Q2, 0.975)
Fa_3L_col_upper<-K*Fa_3L_col/quantile(Q2, 0.025)
c(Fa_3L_col_lower, Fa_3L_col, Fa_3L_col_upper)
# C.I. FOR Ne
Na_3L_col_lower<-20/(2*(Fa_3L_col_upper-0.5/s0-0.5/st))
Na_3L_col_upper<-20/(2*(Fa_3L_col_lower-0.5/s0-0.5/st))
c(Na_3L_col_lower, Na_3L_col, Na_3L_col_upper)

# CALCULATE R MATRIX AND EIGENVALUES FOR Fb
rm(R); rm(Q2); gc(); 
w<-x0*(1-x0)/fsum(x0*(1-x0))
W_half<-diag(sqrt(w))
dyn.load('cpp/r_matrix.dll')
R<-.Call('cal_corr_matrix_unphased', subset_3L_col_2012, subset_3L_col_POS, 0.014, Nb_3L_col, s0, st, 20, 20)
dyn.unload('cpp/r_matrix.dll')
e<-eigen(K*W_half%*%R%*%t(W_half), only.values=TRUE)
gc()

# C.I. WITH Fb
Q2<-gen_Q2(e$values)
# C.I. FOR F
Fb_3L_col_lower<-K*Fb_3L_col/quantile(Q2, 0.975)
Fb_3L_col_upper<-K*Fb_3L_col/quantile(Q2, 0.025)
c(Fb_3L_col_lower, Fb_3L_col, Fb_3L_col_upper)
# C.I. FOR Ne
Nb_3L_col_lower<-20/(2*(Fb_3L_col_upper-0.5/s0-0.5/st))
Nb_3L_col_upper<-20/(2*(Fb_3L_col_lower-0.5/s0-0.5/st))
c(Nb_3L_col_lower, Nb_3L_col, Nb_3L_col_upper)

#################################################
# 3R AND 3L COMBINED FOR COLUZZII
#################################################
rm(K); rm(R); rm(Q2); rm(x0); rm(xt); rm(w); rm(W_half); gc(); 
# COMBINE 3R AND 3L
subset_col_2012<-cbind(subset_3R_col_2012, subset_3L_col_2012)
subset_col_2014<-cbind(subset_3R_col_2014, subset_3L_col_2014)
# COMBINE POS
subset_col_POS<-c(subset_3R_col_POS, 53200195+subset_3L_col_POS)
mode(subset_col_2012)<-'integer'
mode(subset_col_2014)<-'integer'
mode(subset_col_POS)<-'integer'
# TOTAL LOCI NUMBER
K<-ncol(subset_col_2012)
K
x0<-apply(subset_col_2012, 2, fsum)/(2*s0)
xt<-apply(subset_col_2014, 2, fsum)/(2*st)
# THE TWO F STATISTICS
Fa_col<-fsum((x0-xt)^2/(x0*(1-x0)))/K
Na_col<-20/(2*(Fa_col-0.5/s0-0.5/st))
c(Fa_col, Na_col)
Fb_col<-fsum((x0-xt)^2)/fsum(x0*(1-x0))
Nb_col<-20/(2*(Fb_col-0.5/s0-0.5/st))
c(Fb_col, Nb_col)

# CALCULATE R MATRIX AND EIGENVALUES FOR Fa
dyn.load('cpp/r_matrix.dll')
R<-.Call('cal_corr_matrix_unphased', subset_col_2012, subset_col_POS, 0.014, Na_col, s0, st, 20, 20)
dyn.unload('cpp/r_matrix.dll')
e<-eigen(R, only.values=TRUE)
gc()

# C.I. WITH Fa
Q2<-gen_Q2(e$values)
# C.I. FOR F
Fa_col_lower<-K*Fa_col/quantile(Q2, 0.975)
Fa_col_upper<-K*Fa_col/quantile(Q2, 0.025)
c(Fa_col_lower, Fa_col, Fa_col_upper)
# C.I. FOR Ne
Na_col_lower<-20/(2*(Fa_col_upper-0.5/s0-0.5/st))
Na_col_upper<-20/(2*(Fa_col_lower-0.5/s0-0.5/st))
c(Na_col_lower, Na_col, Na_col_upper)

# CALCULATE R MATRIX AND EIGENVALUES FOR Fb
rm(R); rm(Q2); gc(); 
w<-x0*(1-x0)/fsum(x0*(1-x0))
W_half<-diag(sqrt(w))
dyn.load('cpp/r_matrix.dll')
R<-.Call('cal_corr_matrix_unphased', subset_col_2012, subset_col_POS, 0.014, Nb_col, s0, st, 20, 20)
dyn.unload('cpp/r_matrix.dll')
e<-eigen(K*W_half%*%R%*%t(W_half), only.values=TRUE)
gc()

# C.I. WITH Fb
Q2<-gen_Q2(e$values)
# C.I. FOR F
Fb_col_lower<-K*Fb_col/quantile(Q2, 0.975)
Fb_col_upper<-K*Fb_col/quantile(Q2, 0.025)
c(Fb_col_lower, Fb_col, Fb_col_upper)
# C.I. FOR Ne
Nb_col_lower<-20/(2*(Fb_col_upper-0.5/s0-0.5/st))
Nb_col_upper<-20/(2*(Fb_col_lower-0.5/s0-0.5/st))
c(Nb_3L_col_lower, Nb_col, Nb_col_upper)
