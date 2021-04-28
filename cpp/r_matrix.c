#include <R.h>
#include <math.h>
#include <omp.h>
#include <Rinternals.h>
#include <Rdefines.h>

// FUNCTION TO CALCULATE PAIRWISE RECOMBINATION RATE FROM MORGAN PER MB
double recom_rate(int x, int y, double morgan_per_mb)
{
	int phy_dist=abs(x-y);
	double d=(double) (phy_dist/1e6)*morgan_per_mb;
	return 0.5*(1-exp(-2*d));
}

// ADJUSTED r
double cal_adj_r(double r, double s)
{
	double temp=(r*r-1/s)/(1-1/s);
	if (temp<=0)
	{return 0;}
	if (r>0 && temp>0)
	{return sqrt(temp);}
	if (r<0 && temp>0)
	{return (-1)*sqrt(temp);}
}

// FUNCTIONS TO CALCULATE PAIRWISE r, AND RETURN A MATRIX OF ALL THE r

// r phased, MLE
// x[] AND y[] CONTAIN {0, 1} ALLELE ONLY, LENGTH 2*s
double r_phased(int x[], int y[], int s)
{
	double p1=0; double p2=0; double p3=0;
	int temp=0;
	for (int i=0; i<2*s; i++)
	{
		temp=2*x[i]+y[i];
		if (temp==3) {p1=p1+1;}
		if (temp==2) {p2=p2+1;}
		if (temp==1) {p3=p3+1;}
	}
	double p4=2*s-p1-p2-p3;
	return (p1*p4-p2*p3)/sqrt((p1+p2)*(p3+p4)*(p1+p3)*(p2+p4));
}

// r unphased, EM?
// x[] AND y[] CONTAIN {0, 1, 2}, LENGTH s
double r_unphased_EM(int x[], int y[], int s)
{
	int temp=0;
	double g1=0; double g2=0; double g3=0; 
	double g4=0; double g5=0; double g6=0; 
	double g7=0; double g8=0; double g9=0;
	for (int i=0; i<s; i++)
	{
		// CALCULATE THE GENOTYPE TABLE
		temp=3*x[i]+y[i];
		if (temp==8) {g1=g1+1;}
		if (temp==7) {g2=g2+1;}
		if (temp==6) {g3=g3+1;}
		if (temp==5) {g4=g4+1;}
		if (temp==4) {g5=g5+1;}
		if (temp==3) {g6=g6+1;}
		if (temp==2) {g7=g7+1;}
		if (temp==1) {g8=g8+1;}
		if (temp==0) {g9=g9+1;}
	}
	double p1=0.25; double p2=0.25; double p3=0.25; double p4=0.25;
	double u=0; double ds=(double) s;
	while ((u-g5*p1*p4/(p1*p4+p2*p3))*(u-g5*p1*p4/(p1*p4+p2*p3))>0.001)
	{
		u=g5*p1*p4/(p1*p4+p2*p3);
		p1=(2*g1+g2+g4+u)/(2*ds);
		p2=(2*g3+g2+g6+g5-u)/(2*ds);
		p3=(2*g7+g8+g4+g5-u)/(2*ds);
		p4=(2*g9+g8+g6+u)/(2*ds);
	}
	
	return (p1*p4-p2*p3)/sqrt((p1+p2)*(p3+p4)*(p1+p3)*(p2+p4));
}

// r unphased, Burrows'
// x[] AND y[] CONTAIN {0, 1, 2}, LENGTH s
double r_unphased_burrows(int x[], int y[], int s)
{
	double delta_AB=0; 
	double p_A=0; double p_B=0; 
	double h_A=0; double h_B=0; 
	for (int i=0; i<s; i++)
	{
			delta_AB=delta_AB+(double) x[i]*y[i];
			p_A=p_A+(double) x[i];
			p_B=p_B+(double) y[i];
			h_A=h_A+(double) (x[i]==2);
			h_B=h_B+(double) (y[i]==2);
	}
	double ds=(double) s;
	p_A=p_A/(2*ds); p_B=p_B/(2*ds);
	h_A=h_A/ds; h_B=h_B/ds;
	delta_AB=(delta_AB/(2*ds)-2*p_A*p_B)*(ds/(ds-1));
	return (delta_AB*delta_AB)/((p_A-2*p_A*p_A+h_A)*(p_B-2*p_B*p_B+h_B));
}

// CALCULATE r MATRIX, PHASED
SEXP cal_r_matrix_phased(SEXP haplotype, SEXP ncpu)
{
	int num_thread=asInteger(ncpu);
	// OMP STUFF
	omp_set_dynamic(0);
	omp_set_num_threads(num_thread);
	Rprintf("Using %d cores. \n", num_thread);
	// DEFINE OUTPUT MATRIX
	int loci=INTEGER(GET_DIM(haplotype))[1];
	int s=INTEGER(GET_DIM(haplotype))[0]/2;
	SEXP result=PROTECT(allocMatrix(REALSXP, loci, loci));
	// THE POINTERS FOR EASY REFERENCE
	double *result_ptr; int *haplotype_ptr;
	result_ptr=REAL(result);
	haplotype_ptr=INTEGER(haplotype);
	int **x; int **y;
	x=(int **) R_alloc(num_thread, sizeof(int *));
	y=(int **) R_alloc(num_thread, sizeof(int *));
	for (int i=0; i<num_thread; i++)
	{
		x[i]=(int *) R_alloc(2*s, sizeof(int));
		y[i]=(int *) R_alloc(2*s, sizeof(int));
	}
	#pragma omp parallel
	{
		int tid=omp_get_thread_num();
		double temp_r=0;
		#pragma omp for
		for (int i=0; i<loci; i++)
		{
			for (int k=0; k<2*s; k++) {x[tid][k]=haplotype_ptr[i*2*s+k];}
			for (int j=i; j<loci; j++)
			{
				if (i==j)
				{
					result_ptr[i*loci+j]=1;
				}
				else
				{	
				for (int k=0; k<2*s; k++) {y[tid][k]=haplotype_ptr[j*2*s+k];}
				temp_r=r_phased(x[tid], y[tid], s);
				result_ptr[i*loci+j]=temp_r;
				result_ptr[j*loci+i]=temp_r;
				}
			}
		}
	}
	UNPROTECT(1);
	return result;
}

// CALCILATE r MATRIX, UNPHASED,, V2
// R_alloc SHOULD BE OUTSIDE THE omp LOOP
SEXP cal_r_matrix_unphased_EM(SEXP genotype, SEXP ncpu)
{
	int num_thread=asInteger(ncpu);
	// OMP STUFF
	omp_set_dynamic(0);
	omp_set_num_threads(num_thread);
	Rprintf("Using %d cores. \n", num_thread);
	// DEFINE OUTPUT MATRIX
	int loci=INTEGER(GET_DIM(genotype))[1];
	int s=INTEGER(GET_DIM(genotype))[0];
	SEXP result=PROTECT(allocMatrix(REALSXP, loci, loci));
	// THE POINTERS FOR EASY REFERENCE
	double *result_ptr; int *genotype_ptr;
	result_ptr=REAL(result);
	genotype_ptr=INTEGER(genotype);
	int **x; int **y;
	x=(int **) R_alloc(num_thread, sizeof(int *));
	y=(int **) R_alloc(num_thread, sizeof(int *));
	for (int i=0; i<num_thread; i++)
	{
		x[i]=(int *) R_alloc(s, sizeof(int));
		y[i]=(int *) R_alloc(s, sizeof(int));
	}
	#pragma omp parallel
	{
		int tid=omp_get_thread_num();
		double temp_r=0;
		#pragma omp for
		for (int i=0; i<loci; i++)
		{
			#pragma omp simd
			for (int k=0; k<s; k++) {x[tid][k]=genotype_ptr[i*s+k];}
			for (int j=i; j<loci; j++)
			{
				if (i==j)
				{
					result_ptr[i*loci+j]=1;
				}
				else
				{
					#pragma omp simd
					for (int k=0; k<s; k++) {y[tid][k]=genotype_ptr[j*s+k];}
					temp_r=r_unphased_EM(x[tid], y[tid], s);
					result_ptr[i*loci+j]=temp_r;
					result_ptr[j*loci+i]=temp_r;
				}
			}
		}
	}
	UNPROTECT(1);
	return result;
}

// CALCILATE r MATRIX, UNPHASED,, V2
// R_alloc SHOULD BE OUTSIDE THE omp LOOP
SEXP cal_r_matrix_unphased_burrows(SEXP genotype, SEXP ncpu)
{
	int num_thread=asInteger(ncpu);
	// OMP STUFF
	omp_set_dynamic(0);
	omp_set_num_threads(num_thread);
	Rprintf("Using %d cores. \n", num_thread);
	// DEFINE OUTPUT MATRIX
	int loci=INTEGER(GET_DIM(genotype))[1];
	int s=INTEGER(GET_DIM(genotype))[0];
	SEXP result=PROTECT(allocMatrix(REALSXP, loci, loci));
	// THE POINTERS FOR EASY REFERENCE
	double *result_ptr; int *genotype_ptr;
	result_ptr=REAL(result);
	genotype_ptr=INTEGER(genotype);
	int **x; int **y;
	x=(int **) R_alloc(num_thread, sizeof(int *));
	y=(int **) R_alloc(num_thread, sizeof(int *));
	for (int i=0; i<num_thread; i++)
	{
		x[i]=(int *) R_alloc(s, sizeof(int));
		y[i]=(int *) R_alloc(s, sizeof(int));
	}
	#pragma omp parallel
	{
		int tid=omp_get_thread_num();
		double temp_r=0;
		#pragma omp for
		for (int i=0; i<loci; i++)
		{
			#pragma omp simd
			for (int k=0; k<s; k++) {x[tid][k]=genotype_ptr[i*s+k];}
			for (int j=i; j<loci; j++)
			{
				if (i==j)
				{
					result_ptr[i*loci+j]=1;
				}
				else
				{
					#pragma omp simd
					for (int k=0; k<s; k++) {y[tid][k]=genotype_ptr[j*s+k];}
					temp_r=r_unphased_burrows(x[tid], y[tid], s);
					result_ptr[i*loci+j]=temp_r;
					result_ptr[j*loci+i]=temp_r;
				}
			}
		}
	}
	UNPROTECT(1);
	return result;
}

SEXP cal_recom_rate_matrix(SEXP POS, SEXP morgan_per_mb, SEXP ncpu)
{
	int num_thread=asInteger(ncpu);
	// OMP STUFF
	omp_set_dynamic(0);
	omp_set_num_threads(num_thread);
	Rprintf("Using %d cores. \n", num_thread);
	// DEFINE OUTPUT MATRIX
	int loci=length(POS);
	SEXP result=PROTECT(allocMatrix(REALSXP, loci, loci));
	double z=asReal(morgan_per_mb);
	int *POS_ptr; double *result_ptr; 
	POS_ptr=INTEGER(POS);
	result_ptr=REAL(result);
	#pragma omp parallel
	{
		double temp_recom_rate;
		double d;
		int temp_dist;
		#pragma omp for
		for (int i=0; i<loci; i++)
		{
			for (int j=i; j<loci; j++)
			{
				if (i==j)
				{
					result_ptr[i*loci+j]=0;
				}
				else
				{
					temp_dist=abs(POS_ptr[i]-POS_ptr[j]);
					d=(double) temp_dist*z/1e6;
					temp_recom_rate=0.5*(1-exp(-2*d));
					result_ptr[i*loci+j]=temp_recom_rate;
					result_ptr[j*loci+i]=temp_recom_rate;
				}
			}
		}
	}
	UNPROTECT(1);
	return result;
}

// CALCULATE THE CORRELATION OF CHANGE IN ALLELE FREQ IN ONE STEP. MEMORY EFFICEINT
SEXP cal_corr_matrix_unphased(SEXP genotype, SEXP POS, SEXP morgan_per_mb, SEXP N_hat, SEXP s0, SEXP st, SEXP t, SEXP ncpu)
{
	int num_thread=asInteger(ncpu);
	// OMP STUFF
	omp_set_dynamic(0);
	omp_set_num_threads(num_thread);
	Rprintf("Using %d cores. \n", num_thread);
	// DEFINE OUTPUT MATRIX
	int loci=INTEGER(GET_DIM(genotype))[1];
	int s=INTEGER(GET_DIM(genotype))[0];
	SEXP result=PROTECT(allocMatrix(REALSXP, loci, loci));
	// OTHER POINTERS
	int *genotype_ptr; int *POS_ptr; double *result_ptr;
	genotype_ptr=INTEGER(genotype); 
	POS_ptr=INTEGER(POS);
	result_ptr=REAL(result);
	double ds0=asReal(s0); double dst=asReal(st); 
	double dN=asReal(N_hat); double dt=asReal(t);
	int **x; int **y;
	x=(int **) R_alloc(num_thread, sizeof(int *));
	y=(int **) R_alloc(num_thread, sizeof(int *));
	for (int i=0; i<num_thread; i++)
	{
		x[i]=(int *) R_alloc(s, sizeof(int));
		y[i]=(int *) R_alloc(s, sizeof(int));
	}
	double temp2=0.5/ds0+1-pow((1-0.5/dN), dt)*(1-0.5/dst);
	#pragma omp parallel
	{
		int tid=omp_get_thread_num();
		double recom;
		double r; double adj_r;
		double temp1;
		#pragma omp for schedule(static, 4)
		for (int i=0; i<loci; i++)
		{
			#pragma omp simd
			for (int k=0; k<s; k++) {x[tid][k]=genotype_ptr[(unsigned long long int) i*s+k];}
			for (int j=i; j<loci; j++)
			{
				if (i==j)
				{
					result_ptr[(unsigned long long int) i*loci+j]=1;
				}
				else
				{
					#pragma omp simd
					for (int k=0; k<s; k++) {y[tid][k]=genotype_ptr[(unsigned long long int)j*s+k];}
					recom=recom_rate(POS_ptr[i], POS_ptr[j], asReal(morgan_per_mb));
					r=r_unphased_EM(x[tid], y[tid], s);
					adj_r=cal_adj_r(r, ds0);
					temp1=0.5/ds0+(1-recom)/(2*dN*recom+1-recom)-pow(1-0.5/dN, dt)*pow(1-recom, dt)*((1-recom)/(2*dN*recom+1-recom)-0.5/dst);
					result_ptr[(unsigned long long int) i*loci+j]=adj_r*temp1/temp2;
					result_ptr[(unsigned long long int) j*loci+i]=adj_r*temp1/temp2;
				}
			}
		}
	}
	UNPROTECT(1);
	return result;
}

// GENERATE Q2, SINGLE THREADED BUT MEMORY EFFICIENT
SEXP generate_Q2(SEXP adj_eigen, SEXP len)
{
	int n=asInteger(len);
	int len_eigen=length(adj_eigen);
	SEXP result=PROTECT(allocVector(REALSXP, n));
	double *adj_eigen_ptr; double *result_ptr;
	adj_eigen_ptr=REAL(adj_eigen);
	result_ptr=REAL(result);
	GetRNGstate();
	double temp;
	for (int i=0; i<n; i++)
	{
		temp=0;
		for (int j=0; j<len_eigen; j++)
		{
			temp=temp+pow(norm_rand()*sqrt(adj_eigen_ptr[j]), 2);
		}
		result_ptr[i]=temp;
	}
	PutRNGstate();
	UNPROTECT(1);
	return result;
}
