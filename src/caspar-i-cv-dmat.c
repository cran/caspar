#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>`
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Applic.h>
#include <R_ext/RS.h>

/*
ARGUMENTs TO CASPAR:
double *Xpass     : n by p matrix of data
                    NOTE: the first column is assumed to be an intercept
                          and is included in the model always
double *Y         : n vector of response
double *Xtestpass : nt by p matrix of test data
double *Ytest     : nt vector of test response
int *npass        : n
int *ppass        : p
int *ntestpass    : nt
double *plocs     : p by p matrix of predictor distance
double *palpha    : alpha (mixing parameter)
double *ph        : h     (bandwidth parameter)
int *psequ        : RETURN VALUE -- sequence of predictors added
double *psse      : RETURN VALUE -- SSE for test data on predictors
int *ptrace       : logical, print some things?
*/

void caspar(double *Xpass, double *Y, double *Xtestpass, double *Ytest, int *npass, int *ppass, int *ntestpass, double *plocs, double *palpha, double *ph, int *psequ, double *psse, int *ptrace);
double corel(int len, double *X, double *Y);
double depan(double d);

int whichmax(int len, double *vec);
int inintvec(int len, int look, int *vec);

int ivecSum(int *vec, int len);
double dvecSum(double *vec, int len);

void MatrixAdd(double **A, double **B, int n, int m, double **R);
void MatrixCMul(double **A, double c, int n, int m);
void MatrixCAdd(double **A, double c, int n, int m);

double mnnormden(double *x, double covdet, double **covinv, int p);

void pmat(double **mat, int nr, int nc);
void pmati(int **mat, int rows, int cols);
void rpmat(double **mat, int rows, int cols);
double inner(double *A, double *B, int len);

void MatrixMultiply(double **A, double **B, int n, int m, int l, double **R);
void MatrixTMultiplyA(double **A, double **B, int n, int m, int l, double **R);
void MatrixTMultiplyB(double **A, double **B, int n, int m, int l, double **R);

// FORTRAN MATRIX INVERSION

void MatInvF(double **mat, int n, double **mati, int *pcode);

void caspar(double *Xpass, double *Y, double *Xtestpass, double *Ytest, int *npass, int *ppass, int *ntestpass, double *plocs, double *palpha, double *ph, int *psequ, double *psse, int *ptrace){
  double **X, **XT, **Xtest, **gram, **TG, **TGI, **BT, *locs, *kern;
  double h = *ph, alpha = *palpha, dstep = 1.0;
  int n = *npass, p = *ppass, nt = *ntestpass;
  int i,j, *numneigh, ncount, **neigh;
  int *pcode; int verbose = *ptrace;

  pcode = (int *)calloc(1,sizeof(int));

  if(alpha == 0.0) alpha = .0001;

  X = (double **)calloc(n, sizeof(double *));
  for(i=0;i<n;i++) X[i] = (double *)calloc(p,sizeof(double));

  Xtest = (double **)calloc(nt, sizeof(double *));
  for(i=0;i<nt;i++) Xtest[i] = (double *)calloc(p,sizeof(double));

  XT = (double **)calloc(p, sizeof(double *));
  for(i=0;i<p;i++) XT[i] = (double *)calloc(n,sizeof(double));

  gram = (double **)calloc(p, sizeof(double *));
  for(i=0;i<p;i++) gram[i] = (double *)calloc(p,sizeof(double));

  TG = (double **)calloc(p, sizeof(double *));
  for(i=0;i<p;i++) TG[i] = (double *)calloc(p,sizeof(double));

  TGI = (double **)calloc(p, sizeof(double *));
  for(i=0;i<p;i++) TGI[i] = (double *)calloc(p,sizeof(double));

  BT = (double **)calloc(p, sizeof(double *));
  for(i=0;i<p;i++) BT[i] = (double *)calloc(n,sizeof(double));

  neigh = (int **)calloc(p, sizeof(int *));
  for(i=0;i<p;i++) neigh[i] = (int *)calloc(p,sizeof(int));

  locs = (double *)calloc(p*p, sizeof(double));
  kern = (double *)calloc(p, sizeof(double));
  for(i=0;i<p;i++) kern[i] = alpha;
  for(i=0;i<(p*p);i++) locs[i] = plocs[i];
  numneigh = (int *)calloc(p, sizeof(int));


  for(i=0;i<n;i++)
    for(j=0;j<p;j++)
      X[i][j] = Xpass[j*n + i];
  for(i=0;i<nt;i++)
    for(j=0;j<p;j++)
      Xtest[i][j] = Xtestpass[j*nt + i];
  for(i=0;i<n;i++)
    for(j=0;j<p;j++)
      XT[j][i] = X[i][j];

  MatrixTMultiplyB(X,X,p,n,p,gram);

  int step = 0;
  int takestep;
  double **tempdes;
  double **tempdestest;

  tempdes = (double **)calloc(n, sizeof(double *));
  for(i=0;i<n;i++) tempdes[i] = (double *)calloc(p, sizeof(double));
  tempdestest = (double **)calloc(nt, sizeof(double *));
  for(i=0;i<nt;i++) tempdestest[i] = (double *)calloc(p, sizeof(double));

  int *sequ; 
  sequ = (int *)calloc(p, sizeof(int));
  for(i=0;i<p;i++) sequ[i] = -1;
  double *corz; 
  corz = (double *)calloc(p, sizeof(double));
  double *cres;
  cres = (double *)calloc(n, sizeof(double));
  double *crestest;
  crestest = (double *)calloc(nt, sizeof(double));
  double *beta;
  beta = (double *)calloc(p, sizeof(double));
  double deter;

  for(i=0;i<n;i++) cres[i] = Y[i];
  for(i=0;i<nt;i++) crestest[i] = Ytest[i];

  int cutalg = p;
  if(p > n) cutalg = n;

  while(step < cutalg){
    j=0;

    for(i=0;i<p;i++){
      if(inintvec(step, i, sequ) == 0){
	corz[i] = kern[i] * fabs(corel(n, cres, XT[i]));
      }
      else{
	corz[i] = 0.0;
      }
    }

    takestep = whichmax(p,corz);

    if(step == 0) takestep = 0;

    if(verbose == 1){
      Rprintf(" [%d : %d | %d] ", step, n, p);
      Rprintf("STEP: %d",takestep+1);
    }
    sequ[step]  = takestep;
    psequ[step] = takestep;
    
    if(step > 0){
      for(i=0;i<p;i++) 
	if(kern[i] != alpha && step >0) 
	  kern[i] = (kern[i] - alpha)*(dstep-1.0)/(dstep) + alpha;
      //for(i=0;i<p;i++) kern[i] += (1.0-alpha)*depan(fabs(locs[i] - locs[sequ[step]])/h)/dstep;
      for(i=0;i<p;i++) 
	kern[i] += (1.0-alpha)*depan(locs[sequ[step]*p + i]/h)/dstep;
      /*
      for(i=0;i<p;i++) 
	Rprintf("(%g)",kern[i]);
      */
    }


    for(i=0;i<n;i++)
      tempdes[i][step] = X[i][sequ[step]];
    for(i=0;i<nt;i++)
      tempdestest[i][step] = Xtest[i][sequ[step]];

    for(i=0;i<(step+1);i++)
      for(j=0;j<(step+1);j++)
	TG[i][j] = gram[sequ[i]][sequ[j]];

    MatInvF(TG,step+1,TGI,pcode);
    if(verbose == 1) Rprintf(" [I: %d]\n",pcode[0]);

    if(step + 1 == 1) TGI[0][0] = 1.0 / TG[0][0];

    MatrixTMultiplyA(TGI,tempdes,step+1,step+1,n,BT);
    
    for(i=0;i<(step+1);i++){
      beta[i] = inner(BT[i],Y,n);
    }

    for(i=0;i<n;i++)
      cres[i] = Y[i] - inner(tempdes[i],beta,step+1);

    for(i=0;i<nt;i++)
      crestest[i] = Ytest[i] - inner(tempdestest[i],beta,step+1);

    for(i=0;i<nt;i++)
      psse[step] += crestest[i]*crestest[i];
    
    psse[step] /= nt;
    
    step ++;
    dstep = dstep + 1.0;
  }

  for(i=0;i<n;i++) free(X[i]);
  free(X);
  for(i=0;i<p;i++) free(XT[i]);
  free(XT);
  for(i=0;i<nt;i++) free(Xtest[i]);
  free(Xtest);
  for(i=0;i<p;i++) free(gram[i]);
  free(gram);
  for(i=0;i<p;i++) free(TG[i]);
  free(TG);
  for(i=0;i<p;i++) free(TGI[i]);
  free(TGI);
  for(i=0;i<p;i++) free(BT[i]);
  free(BT);
  for(i=0;i<p;i++) free(neigh[i]);
  free(neigh);
  for(i=0;i<n;i++) free(tempdes[i]);
  free(tempdes);
  for(i=0;i<nt;i++) free(tempdestest[i]);
  free(tempdestest);

  free(locs);
  free(kern);
  free(numneigh);
  free(sequ);
  free(corz);
  free(cres);
  free(crestest);
  free(beta);
}

double corel(int len, double *X, double *Y){
  int i;
  double Xsum = 0.0, Ysum = 0.0, XSsum = 0.0, YSsum = 0.0, XYsum = 0.0;
  for(i=0;i<len;i++){
    Xsum  += X[i];
    Ysum  += Y[i];
    XSsum += X[i]*X[i];
    YSsum += Y[i]*Y[i];
    XYsum += X[i]*Y[i];
  }
  
  double cor = (len * XYsum - Xsum*Ysum) / (sqrt( len * XSsum - Xsum*Xsum) * sqrt( len *  YSsum - Ysum*Ysum));

  return(cor);
}

double depan(double d){
  if( d > 1) return(0);
  else return((3.0/4.0)*(1 - d*d));
}

int inintvec(int len, int look, int *vec){
  int i;
  int ret = 0;
  for(i=0;i<len;i++) if(vec[i] == look) ret = 1;
  return(ret);
}

int whichmax(int len, double *vec){
  int i;
  int ret = 0;
  for(i=1;i<len;i++)
    if(vec[i] > vec[ret])
      ret = i;
  return(ret);
}

int ivecSum(int *vec, int len){
  int i, ret = 0;
  for(i=0;i<len;i++) ret += vec[i];
  return(ret);
}

double dvecSum(double *vec, int len){
  int i;
  double ret = 0.0;
  for(i=0;i<len;i++) ret += vec[i];
  return(ret);
}

void pmat(double **mat, int rows, int cols){
  int i, j;

  for(i = 0; i< rows; i++){
    for(j = 0; j< cols; j++){
      printf("%g ",mat[i][j],i,j);
    }
    printf("\n");
  }
}

void pmati(int **mat, int rows, int cols){
  int i, j;

  for(i = 0; i< rows; i++){
    for(j = 0; j< cols; j++){
      printf("%d ",mat[i][j],i,j);
    }
    printf("\n");
  }
}


void rpmat(double **mat, int rows, int cols){
  int i, j;

  for(i = 0; i< rows; i++){
    for(j = 0; j< cols; j++){
      printf("%g ",round(mat[i][j]*100.0));
    }
    printf("\n");
  }
}

// multiplies n X m matrix A
//       with m x l matrix B
// stored in n x l  matrix R, which is all zeroes to start (use calloc)

void MatrixMultiply(double **A, double **B, int n, int m, int l, double **R){
  int i,j,k;
  for(i = 0; i<n; i++)
    for(j = 0; j<l; j++)
      for(k = 0; k<m; k++)
	R[i][j] += A[i][k]*B[k][j];
}

// for A times B^T

void MatrixTMultiplyA(double **A, double **B, int n, int m, int l, double **R){
  int i,j,k;
  for(i = 0; i<n; i++)
    for(j = 0; j<l; j++){
      R[i][j] = 0.0;
      for(k = 0; k<m; k++)
	R[i][j] += A[i][k]*B[j][k];
    }
}

// for A^T times B

void MatrixTMultiplyB(double **A, double **B, int n, int m, int l, double **R){
  int i,j,k;
  for(i = 0; i<n; i++)
    for(j = 0; j<l; j++){
      R[i][j] = 0.0;
      for(k = 0; k<m; k++)
	R[i][j] += A[k][i]*B[k][j];
    }
}

void MatrixAdd(double **A, double **B, int n, int m, double **R){
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      R[i][j] = A[i][j] + B[i][j];
}

void MatrixCMul(double **A, double c, int n, int m){
  int i, j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      A[i][j] *= c;
}

void MatrixCAdd(double **A, double c, int n, int m){
  int i, j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      A[i][j] += c;
}

double inner(double *A, double *B, int len){
  int i;
  double ret = 0.0;
  for(i = 0; i<len; i++) ret += B[i]*A[i];
  return(ret);
}

// mvnorm density
// assume mean 0 vector

double mnnormden(double *x, double covdet, double **covinv, int p){
  int i,j;

  double cterm = 1.0 / (sqrt(covdet) * pow(M_PI,p/2.0));
  double exponent = 0.0;    
  double termtemp, ret;
  
  for(i=0;i<p;i++){
    termtemp = 0.0;
    for(j=0;j<p;j++)
      termtemp += (x[j]) * covinv[i][j];
    exponent += termtemp * (x[i]);
  }

  exponent *= -0.5;
  ret = exp(exponent);
  ret *= cterm;
  return(ret);
}

// Matrix inversion using LAPACK routine

void MatInvF(double **mat, int n, double **mati, int *pcode){
  double *APTF, *BPTF;
  int *pivot, code, i, j;

  APTF  = (double *) calloc(n*n, sizeof(double));
  BPTF  = (double *) calloc(n*n, sizeof(double));
  pivot = (int *)calloc(n  , sizeof(int));

  for(i=0;i<n;i++)
    for(j=0;j<n;j++){
      APTF[j*n + i] = mat[i][j];
      if(i == j) BPTF[j*n + i] = 1.0;
    }

  dgesv_(&n,&n,APTF,&n,pivot,BPTF,&n,&code);

  //for(i=0;i<n;i++) Rprintf("%d ", pivot[i]);

  for(i=0;i<n;i++)
    for(j=0;j<n;j++){
      mati[i][j] = BPTF[j*n + i];
    }

  //if(code != 0) Rprintf(" INVERSION PROBLEM ",code);

  pcode[0] = code;

  free(APTF);
  free(BPTF);
  free(pivot);
}


