#include <iostream>
#include <cstdio>
#include "Constant.h"

#ifndef Nfunctions
#define Nfunctions

using namespace std;

extern "C" void dsyevx_(char* JOBZ, char* RANGE, char* UPLO, int* N, double* A, int* LDA, double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W, double* Z, int* LDZ, double* WORK, int* LWORK, int* IWORK, int* IFAIL, int *INFO);
extern "C" void dgetrf_(int* M, int* N, double* A, int* LDA, int* IPIV, int* INFO);
extern "C" void dgetrs_(char* TRANS, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
extern "C" void dgetri_(int* N, double* A, int* LDA, int* IPIV, double* WORK, int* LWORK, int* INFO);
extern "C" void dgesvx_(char* FACT, char* TRANS, int* N, int* NRHS, double* AB, int* LDAB, double* AFB, int* LDAFB , int* IPIV, char* EQUED, double* R, double* C, double* B, int* LDB, double* X, int* LDX, double* RCOND, double* FERR, double* BERR, double* WORK, int* IWORK, int* INFO);
// Reffer to lapack library

void setMatA1(int n, double *A); // Set matrix for 1D infinity potential well
int Pois(int N, int points, double *pois); //Poisson's equation with electron density of 1D inf well
int psix(int N, int DN, double *psi, double *E); // N x DN matrix for psi(x)
int Nx(int N, int DN, double *n); // DN matrix for N(x)


int R_sch(int N, int points, double *sch, double *pois, double *E); //Poisson equation with N(x) of Schrodinger equ.
int R_Nx(int N, int points, double *psi, double *n, double *E); //Calculate electron density
int R_Pois(int N, int points, double *pois, double *n); //Calculate electrostatic potential with electron density

int Solver3(int N, double *A, double *b); //Solving Ax=b

int Solver3(int N, double *A, double *b){
	int i, j, *ipiv = new int[N], nrhs = 1, info, lwork = 4*N;
	char trans = 'N';
	double *work = new double[lwork], *temp = new double[N];
	
	dgetrf_(&N, &N, A, &N, ipiv, &info);
	dgetri_(&N, A, &N, ipiv, work, &lwork, &info);
	
	for(i = 0; i < N; i++) {*(temp+i) = *(b+i); *(b+i) = 0;}
	for(i = 0; i < N; i++) for(j = 0; j < N; j++) *(b+i) += A[N*i+j]*temp[j];

	delete ipiv, work, temp;
	return info;
}



int Solver(int N, double *A, double *b){
	int i, j, *ipiv = new int[N], nrhs = 1, info;
	char trans = 'N';
	
	dgetrf_(&N, &N, A, &N, ipiv, &info);
	for(i = 0; i < N; i++) {for(j = 0; j < N; j++) cout << A[N*i+j] << " "; cout << endl;}
	dgetrs_(&trans, &N, &nrhs, A, &N, ipiv, b, &N, &info);

	delete ipiv;
	return info;
}


void setMatA1(int n, double *A){

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) A[n*i+j] = 0;

		A[n*i+i] += 2;

		if(!i)  A[n*i+1] += -1;
		else if(i==(n-1)) A[n*i+n-2] += -1;
		else{
			A[n*i+i-1] += -1;
			A[n*i+i+1] += -1;
		}
	}
}

int psix(int N, int DN, double *psi, double *E){
	
	double a = 5e-9, sum = 0, Eigen = N*N*PI*PI/(DN+1)/(DN+1), term = a/(DN+1);
	int judge = 0;

	double *A;
        A = new double[DN*DN];

	setMatA1(DN, A);
	
	char jobz = 'V', range = 'V', uplo = 'U';
	int n, lda, info, ldz, ldvr, lwork, il, iu=N, m;
	double vl=0, vu=Eigen;
	n=lda=ldz=DN;

	lwork=8*DN;
	int *iwork = new int[DN*5];
	int *ifail = new int[DN];
	double *w = new double[DN];
	double *work = new double[lwork];

	dsyevx_(&jobz, &range, &uplo, &n, A, &lda, &vl, &vu, &il, &iu, &vl, &m, w, psi, &ldz, work, &lwork, iwork, ifail, &info);

	for(int j = 0; j < N; j++){
  		for(int i = 0; i < DN; i++) sum += psi[DN*j+i]*psi[DN*j+i]*term;
		for(int i = 0; i < DN; i++) psi[DN*j+i] /= sqrt(sum);
		sum = 0;
	}

	for(int i = 0; i < N; i++){ judge += (ifail[i] || (E <= 0)); *(E+i) = w[i]*(h_bar/(2*m_e*term*term)*h_bar);}

	
	delete A, iwork, ifail, w, work;
	return judge;

}

int Nx(int N, int DN, double *n){
	int T, i, j;
	double *temp1 = new double[DN], *temp2 = new double[DN*DN], temp3;

	T = psix(N, DN, temp2, temp1);
	for(i = 0; i < DN; i++) for(j = 0; j < i; j++){
		temp3 = *(temp2+DN*i+j);
		*(temp2+DN*i+j) = *(temp2+DN*j+i);
		*(temp2+DN*j+i) = temp3;
	}

	for(i = 0; i < DN; i++) *(n+i) = 0;
	for(i = 0; i < DN; i++) for(j = 0; j < N; j++) *(n+i) += (*(temp2+DN*i+j))*(*(temp2+DN*i+j));
	delete temp1, temp2;
	return T;
}
// The above code is used when the potential is zero.

int Pois(int N, int points, double *pois){

	int i, j, info;
	double *temp1 = new double[points*points], a = 5e-9, term = a/(points-1);

	int T = Nx(N, points-2, pois+1);
	*(pois) = 0;
	*(pois+points-1) = 0;

	for(i = 0; i < points; i++){
	       for(j = 0; j < points; j++) *(temp1+i*points+j) = 0;
	       
	       if(!i) {*(temp1+i*points) = 1; continue;}
	       if(i == points -1) {*(temp1+i*points+points-1) = 1; continue;}

	       *(temp1+points*i+i) = -2;
	       *(temp1+points*i+i-1) = 1;
	       *(temp1+points*i+i+1) = 1;
	}

	info = Solver3(points, temp1, pois);

	for(i=1; i<points-1; i++) *(pois+i) = -(*(pois+i))*q_e/es_r(11.68)*term*term/1e-19;
	delete temp1;

	return (info == 0)? true:false;
}

int R_Sch(int N, int points, double *pois, double *sch, double *E){

	int i, j, k, R = 0, judge = 0;

	double *temp1 = new double[(points-2)*(points-2)], sum = 0, term = 5e-9/(points-1);

	for(i = 0; i < points-2; i++){
		for(j = 0; j < points-2; j++) *(temp1+i*(points-2)+j) = 0;
		*(temp1+(points-2)*i+i) = 2 + *(pois+i+1)*q_e/h_bar*2*m_e/h_bar*term*term;
		if(i != 0) *(temp1+(points-2)*i+i-1) = -1;
		if(i != points-3) *(temp1+(points-2)*i+i+1) = -1;

	}


	char jobz = 'V', range = 'V', uplo = 'U';
	int n, lda, info, ldz, ldvr, lwork, il=1, iu=N, m;
	double vl=0, vu=1e30;
	n=lda=ldz=points-2;

	lwork=8*(points-2);
	int *iwork = new int[(points-2)*5];
	int *ifail = new int[points-2];
	double *w = new double[points-2];
	double *work = new double[lwork];

	dsyevx_(&jobz, &range, &uplo, &n, temp1, &lda, &vl, &vu, &il, &iu, &vl, &m, w, sch, &ldz, work, &lwork, iwork, ifail, &info);

	for(j = 0; j < points-2;j++){
		for(i = 0; i < points-2; i++) sum += sch[(points-2)*j+i]*sch[(points-2)*j+i]*term;
		for(i = 0; i < points-2; i++) sch[(points-2)*j+i] /= sqrt(sum);
		sum = 0;
	}

	for(i = 0; i < points-2; i++){ judge += (ifail[i] || (E <= 0)); *(E+i) = w[i]*(h_bar/(2*m_e*term*term)*h_bar);}
	
	delete temp1, iwork, ifail, w, work;
	return m;
}


int R_Nx(int N, int points, double *psi, double *n, double *E){
	int i, j, repeat = N;
	double temp3;

	for(i = 0; i < points; i++) *(n+i) = 0;
	for(i = 0; i < N; i++) for(j = 0; j < points-2; j++) *(n+j+1) += (*(psi+(points-2)*i+j))*(*(psi+(points-2)*i+j));
	
	return repeat;
}

int R_Pois(int N, int points, double *pois, double *n){

	int i, j, info;
	double *temp1 = new double[points*points], a = 5e-9, term = a/(points-1);

	for(i = 0; i < points; i++) *(pois+i) = *(n+i);
	for(i = 0; i < points; i++){
	       for(j = 0; j < points; j++) *(temp1+i*points+j) = 0;
	       
	       if(!i) {*(temp1+i*points) = 1; continue;}
	       if(i == points -1) {*(temp1+i*points+points-1) = 1; continue;}

	       *(temp1+points*i+i) = -2;
	       *(temp1+points*i+i-1) = 1;
	       *(temp1+points*i+i+1) = 1;
	}

	info = Solver3(points, temp1, pois);

	for(i=1; i<points-1; i++) *(pois+i) = -(*(pois+i))*q_e/es_r(11.68)*term*term/1e-19;
	delete temp1;

	return (info == 0)? true:false;
}
// This code is used for self-consistent

#endif
