#include <iostream>
#include <cstdio>
#include "Constant.h"

#ifndef MOSfunctions
#define MOSfunctions

using namespace std;

extern "C" {
	void dgetrf_(int* M, int* N, double* A, int* LDA, int* IPIV, int* INFO);
	void dgetrs_(char* TRANS, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
	void dgetri_(int* N, double* A, int* LDA, int* IPIV, double* WORK, int* LWORK, int* INFO);
} // refer to lapack library


int Solver(int N, double *A, double *b){

	int i, j, *ipiv = new int[N], nrhs = 1, info, lwork = 4*N;
        char trans = 'N';
	double *work = new double[lwork], *temp = new double[N];

	dgetrf_(&N, &N, A, &N, ipiv, &info);
	dgetri_(&N, A, &N, ipiv, work, &lwork, &info);

	for(i = 0; i < N; i++) {*(temp+i) = *(b+i); *(b+i) = 0;}
        for(i = 0; i < N; i++) for(j = 0; j < N; j++) *(b+i) += A[N*i+j]*temp[j];

	delete ipiv, work, temp;
        return info;
} // Linear Solver


int Pois(int points, double *pois, double V0, double rho){

	int i, j, info;
	double *temp1 = new double[points*points], d1 = 1e-9, d2 = 1e-9, a = 5e-9, term = (d1+d2+a)/(points-1);
	double L, R, dist1L, dist1R, dist2L, dist2R;
	rho *= q_e;
	for(i = 0; i < points; i++){
		if(i == 0 || i == points-1) *(pois+i) = V0;
		else if(term*i < d1 || term*i > a+d1) *(pois+i) = 0;
		if(term*i > d1 && term*i < a+d1) *(pois+i) = -rho/es_r(11.68);

		if(term*i < d1 && term*(i+1) > d1) { L=i; dist1L = d1 - term*i; dist1R = term*(i+1) - d1; }
		if(term*i < a+d1 && term*(i+1) > a+d1) { R=i; dist2L = (d1+a) - term*i; dist2R = term*(i+1) - (d1+a); }

	}

	for(i = 0; i < points; i++){
		for(j = 0; j < points; j++) *(temp1+i*points+j) = 0;
		if(i == 0 || i == points-1) {*(temp1+i*points+i) = 1; continue;}
		else if(i != L && i != L+1 && i != R+1 && i != R){
			*(temp1+points*i+i) = -2/term/term;
			*(temp1+points*i+i-1) = 1/term/term;
			*(temp1+points*i+i+1) = 1/term/term;
		}
		else if(i == L) {
			*(temp1+points*i+i-1) = 1/term; 
			*(temp1+points*i+i) = -1/term-1/dist1L+1/dist1L*dist1R*es_r(3.9)/(dist1R*es_r(3.9)+dist1L*es_r(11.68));
			*(temp1+points*i+i+1) = es_r(11.68)/(dist1R*es_r(3.9)+dist1L*es_r(11.68));
		}
		else if(i == L+1) {
			*(temp1+points*i+i-1) = es_r(3.9)/(dist1R*es_r(3.9)+dist1L*es_r(11.68))/(dist1R/2+term/2);
			*(temp1+points*i+i) = (-1/term-1/dist1R+1/dist1R*dist1L*es_r(11.68)/(dist1R*es_r(3.9)+dist1L*es_r(11.68)))/(dist1R/2+term/2);
			*(temp1+points*i+i+1) = 1/term/(dist1R/2+term/2);
		}
		
		else if(i == R) {
			*(temp1+points*i+i-1) = 1/term/(dist2L/2+term/2); 
			*(temp1+points*i+i) = (-1/term-1/dist2L+1/dist2L*dist2R*es_r(11.68)/(dist2R*es_r(11.68)+dist2L*es_r(3.9)))/(dist2L/2+term/2);
			*(temp1+points*i+i+1) = es_r(3.9)/(dist2R*es_r(11.68)+dist2L*es_r(3.9))/(dist2L/2+term/2);
		}
		else if(i == R+1) {
			*(temp1+points*i+i-1) = es_r(11.68)/(dist2R*es_r(11.68)+dist2L*es_r(3.9));
			*(temp1+points*i+i) = -1/term-1/dist2R+1/dist2R*dist2L*es_r(3.9)/(dist2R*es_r(11.68)+dist2L*es_r(3.9));
			*(temp1+points*i+i+1) = 1/term;
		}
	}
	
	info = Solver(points, temp1, pois);

	delete temp1;

	return (info == 0)? true:false;
} // Solve Poisson's equation to find out electrostatic potential

#endif
