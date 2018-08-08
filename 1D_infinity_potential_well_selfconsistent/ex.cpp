#include <iostream>
#include <cstdlib>
#include "Self_Consistent_functions.h"

using namespace std;

int main(int argc, char* argv[]){

	if( argc != 4 ) {cout << "There should be 3 arguments" << endl; return 1;}
	FILE *fp;
	char filename[30];
	int m, R, M = atoi(argv[3]), N = atoi(argv[1]), points = atoi(argv[2]);

	double *sch = new double[points*points], *n = new double[points], *E = new double[points], *pois = new double[points], a = 5e-9, term = a/(points-1);
	m = Pois(N, points, pois);
	for(int i = 0; i < M; i++){
		sprintf(filename, "N%d_P%d_R%d_Nx.txt", N, points, i);
		fp = fopen(filename, "w");
		m = R_Sch(N, points, pois, sch, E);
		R = R_Nx(N, points, sch, n, E);
		R = R_Pois(N, points, pois, n);
		for(int j = 0; j < points; j++) fprintf(fp, "%e, %e\n", j*term, *(n+j));
		fclose(fp);
	}
	
	m = R_Sch(N, points, pois, sch, E);
	for(int i = 0; i < points-2; i++) cout << (i+1)*term << ", " << *(sch+i+points-2) << endl;
	return 0;
}

