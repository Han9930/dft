#include "MOS_function.h"
using namespace std;

int main(int argc, char *argv[]){
	
	int points = atoi(argv[1]);
	double *pois = new double[points], d1 = 1e-9, d2 = 1e-9, a = 5e-9, V_n;
	double V0 = atof(argv[2]), rho = atof(argv[3]);
	double termx = (d1+d2+a)/(points-1), termy = V0/(points-1);
	
	V_n = 0;
	Pois(points, pois, V0, rho);
	
	for(int j = 0; j < points; j++){
		Pois(points, pois, V_n, rho);
		for(int i = 0; i < points; i++) cout << i*termx <<  ", " << V_n << ", " << *(pois+i) << endl;
		V_n += termy;
	}

//	for(int i = 0; i < points; i++) cout << i*termx << ", " << *(pois+i) << endl;
	return 0;
}

// program for calculation of Poisson's equation
