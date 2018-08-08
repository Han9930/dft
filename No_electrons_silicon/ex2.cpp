#include "MOS_function.h"
using namespace std;

double Abs(double N){ return (N<0)? -N:N;}

int main(int argc, char *argv[]){
	
	int points = atoi(argv[1]);
	double *pois = new double[points], d1 = 1e-9, d2 = 1e-9, a = 5e-9, V_n;
	// d1 is the length of first oxide layer, d2 is second, and a is length of silicon layer
	double V0 = atof(argv[2]), rho = atof(argv[3]);
	double termx = (d1+d2+a)/(points-1), termy = V0/(points-1);
	
	V_n = 0;
	Pois(points, pois, V0, rho);
	for(int i = 0; i < points; i++) cout << i*termx << ", " << Abs(test_1(termx*i, V0, rho)-*(pois+i))/test_1(termx*i, V0, rho)*100 << "%" << endl;
	return 0;
}

