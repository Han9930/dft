//fucntions for learning DFT
#include <cmath>

#ifndef E_function
#define E_function

#define PI 3.14159265359
#define h_bar 1.05457180e-34
#define h_bar_ev (h_bar/q_e)
#define m_e 9.10938356e-31
#define q_e 1.6021766208e-19
#define es 8.854187817e-12

double E_1D(int n, double a){return n*n*PI*PI*h_bar*h_bar/2/m_e/(a*a);}
double psi_1D(int n, double x, double a){return sqrt(2/a)*sin(n*PI*x/a);}
double psi_1D_d2(int n, double x, double a){return -(double)n*n*PI*PI/(a*a)*sqrt(2/a)*sin(n*PI*x/a);}
double N_x(int n, double x, double a){return (n+1)/a-sin((n+1)*(PI*x/a))/(a*sin(PI*x/a))*cos(n*PI*x/a);}
double test_f(double x, double a){return x/a;}
double es_r(double r){return es*r;}
double test_1(double x, double V0, double rho){
	double a=5e-9, d=1e-9;
	if(x < 1e-9) return (a*q_e*rho*x/2/es_r(3.9)+V0);
	else if(x >= 1e-9 && x <= 6e-9) return -q_e*rho*(x-d)*(x-d)/2/es_r(11.68)+a*q_e*rho/2/es_r(11.68)*(x-d)+a*d*q_e*rho/2/es_r(3.9)+V0;
	else if(x > 6e-9 && x <= 7e-9) return -a*q_e*rho/2/es_r(3.9)*(x-a-2*d)+V0;
	else std::cout << "error!" << std::endl;
}

#endif
