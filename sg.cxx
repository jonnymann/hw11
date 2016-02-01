#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);

void step(cmplx* v0, cmplx* v1, const double dx, const int Nx, const double xmin, const double dt,const double k){
    cout << "dt=" << dt << '\t' << "Nx=" << Nx  << '\t' << "dx=" << dx << endl;
    double x;
    cmplx ds;
    cmplx as = cmplx(0,+dt/(4*dx*dx));
    cmplx temp1 = 0;
    cmplx temp2 = v0[0];
   
    for(int i=0;i<Nx-1;i++) // Berechnung des Vektors (A*) mal v1
    {
        
        x=xmin+i*dx;
        ds=cmplx(1,-dt/(2*dx*dx)-dt*k/4*x*x);
        v0[i]=as*temp1+ds*temp2+as*v0[i+1];
        temp1=temp2;
        temp2=v0[i+1];
        
    }
    x=x+dx;
    ds=cmplx(1,-dt/(2*dx*dx)-dt*k/4*x*x);
    v0[Nx-1]=as*temp1+ds*temp2;
    
    cmplx* d = new cmplx[Nx];    //Diagnole
    cmplx* au = new cmplx[Nx];  // Obere Diagonalenparallele
    cmplx* al = new cmplx[Nx];  // Untere Diagonalenparallele
    
    for(int i=0; i<Nx;i++) // Initialisierung der Matrix A
    {
        x=xmin+i*dx;
        
        d[i]=cmplx(1,dt/(2*dx*dx)+dt*k/4*x*x);
        au[i]=cmplx(0,-dt/(4*dx*dx));
        al[i]=cmplx(0,-dt/(4*dx*dx));
        
    }
    
    for(int i=0;i<Nx-1;i++) //Transformation in obere Dreiecksmatrix, die Nummerierung ist Spaltenweise
    {
        d[i+1]=d[i+1]-al[i]/d[i]*au[i+1];
        v0[i+1]=v0[i+1]-al[i]/d[i]*v0[i];
    }
    
    v1[Nx-1]=v0[Nx-1]/d[Nx-1]; //Rückeinsetzen zur Bestimmung von v1
    for(int i=Nx-2;i>-1; i--)
    {
        //cout << v0[i] << '\t';
        v1[i]=(v0[i]-au[i+1]*v1[i+1])/d[i];
        //cout << v1[i] << endl;
    }
    delete[] d;
    delete[] au;
    delete[] al;
};
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40.;
        const double xmax = 40.;
	const double Tend = 12*3.1415;
	const double dx = (xmax-xmin)/Nx;
	const double dt = 0.1;
        double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);


	const double lambda = 10.;
        const double omega = 0.2; //-> k = 0.04
        const double k = pow(omega,2);
        const double alpha = sqrt(0.2);

  stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
        cmplx* psi1 = new cmplx[Nx];
        cmplx* temp;

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0",dx,Nx,xmin,alpha,lambda, omega,t);
        
        
	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
                    step(psi0, psi1,dx,Nx,xmin,dt,k);
                    temp = psi0;
                    psi0 = psi1;
                    psi1 = temp;
                    t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;

	return 0;
}
//-----------------------------------

//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
            x = xmin + i * dx;
            xi = alpha * x;
            xil = alpha * lambda;
            h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}
