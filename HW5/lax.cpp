#include<iostream>
#include<cmath>
#include<fstream>
#include<stdlib.h>
#include<ctime>
#include<string>
#include<sstream>

using namespace std;

void timestamp();
void initialize(int Nx, double *x, double *u,double a, double b);
void lax_wendroff_step(int Nx, double l, double* x, double* u_new, double *u_old);
string IntToStr(int n);

/*************************************************************/
int main () {
  //  timestamp();
  
  ofstream out_file;
  string filename;

  int counter=0;
  double a; double b; double dx; double dt,l;
  double *x; double *u_new; double*u_old;
  double t_init,t_end,t;
  int Nx = 2000+1;
  int Nt = 5000+1;
  
  x = new double[Nx];
  u_new = new double [Nx];
  u_old = new double [Nx];

  a=-1.0;b=1.0; t_init = 0.0; t_end = 1.0;
  dt = (t_end-t_init)/(double)(Nt-1);t= t_init;
  dx = (b-a)/(double) (Nx-1);
  l = dt/dx;
  initialize(Nx, x, u_old, a, b);
  for(int i=0;i<Nx;i++)
    u_new[i] = u_old[i];  

  while ( t<t_end){

    if(counter%5==0){
      filename ="data_" + IntToStr(counter) + ".dat";
      out_file.open(filename.c_str());
      for(int i=0;i<Nx;i++)
	out_file<<x[i]<<"\t"<<u_old[i]<<"\n";
      out_file.close(); 
    }

    lax_wendroff_step(Nx,l,x,u_new,u_old);
    
    for(int i=0;i<Nx;i++)
      u_old[i] = u_new[i];  
    
    t = t+dt;
    counter = counter +1;
  }

  timestamp();

  delete [] u_old;
  delete [] u_new;
  delete [] x;
  return 0;
}
/*************************************************************/

void timestamp() {
# define TIME_SIZE 40
  
  static char time_buffer[TIME_SIZE];
  const struct tm *tm_ptr;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm_ptr = localtime ( &now );
  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );
  cout << time_buffer << "\n";  
  return;
# undef TIME_SIZE
}
/****************************************************************************/
void initialize(int Nx, double *x, double *u,double a, double b){
  int i;
  double dx = (b-a)/(double)(Nx-1);
  for (i = 0; i<Nx; i++)
    x[i] = a + dx*(double)i;
  
  /*
  for (i = 0; i<Nx; i++){
    //   if((x[i] >= -0.4) && (x[i] <= 0.4))
    //u[i] = 1.0;
    u[i] = exp(-x[i]*x[i]*15.0);
  }
  */
  /*  
  double q,r,s;
  double pi = 3.141592653589793;
  q = 2.0 * (  1.0 - (-1.0) ) / pi;
  r = (1.0 + (-1.0) ) / 2.0;
  //
  //  S can be varied.  It is the slope of the initial condition at the midpoint.
  //
  s = 1.0;

  for ( int i = 0; i < Nx; i++ )
  {
    u[i] = - q * atan ( s * ( 2.0 * x[i] - x[0] - x[Nx-1] ) 
      / ( x[Nx-1] - x[0] ) ) + r;
  }
  */
    
  //case I
  for(i=0;i<Nx/2;i++)
    u[i] = 0.00;
  for(i=Nx/2+1;i<Nx;i++)
    u[i] = 1.00;
  

  //case II
  /*
  for(i=0;i<Nx/2;i++)
    u[i] = 1.00;
  for(i=Nx/2+1;i<Nx;i++)
    u[i] = -1.00;
  */

  return;
}
//*******************************************************************************
void lax_wendroff_step(int Nx, double l, double* x, double* un, double* uo){
  
  for ( int i=1;i<Nx-1;i++){
    un[i] = uo[i] 
      - ( l ) * 0.5*0.5*( pow ( uo[i+1], 2 ) - pow ( uo[i-1], 2 ) ) 
      +  0.5*0.5* ( l*l ) * 0.5*
      (   ( uo[i+1] + uo[i] ) * ( pow ( uo[i+1], 2 ) - pow ( uo[i], 2 ) ) 
          - ( uo[i] + uo[i-1] ) * ( pow ( uo[i], 2 ) - pow ( uo[i-1], 2 ) ) );
    
    /*
    un[i] = u[i] - l*0.5*0.5*(pow(u[i+1],2)-pow(u[i-1],2))
           + pow(l,2)*0.5*0.5*0.5*((u[i]+u[i+1])*(pow(u[i+1],2)-pow(u[i],2))-
				   (u[i]+u[i-1])*(pow(u[i],2)-pow(u[i-1],2)));
    */
  }


}

//*****************************************************************************
string IntToStr(int n){
  stringstream result;
  result<<n;
  return result.str();
}
//*****************************************************************************
