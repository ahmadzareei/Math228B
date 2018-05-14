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
void glimms_step(int Nx, double l, double* x, double* u_new, double *u_old);
string IntToStr(int n);
void error(double *u,double *x, int Nx){
  double err=0.0;
  double * exact;
  exact = new double [Nx];
  /*
  //for case I
  for(int i=0;i<Nx;i++){
    if(x[i]<0.5)
      exact[i] = 1;
    else
      exact[i] = 0.0;
  }
  */

  //for case II
  for(int i=0;i<Nx;i++){
    if(x[i]<=0) 
      exact[i] = 0;
    else if(x[i] >1.0)
      exact[i] = 1.0;
    else
      exact[i] = x[i];
  }
  for(int i=0;i<Nx;i++)
    err += abs(exact[i] - u[i]);

  err = err /(double)(Nx);

  cout<<"Error is:"<<err<<"\n";
}

/*************************************************************/
int main () {
  //  timestamp();
  
  ofstream out_file;
  string filename;

  int counter=0;
  double a; double b; double dx; double dt,l;
  double *x; double *u_new; double*u_old;
  double t_init,t_end,t;
  int Nx = 1000+1;
  int Nt = 5000+1;
  
  x = new double[Nx];
  u_new = new double [Nx];
  u_old = new double [Nx];


  srand(time(NULL));

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

    glimms_step(Nx,l,x,u_new,u_old);
    
    for(int i=0;i<Nx;i++)
      u_old[i] = u_new[i];  
    
    t = t+dt;
    counter = counter +1;
  }

  //  error(u_new, x, Nx);

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

  
  //exp case
  for (i = 0; i<Nx; i++){
    u[i] = exp(-x[i]*x[i]*15.0);
  }
  

  /*
  //case I
  for(i=0;i<Nx/2;i++)
    u[i] = 1.00;
  for(i=Nx/2+1;i<Nx;i++)
    u[i] =0.00;
  */

  /*  
  //case II
  for(i=0;i<Nx/2;i++)
    u[i] = 0.00;
  for(i=Nx/2+1;i<Nx;i++)
    u[i] =1.00;
  */
  return;
}
//*******************************************************************************
void glimms_step(int Nx, double l, double* x, double* un, double* u){
  double ur=0.0;
  double ul=0.0;
  double result;
  double random;
  double x0,s;

  for ( int i=1;i<Nx-1;i++){
    random = ((double) rand() / (RAND_MAX))-0.5;
    if(random<0){
      ul=u[i-1];
      ur=u[i];
      x0=-0.5;
    }else{
      ul=u[i];
      ur=u[i+1];
      x0=0.5;
    }
    s=(ul+ur)/2;
    if (ul>ur){
      if((random-x0)/l<s){
	un[i]=ul;
      }else{
	un[i]=ur;
      }
    }else{
      if(4*(random-x0)/l<ul){
	un[i]=ul;
      }else if (4*(random-x0)/l>ur){
	un[i]=ur;
      }else {
	un[i]=4*(random-x0)/l;
      } 
    }
  }
  
}

//*****************************************************************************
string IntToStr(int n){
  stringstream result;
  result<<n;
  return result.str();
}
//*****************************************************************************
