#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

//*************************************
int main(){
  int N_section = 25;
  int N = 4*N_section;
  double s,ds;
  double epsilon = 0.000;
  double dt = 0.0001; // time period
  double *x,*y,*xn,*yn; // x, y are the coordinate of the boundary
  double *k; //k is the curvature
  double Step = 1500;;
  ofstream out,out_m,out_i,out_m2;
  out.open("final.dat");
  out_m.open("middle.dat");
  out_i.open("initial.dat");
  out_m2.open("middle2.dat");
  ds=1.0/(double)(N_section);
  x=new double[N+2];
  y=new double[N+2];
  xn=new double[N+2];
  yn=new double[N+2];
  k=new double[N+2];
  

  for(int i=0;i<N_section;i++){
    s=0.0+ds*(double)i;
    x[i] = s; y[i] = cos(M_PI*2.0*s);

    s=s+1.0;
    x[i+N_section] = 1.0; y[i+N_section] = 1.0-3.0*(s-1.0);

    s=s+1.0;
    x[i+2*N_section] = 3.0-s; y[i+2*N_section] = -2.0;

    s=s+1.0;
    x[i+3*N_section] = 0.0; y[i+3*N_section] = -2.0+3.0*(s-3.0);    
  }
  x[N] = x[0]; y[N] = y[0];
  x[N+1] = x[1]; y[N+1] = y[1];

  for(int i=0;i<N+1;i++){
    xn[i] = x[i];
    yn[i] = y[i];
  }
  for(int i=0;i<N+1;i++){
    x[i]=xn[N-i];
    y[i] = yn[N-i];
  }


  for(int i=0;i<N+1;i++)
    out_i<<x[i]<<"\t"<<y[i]<<"\n";


  for(int step=0;step<Step;step++){
    for( int i=1;i<N+1;i++){
      k[i] = 4.0* ( (y[i+1]-2.0*y[i]+y[i-1])*(x[i+1]-x[i-1]) - (x[i+1]-2.0*x[i]+x[i-1])*(y[i+1]-y[i-1]))/(pow((x[i+1]-x[i-1])*(x[i+1]-x[i-1])+(y[i+1]-y[i-1])*(y[i+1]-y[i-1]) ,1.5)+1e-6);
      xn[i] = x[i] + dt*(1.0-epsilon*k[i])*(y[i+1]-y[i-1])/(sqrt((x[i+1]-x[i-1])*(x[i+1]-x[i-1]) + (y[i+1]-y[i-1])*(y[i+1]-y[i-1]))+1e-6);
      yn[i] = y[i] + dt*(1.0-epsilon*k[i])*(-x[i+1]+x[i-1])/(sqrt((x[i+1]-x[i-1])*(x[i+1]-x[i-1]) + (y[i+1]-y[i-1])*(y[i+1]-y[i-1]))+1e-6);
    }
    xn[0] = xn[N]; yn[0] = yn[N];
    xn[N+1] = xn[1]; yn[N+1] = yn[1];

    for (int i = 0; i<N+2;i++){
      x[i] = xn[i]; y[i] = yn[i];
    }
    if(step ==500)
      for(int i=0;i<N+1;i++)
	out_m<<x[i]<<"\t"<<y[i]<<"\n";

    if(step ==1000)
      for(int i=0;i<N+1;i++)
	out_m2<<x[i]<<"\t"<<y[i]<<"\n";

  }

    for(int i=0;i<N+1;i++)
    out<<x[i]<<"\t"<<y[i]<<"\n";


  return 0;
}
//***************************************


