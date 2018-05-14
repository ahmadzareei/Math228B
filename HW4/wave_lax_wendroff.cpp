#include<iostream>
#include<fstream>
#include<cmath>
#include<stdlib.h>
#include<string>
#include<sstream>

using namespace std;

string IntToStr(int n){
  stringstream result;
  result<<n;
  return result.str();
}


void Lax_Wandroff_Step(double*w1, double* w2, double* w1_new, double* w2_new,double dt,double dx,double *a,int N);
double Energy_Compute(double *w1,double *w2,double * u, double *a, double dt,double dx,int N);
void Integrator(double *w1, double *w2, double * u,double *a,  int N);

int main(){
  //output files:
  ofstream out,out_energy;
  out.open("initial_shape.dat");
  out_energy.open("energy_time.dat");
  
  ofstream test; test.open("test.dat");

  //defining variables
  double *x,*w1,*a,*w1_new,*w2,*w2_new,*u;
  int N=2048;
  
  double t = 0;
  double T = 10;
  double dx = 1.0/(double)(N-1);
  double dt = dx/2.0;
  
  x=new double[N];
  a=new double[N];
  u=new double[N];
  w1=new double[N];
  w1_new=new double[N];
  w2=new double[N];
  w2_new=new double[N];
  

  // setting initial condition
  for(int i=0;i<N;i++){
    x[i]=(double)i*1.0/(double)N;
    w1[i] = 0.0;
    a[i] = 1.0;
    w1_new[i]=0.0;w2[i]=0.0;w2_new[i] = 0.0;
  }
  // setting initial condition for u
  for(int i=N/3+1;i<N/2;i++)
    u[i] = 3.0*(x[i]-1.0/3.0);
  for(int i=N/2;i<2*N/3+1;i++)
    u[i] = 3.0*(2.0/3.0-x[i]);
  
  for(int i=0;i<N;i++)
    out<<x[i]<<"\t"<<u[i]<<"\n";

  //setting initial conditions for w1 and w2
  for(int i=1;i<N-1;i++){
    w1[i] = sqrt(a[i]+1.)/2.0*(u[i+1]-u[i])/(dx); //sqrt(a[i]+1.)/2.0*(u[i+1]-u[i-1])/(2.0*dx);
    w2[i]= -sqrt(a[i]+1.)/2.0*(u[i+1]-u[i])/(dx);//-sqrt(a[i]+1.)/2.0*(u[i+1]-u[i-1])/(2.0*dx);
  }

  w1[0] = sqrt(a[0]+1.)/2.0*(u[1]-u[0])/dx;
  w1[N-1] = sqrt(a[N-1]+1.)/2.0*(u[N-1]-u[N-2])/dx;
  w2[0]=-sqrt(a[0]+1.)/2.0*(u[1]-u[0])/dx;
  w2[N-1]=-sqrt(a[N-1]+1.)/2.0*(u[N-1]-u[N-2])/dx;
  
  ofstream out_u;
  string filename;


  while(t<T){
    Lax_Wandroff_Step(w1,w2,w1_new,w2_new,dt,dx,a,N);
    out_energy<<t<<"\t"<<Energy_Compute(w1,w2,u,a,dt,dx,N)<<"\n";

    if((int)(1000*t)%10 == 0){
      filename ="data_" + IntToStr((int)(10*t)) + ".dat";
      out_u.open(filename.c_str());
      Integrator(w1,w2,u,a,N);
      for(int i=0;i<N;i++)
	out_u<<x[i]<<"\t"<<u[i]<<"\n";
      out_u.close();
    }
    
    t=t+dt;    
  }
  
  out.close();
  out_energy.close();
  delete[] x;
  delete[] a,u, w1, w2,w1_new, w2_new;

  return 0;

}
/****************************************************************************************/
void Lax_Wandroff_Step(double*w1, double* w2, double* w1_new, double* w2_new,double dt,double dx,double *a,int N){
  
  //scheme for the whole domain
  for(int i=1;i<N-1;i++){
    w1_new[i] = w1[i] + 0.5*dt/dx*sqrt(a[i])*(w1[i+1]-w1[i-1]) + 0.5*a[i]*dt*dt/dx/dx*(w1[i+1] -2.0*w1[i]+ w1[i-1]);
    w2_new[i] = w2[i] - 0.5*dt/dx*sqrt(a[i])*(w2[i+1]-w2[i-1]) + 0.5*a[i]*dt*dt/dx/dx*(w2[i+1] -2.0*w2[i]+ w2[i-1]);
  }
  // w2 is going right and w1 is going left
  // we solve w1 at first point and w2 at end point 
  // Why is it not working :)
  w1_new[0] = w1[0] + dt/dx*sqrt(a[0])*(w1[1]-w1[0]);// + 0.5*a[0]*dt*dt/dx/dx*(w1[2] -2.0*w1[1]+ w1[0]);
  w2_new[N-1] = w2[N-1] - dt/dx*sqrt(a[N-1])*(w2[N-1]-w2[N-2]);// + 0.5*a[N-1]*dt*dt/dx/dx*(w2[N-1] -2.0*w2[N-2]+ w2[N-3]);

  // because of first and end condition we have 
  // w2 = -w1 at x=0 or i=0
  // w1 = w2 at x=1  or i=N-1
  w2_new[0] = - w1[0];
  w1_new[N-1] = w2[N-1];
  
  for(int i=0;i<N;i++){
    w1[i] = w1_new[i];
    w2[i] = w2_new[i];
  }
  
}
/******************************************************************************************/
double Energy_Compute(double *w1,double *w2,double *u,double *a, double dt,double dx,int N){
  //computing energy using E=1/2 u_t^2 + 1/2*u^2 
  double energy = 0;
  Integrator(w1,w2,u,a,N);
  for(int i=3;i<N-3;i++){
    energy = energy + 0.5*a[i]/(a[i]+1.)*(w1[i]+w2[i])*(w1[i]+w2[i])*dx + 0.5*(w1[i]-w2[i])*(w1[i]-w2[i])*1./(a[i]+1.)*dx;//0.5*((u[i]) * (u[i]))*dx + 0.5*a[i]/(a[i]+1)*((w1[i]-w2[i])*(w1[i]-w2[i]))*dx ;
  }
  return energy;
}
/**************************************************************************************/
void Integrator(double *w1, double *w2, double * u,double *a, int N){
  
  double dx = 1.0/(double)N;
  u[0] = 0.0;
  for(int i=1;i<N;i++){
    u[i]=u[i-1] + 1./sqrt(a[i]+1)*dx*(w1[i-1]-w2[i-1]);
  }
  return;
}
