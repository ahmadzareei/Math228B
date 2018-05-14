#include<iostream>
#include<fstream>
#include<cmath>
void heat_solver(int N, double h, double k, double* u);
using namespace std;

int main(){
  const double pi=acos(-1.0);
  //  ofstream myfile;
  //myfile.open("p1.dat");

  double k=1./600.0 ;
  double h=1./10.0;

  int N;
  N=(int)(1.0/h);
  N=N+1;
  double* u; 
  u=new double [N];
  for(int i=0;i<N;i++)
    u[i]=sin(pi*h*(double)(i));
  u[0]=0.0;
  u[N-1]=0.0;
  
  heat_solver(N,h,k,u);
  for(int i=0;i<N;i++)
    cout<<i*h<<"\t"<<u[i]<<"\n";
  // myfile.close();
}


void heat_solver(int N, double h, double k, double* u){
  double landa = k/(h*h);
  double t=0.0;
  double* u_new;
  u_new = new double[N];
  while(t<=1.0){
    for(int i=1;i<(N-1)/2;i++)
      u_new[i]=landa*u[i-1] + (1.0-2.0*landa)*u[i] + landa*u[i+1];
    for(int i=(N-1)/2;i<(N-1);i++)
      u_new[i]=landa/2.*u[i-1] + (1.0-2.0*landa/2.)*u[i] + landa/2.*u[i+1];    
    for(int i=1;i<N-1;i++)
      u[i]=u_new[i];
    t=t+k;
  }
}
