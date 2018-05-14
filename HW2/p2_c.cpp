#include<iostream>
#include<fstream>
#include<cmath>

void diffusion_solver(int N, double h, double k, double* u);
using namespace std;


int main(){
  const double pi=acos(-1.0);
  ofstream myfile;
  myfile.open("p2_c_p_07_2.dat");
  
  ofstream initial;
  initial.open("p2_c_initial.dat");

  double k=1./600.0/16.0/2.0 ;
  double h=1./10.0/4.0;

  int N;
  N=(int)(1.0/h);
  N=N+1;
  N=N+1;
  double* u; 
  u=new double [N];
  for(int i=0;i<N;i++)
    u[i]=1.0;
  u[N-1]=0.0;
  
  for(int i=1;i<N;i++)
    initial<<(i-1)*h<<"\t"<<u[i]<<"\n";
  initial.close();


  diffusion_solver(N,h,k,u);

  for(int i=1;i<N;i++)
    myfile<<(i-1)*h<<"\t"<<u[i]<<"\n";
  myfile.close();


}


void diffusion_solver(int N, double h, double k, double* u){
  double landa = k/(h*h);
  double t=0.0;
  double* u_new;
  u_new = new double[N];
  while(t<=0.7){
    for(int i=2;i<N-1;i++)
      u_new[i]=landa*(1.0-1.0/((double)i-1.0))*u[i-1] + (1.0-2.0*landa)*u[i] + landa*(1.0+1.0/((double)i-1.0))*u[i+1];
    u_new[0] = u_new[2];
    u_new[1] = 3.0*landa*u[0] + (1.0-2.0*3.0*landa)*u[1] + 3.0*landa*u[2];
    u_new[N-1]=0.0;
    for(int i=0;i<N;i++)
      u[i]=u_new[i];
    t=t+k;
  }
}
