#include<iostream>
#include<fstream>
#include<cmath>
void heat_solver(int N, double h, double k, double* u);
double norm(int N,double* u,double h, double* u_exact,int N_exact);
using namespace std;
double pi=acos(-1.0);


int main(){
  ofstream myfile;
  myfile.open("p133.dat");

  //solving equation for finest one
  double h=1./512.0;
  double k=1.0/6.0*h*h;

  int N;
  N=(int)(1.0/h);
  N=N+1;
  double* u_exact;
  double* u;
  int N_exact = N;
  u=new double[N];
  u_exact=new double [N];
  for(int i=0;i<N;i++)
    u_exact[i]=sin(pi*h*(double)(i));
  u_exact[0]=0.0;
  u_exact[N-1]=0.0;
  heat_solver(N,h,k,u_exact);
  
  double landa[2];
  landa[0]=1.0/2.0;
  landa[1]=1.0/6.0;
  int MM=5;
  double error[2][MM];
  for(int i=0;i<2;i++){
    h=1.0/8.0;
    for(int j=0;j<MM;j=j++){
      h=h/2.0;
      k=landa[i]*h*h;
      N=(int)(1.0/h);
      N=N+1;
      for(int kk=0;kk<N;kk++)
	u[kk]=sin(pi*h*(double)(kk));
      u[0]=0.0;
      u[N-1]=0.0;
      heat_solver(N,h,k,u);
      error[i][j]=norm(N,u,h,u_exact,N_exact);
    }
  }
  for(int i=0;i<MM;i++)
    myfile<<1.0/8.0/pow(2.0,i)<<"\t"<<error[0][i]<<"\t"<<error[1][i]<<"\n";
  myfile.close();
}


double norm(int N,double* u,double h, double* u_exact,int N_exact){
  double diff;
  double kk = (double)(N_exact-1)/(N-1);
  double max =0.0;
  //  max = abs (u[N/2]- u_exact[N_exact/2]);
  //max = abs(u[N/2] - exp(-pi*pi)*sin(pi/2));
    for(int i=0;i<N;i++){
      diff= abs(u[i]-exp(-pi*pi)*sin(pi*i*h));
      //diff= abs(u[i]-u_exact[(int)(i*kk)]);
      if(diff>max)
	max = diff;
    }
  return max;
}


void heat_solver(int N, double h, double k, double* u){
  double la = k/(h*h);
  double t=0.0;
  double* u_new;
  u_new = new double[N];
  while(t<=1.0){
    for(int i=1;i<N-1;i++)
      u_new[i]=la*u[i-1] + (1.0-2.0*la)*u[i] + la*u[i+1];
    for(int i=1;i<N-1;i++)
      u[i]=u_new[i];
    t=t+k;
  }
}
