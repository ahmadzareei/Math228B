#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>

void heat_solver(int N, double h, double k, double* u);
double norm(int,double,double*);

using namespace std;

double pi=acos(-1.0);

int main(){
  int t_start, t_stop;
  //ofstream myfile;
  //myfile.open("p1.dat");
  cout.precision(15);
  double x;
  double k=1./200.0 ;
  double h=1./10.0;
  double error;
  int N;
  N=(int)(1.0/h);
  N=N+1;
  double* u; 
  u=new double [N];
  t_start = clock();
  for(int jj=0;jj<100;jj++){
    for(int i=0;i<N;i++){
      x= h*(double)(i);
      u[i]=x*(1-x);
    }
    u[0]=0.0;
    u[N-1]=0.0;  
    heat_solver(N,h,k,u);
  }
  t_stop = clock();
  error = norm(N,h,u);
  cout<<"Error: "<<error<<"\time: "<<(double)(t_stop-t_start)/(CLOCKS_PER_SEC)/100.0<<"\n";

  t_start = clock();
  for(int jj=0;jj<100;jj++){
    for(int i=0;i<N;i++){
      x=h*double(i);
      u[i] = (8.0/pi/pi)*(exp(-pi*pi)*sin(pi*x) + exp(-9.0*pi*pi)/9.0*sin(3*pi*x));
    }
  }
  
  t_stop = clock();
  cout<<"Time:"<<(double)(t_stop-t_start)/(CLOCKS_PER_SEC)/1000.<<"\n";
  //  for(int i=0;i<N;i++)
  //  myfile<<i*h<<"\t"<<u[i]<<"\n";
    
  //  myfile.close();
}

double norm(int N, double h, double *u){
  double diff;
  double x;
  double max =0.0;
    for(int i=0;i<N;i++){
      x=h*(double)(i);
      diff= abs(u[i]-(8.0/pi/pi/pi)*(exp(-pi*pi)*sin(pi*x) + exp(-9.0*pi*pi)/9.0*sin(3*pi*x)+exp(-25*pi*pi)/25.0*sin(5*pi*x)));
      if(diff>max)
	max = diff;
    }
  return max;
}

void heat_solver(int N, double h, double k, double* u){
  double landa = k/(h*h);
  double t=0.0;
  double* u_new;
  u_new = new double[N];
  while(t<=1.0){
    for(int i=1;i<N-1;i++)
      u_new[i]=landa*u[i-1] + (1.0-2.0*landa)*u[i] + landa*u[i+1];
    for(int i=1;i<N-1;i++)
      u[i]=u_new[i];
    t=t+k;
  }
}
