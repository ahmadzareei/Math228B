#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdio>
#include<sstream>
#include<string>

using namespace std;

void diffusion_solver(int N, double h, double k, double* u);

string IntToStr(int n){
  stringstream result;
  result<<n;
  return result.str();
}



int main(){
  ofstream myfile;
  string filename;
  
  for(int kk=1;kk<5;kk++){
    filename="p2_d_t_02_h_" + IntToStr(kk) +".dat";
    myfile.open(filename.c_str());

    
    double k=1./600.0/(double)kk/(double)kk ;
    double h=1./10.0/(double)kk;
    
    int N;
    N=(int)(1.0/h);
    N=N+1;
    double* u; 
    u=new double [N];

    for(int i=0;i<N;i++)
      u[i]=1.0;
    u[N-1]=0.0;
    
  //  for(int i=1;i<N;i++)
  //  initial<<(i-1)*h<<"\t"<<u[i]<<"\n";
  //initial.close();
    
    
    diffusion_solver(N,h,k,u);
    
    for(int i=1;i<N;i++)
      myfile<<(i)*h<<"\t"<<u[i]<<"\n";
    myfile.close();

    delete[] u;
    
   }

}


void diffusion_solver(int N, double h, double k, double* u){
  double landa = k/(h*h);
  double t=0.0;
  double* u_new;
  u_new = new double[N];

  while(t<=0.2){
    for(int i=1;i<N-1;i++)
      u_new[i]=landa*(1.0-1.0/((double)i))*u[i-1] + (1.0-2.0*landa)*u[i] + landa*(1.0+1.0/((double)i))*u[i+1];
    u_new[0] = 1;
    u_new[N-1]=0.0;
    for(int i=0;i<N;i++)
      u[i]=u_new[i];
    t=t+k;
  }
  delete[] u_new;
}
