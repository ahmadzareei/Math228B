#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdio>
#include<sstream>
#include<string>

using namespace std;


void step_diffusion_solver(int N, double*** u,double*** u_new, double landa,int T);
double mean_value(int N, double*** u);

string IntToStr(int n){
  stringstream result;
  result<<n;
  return result.str();
}



int main(){
  ofstream myfile;
  string filename;
  
  for(int kk=1;kk<5;kk++){
    
    filename="p3_average_u_" + IntToStr(kk) +".dat";
    myfile.open(filename.c_str());

    double fraction = 0.1*kk;

    double k=1./600.0;
    double h=1./10.0;
    
    int N;
    N=(int)(1.0/h);
    N=N+1;

    
    //setting array for  u 
    double*** u; 
    u=new double** [N];
    for(int i=0;i<N;i++){
      u[i]=new double* [N];
      for(int l=0;l<N;l++)
	u[i][l]=new double [N];
    }

    // initial data 
    for(int r=0;r<N;r++)
      for(int s=0;s<N;s++)
	for(int t=0;t<N;t++)
	  u[r][s][t] = 0.0;

    //top boundary conditon
    for(int r=0;r<N;r++)
      for(int s=0;s<N;s++)
	u[r][s][N-1] = 1.0;
    

    double landa = k/(h*h);
    double t=0.0;
    double average;
    double*** u_new;
    u_new=new double** [N];
    for(int i=0;i<N;i++){
      u_new[i]=new double* [N];
      for(int l=0;l<N;l++)
	u_new[i][l]=new double [N];
    }
    for(int r=0;r<N;r++)
      for(int s=0;s<N;s++)
	for(int t=0;t<N;t++)
	  u_new[r][s][t] = u[r][t][s];

    double tt;
    int TT=0;
    while(t<5.0){
      tt = t - floor(t);
      if(tt<=fraction)
	TT=1;
      else
	TT=0;
      step_diffusion_solver(N,u,u_new,landa,TT);
      average= mean_value(N,u);
      myfile<<t<<"\t"<<average<<"\n";
      t=t+k;
    }

    myfile.close();

    delete[] u;
    delete[] u_new;
    
  }

}

//***********************************************************************************
void step_diffusion_solver(int N, double*** u,double*** u_new, double landa, int T){

  for(int i=1;i<N-1;i++)
    for(int l=1;l<N-1;l++)
      for(int m=1;m<N-1;m++)
	u_new[i][l][m]=(1.0-6.0*landa)*u[i][l][m] + landa*(u[i+1][l][m] + u[i-1][l][m] + u[i][l+1][m] +u[i][l-1][m] + u[i][l][m+1] + u[i][l][m-1]);
  
  for(int r=1;r<N-1;r++)
    for(int s=1;s<N-1;s++){
      // BC on x=0 and x=1
      u_new[0][r][s] = u_new[2][r][s];
      u_new[N-1][r][s] = u_new[N-3][r][s];
      //BC on y=0 and y=1
      u_new[r][0][s] = u_new[r][2][s];
      u_new[r][N-1][s] = u_new[r][N-3][s];
      //BC on z=0
      u_new[r][s][0] = u_new[r][s][2];
      }
  if(T==1)
    for(int r=1;r<N-1;r++)
      for(int s=1;s<N-1;s++)
	u_new[r][s][N-1] = 1.0;  
  else 
    for(int r=1;r<N-1;r++)
      for(int s=1;s<N-1;s++)
	u_new[r][s][N-1] = 0.0;  
    
  for(int i=0;i<N;i++)
    for(int l=0;l<N;l++)
      for(int m=0;m<N;m++)
	u[i][l][m]=u_new[i][l][m];
  
}
//***********************************************************************************
double  mean_value(int N, double*** u){
  double ave = 0.0;
  for(int i=1;i<N-1;i++)
    for(int l=1;l<N-1;l++)
      for(int m=1;m<N-1;m++)
	ave =ave + u[i][l][m];
  ave = ave/((double)N-2.0)/((double)N-2.0)/((double)N-2.0);
  return ave;
}
