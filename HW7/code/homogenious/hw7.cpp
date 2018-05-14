#include<iostream>
#include<cmath>
#include<fstream>
#include<stdlib.h>
#include<ctime>
#include<string>
#include<sstream>

using namespace std;

const int Nx = 100;
const int Ny = 100;
const double Lx = 1.0;
const double Ly = 1.0;


int I(int i, int j){
  return i + j*Nx;
}

int main(){
  double * c;
  double* T;
  double TR,TU,TD,TL;
  double T1, T2;
  double dx1, dx2;
  double dx = Lx/(double) Nx;
  double dy = Ly/(double) Ny;
  double newT;
  double A,B,C;
  
  ofstream out;
  out.open("T.dat");

  double min;
  int index;

  T = new double[Nx*Ny];
  c = new double[Nx*Ny];

  for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++){
      T[I(i,j)] = 1000.0;
      c[I(i,j)] = 1.0;
    }
  }
  
  T[I(0,0)] = 0.0;
  
  
  for(int step =0; step <10;step++){
    
    
    for(int j=0;j<Ny;j++){
      for(int i=0;i<Nx;i++){
	if( i!=0 | j!=0){
	  
	  TL = ( i == 0 ? 100000 : T[I(i-1,j)]); //left
	  TD = ( j == 0 ? 100000 : T[I(i,j-1)]); //down
	  //TR = ( i == Nx-1 ? 100000 : T[I(i+1,j)]); //right
	  //TU = ( j == Ny-1 ? 100000 : T[I(i,j+1)]); //up
	  
	  /*
	    TL =  T[I(i-1,j)]; //left
	    TD = T[I(i,j-1)]; //down
	  TR =  T[I(i+1,j)]; //right
	  TU = T[I(i,j+1)]; //up
	  */
	  
	  if(TL<TD){
	    newT = TL + dx*c[I(i,j)];
	    if(newT <=TD)
	      T[I(i,j)] = newT;
	    else{
	      A = 2.0;
	      B = -2.0*(TL+TD);
	      C = TL*TL + TD*TD - c[I(i,j)]*c[I(i,j)]*dx*dx;
	      newT = (-B+sqrt(B*B-4.0*A*C) )/(2.0*A);
	      T[I(i,j)] = newT;
	    }
	  }else{
	    newT = TD + dx*c[I(i,j)];
	    if(newT <=TL)
	      T[I(i,j)] = newT;
	    else{
	      A = 2.0;
	      B = -2.0*(TL+TD);
	      C = TL*TL + TD*TD - c[I(i,j)]*c[I(i,j)]*dx*dx;
	      newT = (-B+sqrt(B*B-4.0*A*C) )/(2.0*A);
	      T[I(i,j)] = newT;
	    }
	  }
	}
      }
    }
  }




  
  
  
  
  
  
  
  for(int j=0;j<Ny;j++){
    //out<<"j"<<"\t";
    for(int i=0;i<Nx;i++){
      out<<"\t"<<T[I(i,j)]<<"\t";
    }
    out<<"\n";
  }
  
  return 0;

  }
