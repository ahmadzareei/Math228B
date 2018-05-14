#include<iostream>
#include<cmath>
#include<fstream>
#include<stdlib.h>
#include<ctime>
#include<string>
#include<sstream>
#include<vector>

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
  double dx = 0.01;//1.0;
  double dy = 0.01;//1.0;
  double newT;
  double A,B,C;
  int indexi, indexj;
  double der_x, der_y;
  double xx,yy;

  vector<double> pathx;
  vector<double> pathy;
  

  pathx.push_back(75.0);
  pathy.push_back(75.0);


  ofstream out;
  out.open("T.dat");
  
  ofstream out2;
  out2.open("path.dat");


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

  /*
  for(int i=Nx/2;i<Nx;i++){
    for(int j=0;j<Ny;j++){
      c[I(i,j)] = 2.0;
    }
  }
  */
  

  T[I(0,0)] = 0.0;
  /*
  for(int i= 20;i<30;i++){
    for(int j=5;j<25;j++){
	c[I(i,j)] = 3.0;
    }
  }
  */
  
  for(int j=40;j<60;j++){
    for(int i=45;i<55;i++){
      c[I(i,j)] = 3.0;//000;
      //      T[I(i,j)] = 1000000.;
    }
  }
  
  

  
  for(int step =0; step <10;step++){
    
    
    for(int j=0;j<Ny;j++){
      for(int i=0;i<Nx;i++){
	if( i!=0 | j!=0){
	  
	  TL = ( i == 0 ? 100000 : T[I(i-1,j)]); //left
	  TD = ( j == 0 ? 100000 : T[I(i,j-1)]); //down
	  //TR = ( i == Nx-1 ? 100000 : T[I(i+1,j)]); //right
	  //TU = ( j == Ny-1 ? 100000 : T[I(i,j+1)]); //up
	  
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
	/*
	for(int i= 20;i<30;i++){
	  for(int j=5;j<25;j++){
	    T[I(i,j)] = 1000000.0;
	  }
	}
	*/
	
	for(int j=40;j<60;j++){
	  for(int i=45;i<55;i++){
	    //	    T[I(i,j)] = 100.0;
	  }
	}
	
      }
      
    }
  }
  double temp;
  indexi = 75;
  indexj = 75;
  
  while(indexi !=0 & indexj !=0){
    
    indexi = round(pathx.at(pathx.size()-1));
    indexj = round(pathy.at(pathy.size()-1));
    
    der_x = (T[I(indexi+1,indexj)]-T[I(indexi,indexj)])/dx;
    
    der_y = (T[I(indexi,indexj+1)]-T[I(indexi,indexj)])/dx;
    //temp  = der_x;
        
    //    cout<<der_x<<endl;
    //cout<<der_y<<endl;

    
    
    xx = pathx.at(pathx.size()-1) - der_x * dx;
    yy = pathy.at(pathx.size()-1) - der_y * dy;
  
    
    
    pathx.push_back(xx);
    pathy.push_back(yy);
    //cout<<der_x<<"\t"<<der_y<<endl;
    cout<<pathx.at(pathx.size()-1)<<"\t"<<pathy.at(pathx.size()-1)<<endl;
  }
  

  for(int i=0;i<pathx.size();i++){
    out2<<pathx.at(i)<<"\t"<<pathy.at(i)<<endl;
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
