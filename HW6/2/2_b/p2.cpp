#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<sstream>

using namespace std;
int M=601;
double DX =  6.0;
double Xmin = -2.5;
double Ymin = -2.5;
double DY =  6.0;
double dx = DX/double(M-1);
double dy = DY/double(M-1);
double dt = 0.0001;
double A = 1.0;
double epsilon = 0.0;
int Steps = 7000;
int MM = 500;

int I(int i, int j){ return j*M+i;}
double Max ( double a, double b) { return a>b ? a : b;}
double Min(double a, double b){return a>b ? b : a; }
double distance(double X,double Y,double *curve_x,double *curve_y,int N);
double inside(double X,double Y,double * curve_x,double * curve_y,int N);
string IntToStr(int n){
  stringstream result;
  result<<n;
  return result.str();
}
//*************************************
int main(){
  cout<<"Initializing data..."<<endl;
  int N = 280;
  double *curve_x,*curve_y;
  double s,ds;
  double px,py,pxx,pyy,pxy,delp,delm;
  ofstream out;
  //temporary variables
  double temp, temp2;
  out.open("final.dat");
  ofstream out2;
  out2.open("initial.dat");

  ofstream out_file;
  string filename;
  
  
  
  int iNext,iPrev,jNext,jPrev;  

  ds=2.0*M_PI/(double)(N-1);
  curve_x=new double[N];
  curve_y=new double[N];

  double X,Y,*phi,*phin;
  phi = new double[(M)*(M)];
  phin = new double[(M)*(M)];

  cout<<"Writing the intial curve..."<<endl;
  
  for(int i=0;i<N;i++){
    s=0.0+ds*(double)i;
    curve_x[i] =1.0*(sin(s)*sin(s)*sin(s)) ;
    curve_y[i] = (13.0*cos(s)-5.0*cos(2.0*s)-2.0*cos(3.0*s)-cos(4.0*s))/16.0;
  }
  /*
  for(int i=0;i<N;i++){
    s = 0.0 + (double)i *2*M_PI/(double)N;
    curve_x[i] = cos(s);
    curve_y[i] = sin(s);
  }
  */
  
  
  for(int i=0;i<N+1;i++)
    out2<<curve_x[i]<<"\t"<<curve_y[i]<<"\n";

  
  /* This part is for testing the sign method
  cout<<"sign of  of 0.5, -1.5 is "<<inside(0.5,-1.5,curve_x,curve_y,N)<<endl;
  cout<<"sign of  of 1.5, -1 is "<<inside(1.5,-1.0,curve_x,curve_y,N)<<endl;
  cout<<"sign of  of 0.5, +1.0 is "<<inside(0.5,1.0,curve_x,curve_y,N)<<endl;
  cout<<"sign of  of 0.0, -2.0 is "<<inside(0.0,-2.0,curve_x,curve_y,N)<<endl;
  cout<<"sign of  of 0.75, -1.5 is "<<inside(0.75,-1.5,curve_x,curve_y,N)<<endl;
  */
  cout<<"Computing value of phi at each point ...\n";
  for(int i=0;i<M;i++){
    for(int j=0;j<M;j++){
      X = Xmin + (double)i/(double)M *DX;
      Y = Ymin + (double)j/(double)M *DY;
      phi[I(i,j)] = inside(X,Y,curve_x,curve_y,N)*distance(X,Y,curve_x,curve_y,N);
      phin[I(i,j)] = phi[(i,j)];
    }
  }


  filename ="data_" + IntToStr(0) + ".dat";
  out_file.open(filename.c_str());
  for(int j = 0;j<M;j++){
    for(int i=0;i<M;i++){
      X = Xmin + (double)i/(double)M *DX;
      Y = Ymin + (double)j/(double)M *DY;
      //inside(X,Y,curve_x,curve_y,N) * distance(X,Y,curve_x,curve_y,N);
      out_file<<X<<"\t"<<Y<<"\t"<<phi[I(i,j)]<<"\n";
    }
    out_file<<"\n";
  }
  
    out_file.close();

  cout<<"Starting the loop over time ....\n";
  //main loop for propagating the level set method
  for(int step = 0; step<Steps;step++){
    cout<<"Iteration at step :"<<step<<"\n";
    for(int i=0;i<M;i++){
      for(int j=0;j<M;j++){
	iNext = (i==M-1) ? 0 : i+1;
	iPrev = (i==0) ? M-1 : i-1;
	jNext = (j==M-1) ? 0 : j+1;
	jPrev = (j==0) ? M-1 : j-1;
	pxx = (phi[I(iNext,j)]-2.0*phi[I(i,j)]+phi[I(iPrev,j)])/(dx*dx);
	pyy = (phi[I(i,jNext)]-2.0*phi[I(i,j)]+phi[I(i,jPrev)])/(dy*dy);
	py = (phi[I(i,jNext)]-phi[I(i,jPrev)])/(2.0*dy);
	px = (phi[I(iNext,j)]-phi[I(iPrev,j)])/(2.0*dx);
	pxy = (phi[I(iNext,jNext)]+phi[I(iPrev,jPrev)]-phi[I(iPrev,jNext)]-phi[I(iNext,jPrev)])/(4.0*dy*dx);
	temp2 = pxx*py*py -2.0*px*py*pxy+ pyy*px*px;
	temp2 = temp2/(pow(pow(px,2)+ pow(py,2),1.5)+1e-6);

	temp2 = epsilon * temp2*sqrt(px*px+py*py);

	delp = sqrt(pow(fmax((phi[I(i,j)] - phi[I(iPrev,j)])/dx,0.0),2)
		    + pow(fmax((phi[I(i,j)] - phi[I(i,jPrev)])/dy,0.0),2)
		    + pow(fmin((phi[I(iNext,j)] - phi[I(i,j)])/dx,0.0),2)
		    + pow(fmin((phi[I(i,jNext)] - phi[I(i,j)])/dy,0.0),2));
	delm = sqrt(pow(fmin((phi[I(i,j)] - phi[I(iPrev,j)])/dx,0.0),2)
		    + pow(fmin((phi[I(i,j)] - phi[I(i,jPrev)])/dy,0.0),2)
		    + pow(fmax((phi[I(iNext,j)] - phi[I(i,j)])/dx,0.0),2)
		    + pow(fmax((phi[I(i,jNext)] - phi[I(i,j)])/dy,0.0),2));

	temp = fmax(A,0.0)*delp + fmin(A,0.0)*delm;
		     
	 /*
	temp = sqrt(pow(Max((phi[I(i,j)] - phi[I(iPrev,j)])/dx,0.0),2)
		    + pow(Max((phi[I(i,j)] - phi[I(i,jPrev)])/dy,0.0),2)
		    + pow(Min((phi[I(iNext,j)] - phi[I(i,j)])/dx,0.0),2)
		    + pow(Min((phi[I(i,jNext)] - phi[I(i,j)])/dy,0.0),2));
	 */
	phin[I(i,j)] = phi[I(i,j)] -dt*(temp - temp2);// + epsilon *dt *temp2;;
	
      }
    }
    /*    
    for(int i=0;i<M;i++){
      for(int j=0;j<M;j++){
	phi[I(i,j)] = phin[I(i,j)];
      }
    }
    */



    for(int i=1;i<M-1;i++){
      for(int j=1;j<M-1;j++){
	phi[I(i,j)] = phin[I(i,j)];
      }
    }

    
    for(int i=0;i<M;i++){
      phi[I(i,M-1)] = phi[I(i,M-2)];
      phi[I(i,0)] = phi[I(i,1)];
      phi[I(M-1,i)] = phi[I(M-2,i)];
      phi[I(0,i)] = phi[I(1,i)];
    }




    if(step%MM==0){
      filename ="data_" + IntToStr(step) + ".dat";
      out_file.open(filename.c_str());
      for(int j = 0;j<M;j++){
	for(int i=0;i<M;i++){
	  X = Xmin + (double)i/(double)M *DX;
	  Y = Ymin + (double)j/(double)M *DY;
	  //inside(X,Y,curve_x,curve_y,N) * distance(X,Y,curve_x,curve_y,N);
	  out_file<<X<<"\t"<<Y<<"\t"<<phi[I(i,j)]<<"\n";
	}
	out_file<<"\n";
      }
    }
    out_file.close();
    
  }

  for(int j = 0;j<M;j++){
    for(int i=0;i<M;i++){
      X = Xmin + (double)i/(double)M *DX;
      Y = Ymin + (double)j/(double)M *DY;
      //inside(X,Y,curve_x,curve_y,N) * distance(X,Y,curve_x,curve_y,N);
      out<<X<<"\t"<<Y<<"\t"<<phi[I(i,j)]<<"\n";
    }
    out<<"\n";
  }
    

  delete [] curve_x;
  delete [] curve_y;
  delete [] phi;
  return 0;
}


//************************************** Computing the distance from a point to the curve *******************
double distance(double X,double Y,double *curve_x,double *curve_y,int N){
  double min = 1e10;
  double dist=0.0;
  for(int i=0;i<N+1;i++){
    dist = sqrt(pow(X-curve_x[i],2)+pow(Y-curve_y[i],2));
    if(min>dist)
      min = dist;
  }
  return min;
}
//**************************************** Check if the point is inside or outside ********************
double inside(double X,double Y,double * curve_x,double * curve_y,int N){
  double w = 0;
  
  double r =0;
  for(int i=0;i<N-1;i++){
    if((curve_y[i]-Y)*(curve_y[i+1]-Y) < 0){
      r = (curve_x[i]-X) + (curve_y[i]-Y)*(curve_x[i+1] -curve_x[i]) / (curve_y[i]-curve_y[i+1]+1e-6);
      if(r>0){
	if(curve_y[i]-Y <0)
	  w=w+1;
	else
	  w = w-1;
      }
    }else if ( (curve_y[i] - Y == 0.0) & (curve_x[i] -X > 0)){
      if(curve_y[i+1] > 0)
	w = w+0.5;
      else 
	w = w- 0.5;
    }else if ((curve_y[i+1] - Y == 0.0) & (curve_x[i+1] - X >0.0)){
      if(curve_y[i] < 0.0)
	w = w +0.5;
      else
	w = w -0.5;
    }
  }
  
  if(w ==0.0)
    return +1.0;
  else 
    return -1.0;
  
  /*
  if((X>=0.0) & X<=1.0){
    if((Y>=-2) & Y<=cos(2*M_PI*X)){
      return -1.0;
    }
  }
  
  return +1.0;
  */
  /*
  if(X*X+ Y*Y <=1.0)
    return -1.0;
  else 
    return +1.0;
  */

}
//******************************************************************************************************
