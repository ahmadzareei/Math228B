#include<iostream>
#include<cmath>
#include<fstream>
#include<cstdlib>


using namespace std;
int M=256;
double DX =  5.0;
double Xmin = -2.5;
double Ymin = -2.5;
double DY =  5.0;
double dx = DX/double(M-1);
double dy = DY/double(M-1);
double dt = 0.0001;

int Steps = 5;
double probability=1;

int I(int i, int j){ return j*M+i;}
void redo_intensity(double* intensity){
  for(int i=0;i<M;i++){
    for(int j=0;j<M;j++){
      if(intensity[I(i,j)]<127.0)
	intensity[I(i,j)] = 0.0;
      else 
	intensity[I(i,j)] = 255.0;
    }
  }
}
//*************************************
int main(){
  cout<<"Initializing data..."<<endl;

  double pxx,pyy;
  ofstream out;
  //temporary variables
  double temp, temp2;
  out.open("image_with_noise.dat");
  ofstream out2;
  out2.open("initialshape.dat");
  ofstream out3;
  out3.open("finalshape_after_renorm.dat");
  ofstream out4;
  out4.open("finalshape_before_renorm.dat");

  
  int iNext,iPrev,jNext,jPrev;  


  double X,Y,*intensity,*intensityn;
  intensity = new double[M*M];
  intensityn = new double[M*M];

  cout<<"Writing the shape of the image..."<<endl;
  
  // writing a circle in the middle ...
  for(int i=0;i<M;i++){
    for(int j=0;j<M;j++){
      temp = Xmin + dx*(double)i;
      temp2 = Ymin + dy*(double)j;
      if(((temp-DX/5.0)* (temp-DX/5.0) + temp2*temp2)<=(DX*DX)/16.0 | ((temp+DX/5.0)* (temp+DX/5.0) + temp2*temp2)<=(DX*DX)/16.0 ){
	intensity[I(i,j)] = 0.0;
	intensityn[I(i,j)] = 0.0;}
      else{
	intensity[I(i,j)] = 255.0;
	intensityn[I(i,j)] = 255.0;
      }
    }
  }
  




  
  for(int i=0;i<M;i++){
    for(int j=0;j<M;j++){
      out2<<intensity[I(i,j)]<<"\t";
    }
    out2<<"\n";
  }

  cout<<"Making noise on the initial data ...";
  for(int i=1;i<M-1;i++){
    for(int j=1;j<M-1;j++){
      temp = rand()%101;
      if(temp<probability)
	intensity[I(i,j)] = rand()%256;
    }
  }
    
  redo_intensity(intensity);
  redo_intensity(intensityn);

  for(int i=0;i<M;i++){
    for(int j=0;j<M;j++){
      out<<intensity[I(i,j)]<<"\t";
    }
    out<<"\n";
  }

  cout<<"Starting the loop over time ....\n";
  //main loop for propagating the level set method
  for(int step = 0; step<Steps;step++){
    cout<<"Iteration at step :"<<step<<"\n";
    for(int i=1;i<M-1;i++){
      for(int j=1;j<M-1;j++){
	iNext = (i==M-1) ? 0 : i+1;
	iPrev = (i==0) ? M-1 : i-1;
	jNext = (j==M-1) ? 0 : j+1;
	jPrev = (j==0) ? M-1 : j-1;


	pxx = (intensity[I(iNext,j)]-2.0*intensity[I(i,j)]+intensity[I(iPrev,j)])/(dx*dx);
	pyy = (intensity[I(i,jNext)]-2.0*intensity[I(i,j)]+intensity[I(i,jPrev)])/(dy*dy);

	//solving heat equation
	intensityn[I(i,j)] = intensity[I(i,j)]+dt*(pxx+pyy);
      
      }
    }
    for(int i=1;i<M-1;i++){
      intensityn[I(0,i)] = intensityn[I(1,i)];
      intensityn[I(M-1,i)] = intensityn[I(M-2,i)];
      intensityn[I(i,0)] = intensityn[I(i,1)];
      intensityn[I(i,M-1)] = intensityn[I(i,M-2)];
    }

    for(int i=1;i<M-1;i++){
      for(int j=1;j<M-1;j++){
	intensity[I(i,j)] = intensityn[I(i,j)];
      }
    }
  }
  
  for(int i=0;i<M;i++){
    for(int j=0;j<M;j++){
      out4<<intensity[I(i,j)]<<"\t";
    }
    out4<<"\n";
  }    


    redo_intensity(intensity);
      

  for(int i=0;i<M;i++){
    for(int j=0;j<M;j++){
      out3<<intensity[I(i,j)]<<"\t";
    }
    out3<<"\n";
  }    
  delete [] intensity;
  delete [] intensityn;
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
  /*
  double r =0;
  for(int i=0;i<N+1;i++){
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
  */

  if((X>=0.0) & X<=1.0){
    if((Y>=-2) & Y<=cos(2*M_PI*X)){
      return -1.0;
    }
  }
  
  return +1.0;
  /*
  if(X*X+ Y*Y <=1.0)
    return -1.0;
  else 
    return +1.0;
  */

}
//******************************************************************************************************
