#include<iostream>
#include<cmath>
#include<fstream>
#include<cstdlib>
#include<string>
#include<sstream>

using namespace std;
int M=256;
double DX =  5.0;
double Xmin = -2.5;
double Ymin = -2.5;
double DY =  5.0;
double dx = DX/double(M-1);
double dy = DY/double(M-1);
double dt = 0.001;
double dt_noise = 0.0001;
double epsilon = .0;
int Steps = 2000;
double probability=20;
bool noise_canceling = false;
int Steps_Noinse_cancelling = 200;
double A_noise = 1.0;
double epsilon_noise = 1.0;

int MM=200;
double Radius=0.4;



double distance(double X,double Y,double *curve_x,double *curve_y,int N);
double inside(double X,double Y,double * curve_x,double * curve_y,int N);
int I(int i, int j){ return j*M+i;}
string IntToStr(int n){
  stringstream result;
  result<<n;
  return result.str();
}
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
  int N = 311;
  double *curve_x,*curve_y;
  double s,ds;

  double px,py,pxx,pyy,pxy,delp,delm;
  ofstream out;
  //temporary variables
  double temp,temp2,F;
  out.open("shape_with_noise.dat");
  ofstream out2;
  out2.open("initialshape.dat");
  ofstream out3;
  out3.open("finalshape.dat");
  ofstream out_file;
  string filename;

  
  int iNext,iPrev,jNext,jPrev;  


  double X,Y,*phi,*phin,*intensity,*intensityn;
  phi = new double[M*M];
  phin = new double[M*M];
  intensity = new double [M*M];
  intensityn = new double [M*M];
  cout<<"initializing Phi ..."<<endl;
  ds=2.0*M_PI/(double)(N-1);
  curve_x=new double[N];
  curve_y=new double[N];
  
  for(int i=0;i<N;i++){
    s=0.0+ds*(double)i;
    curve_x[i] =Radius*cos(s);
    curve_y[i] = Radius*sin(s);
  }




  for(int i=0;i<M;i++){
    for(int j=0;j<M;j++){
      X = Xmin + (double)i/(double)M *DX;
      Y = Ymin + (double)j/(double)M *DY;
      phi[I(i,j)] = inside(X,Y,curve_x,curve_y,N)*distance(X,Y,curve_x,curve_y,N);
      phin[I(i,j)] = phi[(i,j)];
    }
  }



  filename ="data_initial.dat";
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
  
  


  /*
  for(int i=0;i<M;i++){
    for(int j=0;j<M;j++){
      X = Xmin + dx*(double)i;
      Y = Ymin + dy*(double)j;
      if(X*X+Y*Y<.25){
	phi[I(i,j)] = 0.0;
	phin[I(i,j)] = 0.0;
      }else{
	phi[I(i,j)] = 1.0;
	phin[I(i,j)] = 1.0;
      }
    }
  }
  */





  cout<<"Writing the shape of the image..."<<endl;
  // writing a circle in the middle ...
  for(int i=0;i<M;i++){
    for(int j=0;j<M;j++){
      X = Xmin + dx*(double)i;
      Y = Ymin + dy*(double)j;
      if(((X-DX/8.0)* (X-DX/8.0) + Y*Y)<=(DX*DX)/32.0 | ((X+DX/5.0)* (X+DX/7.0) + Y*Y)<=(DX*DX)/10. ){
      // if(((X-DX/8.0)* (X-DX/8.0) + Y*Y)<=(DX*DX)/32.0 | (abs(X)<0.750 & abs(Y)<0.750)){
	intensity[I(i,j)] = 0.0;
      }else{
	intensity[I(i,j)] = 255.0;
      }
    }
  }
  
  
  for(int i=0;i<M;i++){
    for(int j=0;j<M;j++){
      out2<<intensity[I(j,i)]<<"\t";
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

  for(int i=0;i<M;i++){
    for(int j=0;j<M;j++){
      out<<intensity[I(j,i)]<<"\t";
    }
    out<<"\n";
  }

  filename ="data_0.dat";
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



  //if noise cancelling is on do this to cancel the noise
  if(noise_canceling){
    ofstream out5;
    out5.open("initial_shape_after_noise_Canceling.dat");
    for(int step = 0; step<Steps_Noinse_cancelling;step++){
      cout<<"Iteration at step :"<<step<<"\n";
      for(int i=1;i<M-1;i++){
	for(int j=1;j<M-1;j++){
	  iNext = (i==M-1) ? 0 : i+1;
	  iPrev = (i==0) ? M-1 : i-1;
	  jNext = (j==M-1) ? 0 : j+1;
	  jPrev = (j==0) ? M-1 : j-1;
	  pxx = (intensity[I(iNext,j)]-2.0*intensity[I(i,j)]+intensity[I(iPrev,j)])/(dx*dx);
	  pyy = (intensity[I(i,jNext)]-2.0*intensity[I(i,j)]+intensity[I(i,jPrev)])/(dy*dy);
	  py  = (intensity[I(i,jNext)]-intensity[I(i,jPrev)])/(2.0*dy);
	  px  = (intensity[I(iNext,j)]-intensity[I(iPrev,j)])/(2.0*dx);
	  pxy = (intensity[I(iNext,jNext)]+intensity[I(iPrev,jPrev)]-intensity[I(iPrev,jNext)]-intensity[I(iNext,jPrev)])/(4.0*dy*dx);
	  temp2 = pxx*py*py -2.0*px*py*pxy+ pyy*px*px;
	  temp2 = temp2/(pow(px*px+ py*py,1.5)+1e-6);
	  
	  temp2 = epsilon_noise*temp2*sqrt(px*px+ py*py);
	  
	  delp = sqrt(pow(fmax((intensity[I(i,j)] - intensity[I(iPrev,j)])/dx,0.0),2)
		      + pow(fmax((intensity[I(i,j)] - intensity[I(i,jPrev)])/dy,0.0),2)
		      + pow(fmin((intensity[I(iNext,j)] - intensity[I(i,j)])/dx,0.0),2)
		      + pow(fmin((intensity[I(i,jNext)] - intensity[I(i,j)])/dy,0.0),2));
	  delm = sqrt(pow(fmin((intensity[I(i,j)] - intensity[I(iPrev,j)])/dx,0.0),2)
		      + pow(fmin((intensity[I(i,j)] - intensity[I(i,jPrev)])/dy,0.0),2)
		      + pow(fmax((intensity[I(iNext,j)] - intensity[I(i,j)])/dx,0.0),2)
		      + pow(fmax((intensity[I(i,jNext)] - intensity[I(i,j)])/dy,0.0),2));
	  
	  temp = fmax(A_noise,0)*delp + fmin(A_noise,0)*delm;
	  intensityn[I(i,j)] = intensity[I(i,j)]+dt_noise*(temp + temp2);
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
   
    redo_intensity(intensity);
    
    
    for(int i=0;i<M;i++){
      for(int j=0;j<M;j++){
	out5<<intensity[I(j,i)]<<"\t";
      }
      out5<<"\n";
    }    
    
  }    
  
  
  cout<<"Starting the loop over time ....\n";
  //main loop for propagating the level set method
  for(int step = 1; step<Steps;step++){
    cout<<"Iteration at step :"<<step<<"\n";
    for(int i=1;i<M-1;i++){
      for(int j=1;j<M-1;j++){
	iNext = (i==M-1) ? 0 : i+1;
	iPrev = (i==0) ? M-1 : i-1;
	jNext = (j==M-1) ? 0 : j+1;
	jPrev = (j==0) ? M-1 : j-1;

	pxx = (phi[I(iNext,j)]-2.0*phi[I(i,j)]+phi[I(iPrev,j)])/(dx*dx);
	pyy = (phi[I(i,jNext)]-2.0*phi[I(i,j)]+phi[I(i,jPrev)])/(dy*dy);
	py  = (phi[I(i,jNext)]-phi[I(i,jPrev)])/(2.0*dy);
	px  = (phi[I(iNext,j)]-phi[I(iPrev,j)])/(2.0*dx);
	pxy = (phi[I(iNext,jNext)]+phi[I(iPrev,jPrev)]-phi[I(iPrev,jNext)]-phi[I(iNext,jPrev)])/(4.0*dy*dx);
	temp2 = pxx*py*py -2.0*px*py*pxy+ pyy*px*px;
	temp2 = temp2/(pow(px*px+ py*py,1.5)+1e-6);
	
	temp2 = epsilon*temp2*sqrt(px*px+ py*py);

	py  = (intensity[I(i,jNext)]-intensity[I(i,jPrev)])/(2.0*dy);
	px  = (intensity[I(iNext,j)]-intensity[I(iPrev,j)])/(2.0*dx);
	F = -1.0/(1.0 + sqrt(px*px+ py*py));


	delp = sqrt(pow(fmax((phi[I(i,j)] - phi[I(iPrev,j)])/dx,0.0),2)
		    + pow(fmax((phi[I(i,j)] - phi[I(i,jPrev)])/dy,0.0),2)
		    + pow(fmin((phi[I(iNext,j)] - phi[I(i,j)])/dx,0.0),2)
		    + pow(fmin((phi[I(i,jNext)] - phi[I(i,j)])/dy,0.0),2));
	delm = sqrt(pow(fmin((phi[I(i,j)] - phi[I(iPrev,j)])/dx,0.0),2)
		    + pow(fmin((phi[I(i,j)] - phi[I(i,jPrev)])/dy,0.0),2)
		    + pow(fmax((phi[I(iNext,j)] - phi[I(i,j)])/dx,0.0),2)
		    + pow(fmax((phi[I(i,jNext)] - phi[I(i,j)])/dy,0.0),2));
	
	temp = fmin(F,0)*delp + fmax(F,0)*delm;
	phin[I(i,j)] = phi[I(i,j)]+dt*(temp+temp2);	
      }
    }
    for(int i=0;i<M;i++){
      phin[I(0,i)] = phin[I(1,i)];
      phin[I(M-1,i)] = phin[I(M-2,i)];
      phin[I(i,0)] = phin[I(i,1)];
      phin[I(i,M-1)] = phin[I(i,M-2)];
    }
    
    for(int i=0;i<M;i++){
      for(int j=0;j<M;j++){
	phi[I(i,j)] = phin[I(i,j)];
      }
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
      out_file.close();      
    }
    



  }


  for(int i=0;i<M;i++){
    for(int j=0;j<M;j++){
      out3<<phi[I(i,j)]<<"\t";
    }
    out3<<"\n";
  }    
  delete [] phi;
  delete [] phin;
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
  */
  
  if(X*X+ Y*Y <=Radius){
      return -1.0;
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
