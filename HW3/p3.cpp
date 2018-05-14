#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdio>
#include<sstream>
#include<string>

using namespace std;

void step_diffusion_solver(int ,int ,int,int,int,double *** ,double *** ,double ***, double, double *** _p, double *** );
bool check_smell_room3(int ,int,int,double ***);
void  rhs_solve (int ,int,int ,double *** , double *** , double );
void  boundary_wall(int Nx_room,int Ny_room,int Nz,double *** u);
string IntToStr(int n){
  stringstream result;
  result<<n;
  return result.str();
}

int main(){
  ofstream myfile,out2,out3,out4;
  string filename;
 
  //  for(int kk=1;kk<4;kk++){

   int kk=1;

  

  //opening file for writing the values
  filename="p1_vlaue_" + IntToStr(kk) +".dat";
  
  myfile.open(filename.c_str());
  out2.open("hall_data.dat");
  out3.open("room_data.dat");
  out4.open("secondroom.dat");
  //increments in x,y,z in h, increments in t is k
  //  double k=1./100.0;
  //double h=1./4.0;
 
  double k=0.2*0.2/6.0;//1./150.0/(double)kk/(double)kk;
  double h=0.2; //1./5.0/(double)kk;
 

 int Nx_room,Ny_room,Nz,Nx_hall,Ny_hall;
 
  int Nodes_per_meter;
  Nodes_per_meter=(int)(1.0/h);
  //  Nodes_per_meter=Nodes_per_meter+1;
  Nx_room = 10*Nodes_per_meter;
  Ny_room = 7*Nodes_per_meter;

  Nx_hall = 26*Nodes_per_meter;
  Ny_hall = 3*Nodes_per_meter;

  Nz = 3*Nodes_per_meter;

  Nx_room+=1;
  Ny_room+=1;
  Nz+=1;
  Nx_hall+=1;
  Ny_hall+=1;
  
  //setting array for  u 
  double *** room9,*** room3,***hall; 
  double *** room_p,***hall_p; 

  room9=new double** [Nx_room];
  room3=new double** [Nx_room];
  for(int i=0;i<Nx_room;i++){
    room9[i]=new double* [Ny_room];
    room3[i]=new double* [Ny_room];
    for(int l=0;l<Ny_room;l++){
      room9[i][l]=new double [Nz];
      room3[i][l]=new double [Nz];
    }
  }
  
  hall = new double** [Nx_hall];
  for(int i=0;i<Nx_hall;i++){
    hall[i] = new double* [Ny_hall];
    for(int l=0;l<Ny_hall;l++)
      hall[i][l] = new double [Nz];
  }
  
  room_p=new double** [Nx_room];
  for(int i=0;i<Nx_room;i++){
    room_p[i]=new double* [Ny_room];
    for(int l=0;l<Ny_room;l++){
      room_p[i][l]=new double [Nz];
    }
  }

  
  hall_p = new double** [Nx_hall];
  for(int i=0;i<Nx_hall;i++){
    hall_p[i] = new double* [Ny_hall];
    for(int l=0;l<Ny_hall;l++)
      hall_p[i][l] = new double [Nz];
  }

  // initial data 
  for(int r=0;r<Nx_room;r++)
    for(int s=0;s<Ny_room;s++)
      for(int t=0;t<Nz;t++){
	room3[r][s][t] = 0.0;
	room9[r][s][t] = 0.0;
	room_p[r][s][t] = 0.0;
      }

  for(int r=0;r<Nx_hall;r++)
    for(int s=0;s<Ny_hall;s++)
      for(int t=0;t<Nz;t++){
	hall[r][s][t] = 0.0;
	hall_p[r][s][t] = 0.0;
      }

  //top boundary conditon
  //for(int r=(int)(0.25*Nz);r<(int)(0.75*Nz);r++)
  //for(int  t=0;t<Nz;t++)
  //  for(int iii=-1;iii<=1;iii++)
  //  for(int jjj=-1;jjj<=+1;jjj++)
  //    for(int kkk=-1;kkk<=1;kkk++)
	room9[Nx_room/2][Ny_room/2][Nz/2] = 100.0;
  //  for(int r=0;r<Nx_room;r++)
  //  for(int s=0;s<Ny_room;s++)
  //    for(int  t=0;t<Nz;t++)
  //	room9[r][s][t] = 100.0;

  double landa = k/(h*h);
  double t=0.0;
  double average;
  
  bool final_state = false;


  while(!final_state){
    step_diffusion_solver(Nx_room,Ny_room,Nz,Nx_hall,Ny_hall,room3,room9,hall,landa,room_p,hall_p);
    final_state= check_smell_room3(Nx_room, Ny_room,Nz,room3); //room3
    t=t+k;
    myfile<<t<<"\t"<<room3[Nx_room/2][Ny_room/2][Nz/2]<<"\n";
  }

  cout<<t<<"\n";
  
  for(int i=0;i<Nx_hall;i++)
    out2<<(double)i *26.0/(double)Nx_hall<<"\t"<<hall[i][Ny_hall/2][Nz/2]<<"\n";

  for(int i=0;i<Nx_room;i++)
    out3<<(double)i *10.0/(double)Nx_room<<"\t"<<room3[i][Ny_room/2][Nz/2]<<"\n";

  for(int i=0;i<Nx_room;i++)
    out4<<(double)i *10.0/(double)Nx_room<<"\t"<<room9[i][Ny_room/2][Nz/2]<<"\n";

  out2.close();
  out4.close();
  out3.close();
  
  //closing file 
  myfile.close();
  
  //free space
  delete [] room3,room_p, room9, hall, hall_p;  
  //}

}

//***********************************************************************************
void step_diffusion_solver(int Nx_room,int Ny_room, int Nz,int Nx_hall,int Ny_hall,double *** room3,double *** room9, double *** hall,double landa, double *** room_p, double *** hall_p){



  //we first solve room number 9
  rhs_solve(Nx_room,Ny_room,Nz,room_p,room9,landa);
  boundary_wall(Nx_room,Ny_room,Nz,room_p);
    
  //implementing Boundary of doors for room 9
  for(int t=0;t<Nz;t++){
    for(int r=(int)(0.2*Nx_room);r<(int)(0.3*Nx_room);r++)
      room_p[r][Ny_room-1][t]= (1.0-6.0*landa)*room9[r][Ny_room-1][t] + landa*(room9[r+1][Ny_room-1][t] + room9[r-1][Ny_room-1][t] + hall[r][1][t] + room9[r][Ny_room-2][t] + room9[r][Ny_room-1][t+1]+room9[r][Ny_room-1][t-1]);
  }
  for(int t=0;t<Nz;t++){
    for(int r=(int)(0.6*Nx_room);r<(int)(0.7*Nx_room);r++)
      room_p[r][Ny_room-1][t]= (1.0-6.0*landa)*room9[r][Ny_room-1][t] + landa*(room9[r+1][Ny_room-1][t] + room9[r-1][Ny_room-1][t] + hall[r][1][t] + room9[r][Ny_room-2][t] + room9[r][Ny_room-1][t+1]+room9[r][Ny_room-1][t-1]);
  }
  
  for(int r=0;r<Nx_room;r++)
    for(int s=0;s<Ny_room;s++)
      for(int t=0;t<Nz;t++)
	room9[r][s][t] = room_p[r][s][t];


  //we first solve room number 9
  rhs_solve(Nx_hall,Ny_hall,Nz,hall_p,hall,landa);
  boundary_wall(Nx_hall,Ny_hall,Nz,hall_p);

  //updating values of hall :D  
  for(int t=0;t<Nz;t++)
    for(int r=(int)(0.2*Nx_room);r<(int)(0.3*Nx_room);r++)
      hall_p[r][0][t] = room_p[r][Ny_room-1][t];
  for(int t=0;t<Nz;t++)
    for(int r=(int)(0.6*Nx_room);r<(int)(0.7*Nx_room);r++)
      hall_p[r][0][t] = room_p[r][Ny_room-1][t];


  //implementing Boundary of doors for room 3
  for(int t=0;t<Nz;t++){
    for(int r=(int)(0.2*Nx_room);r<(int)(0.3*Nx_room);r++)
      hall_p[Nx_hall-r][0][t]= (1.0-6.0*landa)*room3[r][Ny_room-1][t] + landa*(room3[r+1][Ny_room-1][t] + room3[r-1][Ny_room-1][t] + room3[r][Ny_room-2][t] + hall[Nx_hall-r][1][t] + room3[r][Ny_room-1][t+1]+room3[r][Ny_room-1][t-1]);
  }
  for(int t=0;t<Nz;t++){
    for(int r=(int)(0.6*Nx_room);r<(int)(0.7*Nx_room);r++)
      hall_p[Nx_hall-r][0][t]= (1.0-6.0*landa)*room3[r][Ny_room-1][t] + landa*(room3[r+1][Ny_room-1][t] + room3[r-1][Ny_room-1][t] + room3[r][Ny_room-2][t] + hall[Nx_hall-r][1][t] + room3[r][Ny_room-1][t+1]+room3[r][Ny_room-1][t-1]);
  }

  
  for(int r=0;r<Nx_hall;r++)
    for(int s=0;s<Ny_hall;s++)
      for(int t=0;t<Nz;t++)
	hall[r][s][t] = hall_p[r][s][t];



  //we first solve room number 3
  rhs_solve(Nx_room,Ny_room,Nz,room_p,room3,landa);
  boundary_wall(Nx_room,Ny_room,Nz,room_p);

  for(int t=0;t<Nz;t++)
    for(int r=(int)(0.2*Nx_room);r<(int)(0.3*Nx_room);r++)
      room_p[r][Ny_room-1][t] = hall_p[Nx_hall-r][0][t];
  for(int t=0;t<Nz;t++)
    for(int r=(int)(0.6*Nx_room);r<(int)(0.7*Nx_room);r++)
      room_p[r][Ny_room-1][t] = hall_p[Nx_hall-r][0][t];

    
  
  for(int r=0;r<Nx_room;r++)
    for(int s=0;s<Ny_room;s++)
      for(int t=0;t<Nz;t++)
	room3[r][s][t] = room_p[r][s][t];

  //  for(int  t=0;t<Nz;t++)
  //    room9[Nx_room/2][Ny_room/2][Nz/2] = 100.0;
  //  for(int iii=-1;iii<=1;iii++)
  //for(int jjj=-1;jjj<=+1;jjj++)
  //  for(int kkk=-1;kkk<=1;kkk++)
	room9[Nx_room/2][Ny_room/2][Nz/2] = 100.0;


}
//***********************************************************************************
void  rhs_solve (int Nx,int Ny,int Nz,double *** u_new, double *** u, double landa){

  for(int i=1;i<Nx-1;i++)
    for(int l=1;l<Ny-1;l++)
      for(int m=1;m<Nz-1;m++)
	u_new[i][l][m] =(1.0-6.0*landa)*u[i][l][m] + landa*(u[i+1][l][m] + u[i-1][l][m] + u[i][l+1][m] +u[i][l-1][m] + u[i][l][m+1] + u[i][l][m-1]);
  
}

//*********************************************************************************
void  boundary_wall(int Nx_room,int Ny_room,int Nz,double *** u){
  
  //walls in room 9
  for(int r=1;r<Nx_room-1;r++)
    for(int s=1;s<Ny_room-1;s++){
      // BC on z=0 and z=hz
      u[r][s][0] = u[r][s][1];
      u[r][s][Nz-1] = u[r][s][Nz-2];
    }
  for(int r=1;r<Nx_room-1;r++)
    for(int s=1;s<Nz-1;s++){
      //BC on y=0 and y=hy
      u[r][0][s] = u[r][1][s];
      u[r][Ny_room-1][s] = u[r][Nz-2][s];
    }
  for(int r=1;r<Ny_room-1;r++)
    for(int s=1;s<Nz-1;s++){
      // BC on x=0 and x=hx
      u[0][r][s] = u[1][r][s];
      u[Nx_room-1][r][s] = u[Nx_room-2][r][s];
    }
}
//*********************************************************************************
bool check_smell_room3(int Nx,int Ny,int Nz,double ***u){

  bool result = false;
  //  double average=0;
  //for(int i=1;i<Nx-1;i++)
  //   for(int l=1;l<Ny-1;l++)
  //    for(int m=1;m<Nz-1;m++)
  //	average += u[i][l][m];

  //  average = average/((double)Nx-2.0)/((double)Ny-2.0)/((double)Nz-2.0);
  //if(average >= 1)
   
  if(u[(int)(Nx/2)][6][Nz/2] >= 1)
    result = true;

  return result;
  
}










































