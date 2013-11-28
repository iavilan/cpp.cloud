#include <math.h>
#include <time.h>
#include <vector>
#include <string>
#include <cstring>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include <ctime>
#include "electron_manipulations.h"
#include "grid_manipulations.h"
#include "dependencies.h"

#ifndef EPS0 
#define EPS0 8.85e-12
#endif
#define PI 3.141592651
#define CTOM 1.758820088e11
#define EC 1.60217657e-19
#define EM 9.10938291e-31
#define CL 2.99e8

using namespace std;
struct weight;
//THIS FILE CONTAINS DEFINITIONS OF FUNCTIONS THAT ARE USED FOR ELECTRON MANIPULATIONS







//GETTING THE POSITION OF A POINT WHERE ELECTRON CROSSES THE RECTANGULAR WALL.
/*void cross_point(int geoflag, double Rx, double Ry, double x, double y, double Vx, double Vy, double dt, double &xc, double &yc, double &ex, double &ey)
	{
	//GETTING PREVIOUS STEP POSITION
	double x0=x-Vx*dt;
	double y0=y-Vy*dt;
	if(fabs(x0)>Rx || fabs(y0)>Ry)
		{
		printf("Alert: Particle has no cross point %e>%e or %e>%e\n",fabs(x0),Rx,fabs(y0),Ry);
		sleep(1);
		}
	//CROSS POINTS
	xc=0.0;
	yc=0.0;
	//RECTANGULAR GEOMETRY
	if(geoflag==0)
		{
		if(Vx==0)
			{
			xc=x;
			ex=0.0;
			if(Vy>0)
				{
				ey=-1.0;
				yc=Ry;
				}
			else
				{
				ey=1.0;
				yc=-Ry;
				}
			}
		else if (Vy==0)
			{
			yc=y;
			ey=0.0;
			if(Vx>0)
				{
				ex=-1.0;
				xc=Rx;
				}
			else
				{
				ex=1.0;
				xc=-Rx;
				}
			}
		else
			{
			//COEFFICIENTS OF A LINEAR FUNCTION
			double c1=Vy/Vx;
			double c0=y-c1*x;
			double flag=-1;
			//IF A PARTICLE HITS FOR SURE VERTICAL WALLS
			if(fabs(y)<Ry)
				{
				flag=0;
				}
			//IF A PARTICLE HITS FOR SURE HORIZONTAL WALLS
			else if(fabs(x)<Rx)
				{
				flag=1;
				}
			//IF A SITUATION IS NOT THAT OBVIOUS (NEAR THE CORNER)
			else 
				{
				//CALCULATING THE WALL THAT WAS CROSSED
				if(fabs(c1)>(Ry-fabs(y0))/(Rx-fabs(x0)))
					{
					flag=1;
					}
				else
					{
					flag=0;
					}
				}
			//CALCULATING THE CROSS POINT IF A VERTICAL WALL WAS CROSSED
			if(flag==0)
				{
				ey=0.0;
				if(Vx>0)
					{
					ex=-1.0;
					xc=Rx;
					yc=c1*xc+c0;
					}
				else
					{
					ex=1.0;
					xc=-Rx;
					yc=c1*xc+c0;
					}
				}
			//CALCULATING THE CROSS POINT IF A HORIZONTAL WALL WAS CROSSED
			else if(flag==1)
				{
				ex=0.0;
				if(Vy>0)
					{
					ey=-1.0;
					yc=Ry;
					xc=(c0-yc)/c1;
					}
				else
					{
					ey=1.0;
					yc=-Ry;
					xc=(c0-yc)/c1;
					}
				}
			}
		}
	else if(geoflag==1)
		{
		int flagg=0;
		double tang=0.0;
		double tang0=0.0;
		double a=0.0;
		double b=0.0;
		double ratio=Rx/Ry;
		double xr=x/Rx;
		double yr=y/Ry;
		double x0r=x0/Rx;
		double y0r=y0/Ry;
		if(fabs(yr+y0r)<fabs(xr+x0r))
		{
		tang=yr/xr;
		tang0=y0r/x0r;
		a=(yr-y0r)/(xr-x0r);
		b=yr-a*xr;
		}
		else
		{
		flagg=1;
		tang=xr/yr;
		tang0=x0r/y0r;
		a=(xr-x0r)/(yr-y0r);
		b=xr-a*yr;
		}

		double r2=xr*xr+yr*yr;
		double r02=x0r*x0r+y0r*y0r;
		//SOLVING EQUATION FOR yc/xc;
		//double c0=Rx*Rx-b*b;
		double c0=1.0-b*b;
		//double c1=-2.0*a*Rx*Rx;
		double c1=-2.0*a;
		//double c2=Rx*Rx*a*a-b*b;
		double c2=a*a-b*b;
		//DISCRIMINANT
		double disc=c1*c1-4.0*c0*c2;
		double A1=(-c1+sqrt(disc))/(2.0*c0);
		if(   (A1<tang && A1>tang0) || (A1>=tang && A1<=tang0) )
			{
			if(flagg==0)
				{
				xc=Rx/sqrt(1.0+A1*A1);
				yc=xc*A1;
				}
			else
				{
				yc=Ry/sqrt(1.0+A1*A1);
				xc=yc*A1;
				}
			}
		else
			{
			double A2=(-c1-sqrt(disc))/(2.0*c0);
			if(flagg==0)
				{
				xc=Rx/sqrt(1.0+A2*A2);
				yc=xc*A2;
				}
			else
				{
				yc=Ry/sqrt(1.0+A2*A2);
				xc=yc*A2;
				}
			}

		if( (xc<x && xc<x0) || (xc>=x && xc>=x0) )
			{
			xc=-xc;
			}
		if( (yc<y && yc<y0) || (yc>=y && yc>=y0) )
			{
			yc=-yc;
			}
		//CALCULATING THE NORMAL VECTOR
		//printf("x=%e y=%e x0=%e y0=%e xc=%e yc=%e\n",x,y,x0,y0,xc,yc);
		double denumenator=sqrt(pow(Rx,4)*yc*yc+pow(Ry,4)*xc*xc);
		ex=-xc*pow(Ry,2)/denumenator;
		ey=-yc*pow(Rx,2)/denumenator;
		}
	};*/

//CROSSPOINT NICELY DONE FOR ELLIPSE ALSO
void cross_point(int geoflag, double Rx, double Ry, double x, double y, double Vx, double Vy, double dt, double &xc, double &yc, double &ex, double &ey)
	{
	//GETTING PREVIOUS STEP POSITION
	double x0=x-Vx*dt;
	double y0=y-Vy*dt;
	if(fabs(x0)>Rx || fabs(y0)>Ry)
		{
		printf("Alert: Particle has no cross point %e>%e or %e>%e\n",fabs(x0),Rx,fabs(y0),Ry);
		sleep(1);
		}
	//CROSS POINTS
	xc=0.0;
	yc=0.0;
	//RECTANGULAR GEOMETRY
	if(geoflag==2)
		{
		if(Vx==0)
			{
			xc=x;
			ex=0.0;
			if(Vy>0)
				{
				ey=-1.0;
				yc=Ry;
				}
			else
				{
				ey=1.0;
				yc=-Ry;
				}
			}
		else if (Vy==0)
			{
			yc=y;
			ey=0.0;
			if(Vx>0)
				{
				ex=-1.0;
				xc=Rx;
				}
			else
				{
				ex=1.0;
				xc=-Rx;
				}
			}
		else
			{
			//COEFFICIENTS OF A LINEAR FUNCTION
			double c1=Vy/Vx;
			double c0=y-c1*x;
			double flag=-1;
			//IF A PARTICLE HITS FOR SURE VERTICAL WALLS
			if(fabs(y)<Ry)
				{
				flag=0;
				}
			//IF A PARTICLE HITS FOR SURE HORIZONTAL WALLS
			else if(fabs(x)<Rx)
				{
				flag=1;
				}
			//IF A SITUATION IS NOT THAT OBVIOUS (NEAR THE CORNER)
			else 
				{
				//CALCULATING THE WALL THAT WAS CROSSED
				if(fabs(c1)>(Ry-fabs(y0))/(Rx-fabs(x0)))
					{
					flag=1;
					}
				else
					{
					flag=0;
					}
				}
			//CALCULATING THE CROSS POINT IF A VERTICAL WALL WAS CROSSED
			if(flag==0)
				{
				ey=0.0;
				if(Vx>0)
					{
					ex=-1.0;
					xc=Rx;
					yc=c1*xc+c0;
					}
				else
					{
					ex=1.0;
					xc=-Rx;
					yc=c1*xc+c0;
					}
				}
			//CALCULATING THE CROSS POINT IF A HORIZONTAL WALL WAS CROSSED
			else if(flag==1)
				{
				ex=0.0;
				if(Vy>0)
					{
					ey=-1.0;
					yc=Ry;
					xc=(c0-yc)/c1;
					}
				else
					{
					ey=1.0;
					yc=-Ry;
					xc=(c0-yc)/c1;
					}
				}
			}
		}
	else if(geoflag==1)
		{
		int flagg=0;
		double tang=0.0;
		double tang0=0.0;
		double a=0.0;
		double b=0.0;
		double ratio=Ry/Rx;
		double xr=x/Rx;
		double yr=y/Ry;
		double x0r=x0/Rx;
		double y0r=y0/Ry;
		if(fabs(xr-x0r)>fabs(yr-y0r))//fabs(yr+y0r)<fabs(xr+x0r))
		{
		tang=yr/xr;
		tang0=y0r/x0r;
		a=(yr-y0r)/(xr-x0r);
		b=yr-a*xr;
		}
		else
		{
		flagg=1;
		tang=xr/yr;
		tang0=x0r/y0r;
		a=(xr-x0r)/(yr-y0r);
		b=xr-a*yr;
		}

		double r2=xr*xr+yr*yr;
		double r02=x0r*x0r+y0r*y0r;
		//SOLVING EQUATION FOR yc/xc;
		//double c0=Rx*Rx-b*b;
		double c0=b*b-1.0;
		//double c1=-2.0*a*Rx*Rx;
		double c1=2.0*a*b;
		//double c2=Rx*Rx*a*a-b*b;
		double c2=a*a+1.0;
		//DISCRIMINANT
		double disc=sqrt(c1*c1-4.0*c0*c2);
		double A1=(-c1+disc)/(2.0*c2);
		double A2=(-c1-disc)/(2.0*c2);
		//WHEN y=ax+b
		if(flagg==0)
			{
			if((A1<xr && A1>x0r) || (A1>=xr && A1<=x0r))
				{
				xc=A1;
				yc=a*A1+b;
				}
			else if((A2<xr && A2>x0r) || (A2>=xr && A2<=x0r))
				{
				xc=A2;
				yc=a*A2+b;
				}
			}
		else
			{
			if((A1<yr && A1>y0r) || (A1>=yr && A1<=y0r))
			//if(A1/(a*A1+b)>tang && A1/(a*A1+b)<tang0)
				{
				yc=A1;
				xc=a*A1+b;
				}
			else if((A2<yr && A2>y0r) || (A2>=yr && A2<=y0r))
				{
				//double A2=(-c1-disc)/(2.0*c2);
				yc=A2;
				xc=a*A2+b;
				}
			}

		//CALCULATING THE NORMAL VECTOR
		//printf("x0=%e xc=%e x=%e\n",x0r,xc,xr);
		//printf("x=%e y=%e x0=%e y0=%e xc=%e yc=%e\n",x,y,x0,y0,xc,yc);
		
		double denominator=sqrt(xc*xc*ratio*ratio+yc*yc);
		//ex=-ratio*xc/denominator;
		//ey=-yc/denominator;
		ex=-xc*ratio/denominator;
		ey=-yc/denominator;
		xc*=Rx;
		yc*=Ry;
		if(xc==0.0 && yc==0.0)
			{
			printf("Crosspoint in the middle of a pipe\n");
			printf("x0=%e y0=%e\n",x0,y0);
			printf("x=%e y=%e\n",x,y);
			}
		
		}
	};




//CHECKING WETHER THE PARTICLE IS OUTSIDE
bool outside(int geoflag, double Rx, double Ry, double x, double y)
	{
	if(geoflag==1)
		{
		return pow(x/Rx,2)+pow(y/Ry,2)>1.0;
		}
	else if(geoflag==2)
		{
		return (fabs(x)>Rx || fabs(y)>Ry);
		}
	else {return true;}
	}







//KICK ELECTRON VELOCITY
void accelerate_electron(double ax_dt, double ay_dt, double &Vx, double &Vy)
	{
	Vx+=ax_dt;
	Vy+=ay_dt;
	};






//MOVE ELECTRON WITH A GIVEN VELOCITY AND TIME STEP
void shift_electron(double &x, double &y, double Vx, double Vy, double dt)
	{
	x+=Vx*dt;
	y+=Vy*dt;
	};






//ROTATE ELECTRON VELOCITY AROUND THE B FIELD
void rotate_velocity(double &Vx, double &Vz, double cosi, double sini)
	{
	double tVx=Vx*cosi+Vz*sini;
	double tVz=-Vx*sini+Vz*cosi;
	Vx=tVx;
	Vz=tVz;
	}









//GENERATE VELOCITY OF AN ELECTRON (TAKES SEY PARAMETERS AND RANDOM NUMBER GENERATOR)
//sigV - average velocity of new electrons
double gen_velocity(double sigV, gsl_rng *rr)
	{
	double Vborn=0.0;
	//GENERATE VELOCITY WITH TAIL GAUSSIAN DISTRIBUTION
	Vborn=gsl_ran_gaussian_tail(rr,0.0,sigV);
	return Vborn;
	};

double gen_velocity_peak(double sigV, gsl_rng *rr)
	{
	double Vborn1=gsl_ran_gaussian_tail(rr,0.0,sigV);;
	double Vborn2=gsl_ran_gaussian_tail(rr,0.0,sigV);;
	//GENERATE VELOCITY WITH TAIL GAUSSIAN DISTRIBUTION
	//Vborn1=
	return sqrt((Vborn1*Vborn1+Vborn2*Vborn2)/2.0);
	};





//CREATE INITIAL KV CLOUD DISTRIBUTION WITH ZERO VELOCITIES AND GIVEN MACROPARTICLE WEIGHT
void uniform_cloud_ellipse(double Rx, double Ry, vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, vector <double> &lambda_mac, double lambda_mac0, int NUM, gsl_rng *rr)
	{
	//GENERATING A RADIUS-VECTOR

	for(int i=0;i<NUM;++i)
		{
		double r = sqrt(gsl_rng_uniform(rr));
		double phi=8.0*PI*gsl_rng_uniform(rr);
		xe.push_back(Rx*r*sin(phi));
		ye.push_back(Ry*r*cos(phi));
		Vxe.push_back(0.0);
		Vye.push_back(0.0);
		Vze.push_back(0.0);
		lambda_mac.push_back(lambda_mac0);
		}
	};
//CREATE INITIAL KV CLOUD DISTRIBUTION WITH ZERO VELOCITIES AND GIVEN MACROPARTICLE WEIGHT USING ELECTRON STRUCTURE
void uniform_cloud_ellipse(double Rx, double Ry, vector <electron> &ELECTRONS, double lambda_mac0, int NUM, gsl_rng *rr)
	{
	//GENERATING A RADIUS-VECTOR
	electron temp;
	for(int i=0;i<NUM;++i)
		{
		//if(i%300==0)
		//	{
		//	printf("Number of electrons");
		//	}
		//temp.x=0.0;
		double r = sqrt(gsl_rng_uniform(rr));
		double phi=8.0*PI*gsl_rng_uniform(rr);
		temp.x=Rx*r*sin(phi);
		temp.y=Ry*r*cos(phi);	
		temp.Vx=0.0;
		temp.Vy=0.0;
		temp.Vz=0.0;
		temp.lambda_mac=lambda_mac0;
		ELECTRONS.push_back(temp);
		}
	};
//CREATE INITIAL KV CLOUD DISTRIBUTION WITH ZERO VELOCITIES AND GIVEN MACROPARTICLE WEIGHT USING ELECTRON STRUCTURE
void gaussian_cloud_ellipse(double Rx, double Ry, double DX, vector <electron> &ELECTRONS, double lambda_mac0, int NUM, gsl_rng *rr)
	{
	//GENERATING A RADIUS-VECTOR
	electron temp;
	for(int i=0;i<NUM;++i)
		{
		//if(i%300==0)
		//	{
		//	printf("Number of electrons");
		//	}
		//temp.x=0.0;

		temp.x=gsl_ran_gaussian(rr,Rx)+DX;
		temp.y=gsl_ran_gaussian(rr,Ry);	
		while(fabs(temp.x-DX)>3.0*Rx)
			{
			temp.x=gsl_ran_gaussian(rr,Rx)+DX;
			}
		while(fabs(temp.y)>3.0*Ry)
			{
			temp.y=gsl_ran_gaussian(rr,Ry);
			}
		temp.Vx=0.0;
		temp.Vy=0.0;
		temp.Vz=0.0;
		temp.lambda_mac=lambda_mac0;
		ELECTRONS.push_back(temp);
		}
	};
//CREATE AN ELECTRON CLOUD RING
void uniform_cloud_ring(double Rmin, double Rmax, vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, vector <double> &lambda_mac, double lambda_mac0, int NUM, gsl_rng *rr)
	{
	//GENERATING A RADIUS-VECTOR

	for(int i=0;i<NUM;++i)
		{
		double r = Rmax*sqrt(pow(Rmin/Rmax,2)+(1.0-pow(Rmin/Rmax,2))*gsl_rng_uniform(rr));
		double phi=8.0*PI*gsl_rng_uniform(rr);
		xe.push_back(r*sin(phi));
		ye.push_back(r*cos(phi));
		Vxe.push_back(0.0);
		Vye.push_back(0.0);
		Vze.push_back(0.0);
		lambda_mac.push_back(lambda_mac0);
		}
	};
//CREATE AN ELECTRON CLOUD RING USING ELECTRON STRUCTURE
void uniform_cloud_ring(double Rmin, double Rmax, vector <electron> &ELECTRONS, double lambda_mac0, int NUM, gsl_rng *rr)
	{
	//GENERATING A RADIUS-VECTOR
	electron temp;
	for(int i=0;i<NUM;++i)
		{
		double r = Rmax*sqrt(pow(Rmin/Rmax,2)+(1.0-pow(Rmin/Rmax,2))*gsl_rng_uniform(rr));
		double phi=8.0*PI*gsl_rng_uniform(rr);
		temp.x=r*sin(phi);
		temp.y=r*cos(phi);
		temp.Vx=0.0;
		temp.Vy=0.0;
		temp.Vz=0.0;
		temp.lambda_mac=lambda_mac0;
		ELECTRONS.push_back(temp);
		}
	};





//CREATE INITIAL RECTANGULAR DISTRIBUTION
void uniform_cloud_rect(double Rx, double Ry, vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye,vector <double> &Vze, vector <double> &lambda_mac, double lambda_mac0, int NUM, gsl_rng *rr)
	{
	#pragma omp parallel
	{
	#pragma omp for
	for(int i=0;i<NUM;++i)
		{
		double xx = Rx*(gsl_rng_uniform(rr)*2.0-1.0);
		double yy = Ry*(gsl_rng_uniform(rr)*2.0-1.0);
		xe.push_back(xx);
		ye.push_back(yy);
		Vxe.push_back(0.0);
		Vye.push_back(0.0);
		Vze.push_back(0.0);
		lambda_mac.push_back(lambda_mac0);
		}
	}
	};
//CREATE INITIAL RECTANGULAR DISTRIBUTION USING ELECTRON STRUCTURE
void uniform_cloud_rect(double Rx, double Ry, vector <electron> &ELECTRONS, double lambda_mac0, int NUM, gsl_rng *rr)
	{
	electron temp;
	#pragma omp parallel
	{
	#pragma omp for
	for(int i=0;i<NUM;++i)
		{
		temp.x=Rx*(gsl_rng_uniform(rr)*2.0-1.0);
		temp.y=Ry*(gsl_rng_uniform(rr)*2.0-1.0);
		temp.Vx=0.0;
		temp.Vy=0.0;
		temp.Vz=0.0;
		temp.lambda_mac=lambda_mac0;
		ELECTRONS.push_back(temp);
		}
	}
	};





//THIS FUNCTION IS THE SIMPLIEST MERGING FUNCTION. IT REMOVES RANDOMLY PARTICLES FROM THE DISTRIBUTION AND REWEIGHTS ALL OF THE LEFT ONES PROPORTIONALLY
void simple_merge(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, vector <double> &lambda_mac, int MAXNUM, gsl_rng *rr)
	{
	int SIZE=xe.size();
	double ratio=(double)SIZE/MAXNUM;
	if(SIZE>MAXNUM)
		{
		//NEW MACROPARTICLE WEIGHT
		double neweight=lambda_mac.at(0);
		//IF SIZE OF THE CLOUD IS BIGGER THAN THE ALLOWED MAXIMUM

		{

		for(int i=0;i<SIZE-MAXNUM;++i)
		//while(xe.size()>MAXNUM)
			{
			int todelete=floor(xe.size()*gsl_rng_uniform(rr));
			xe.erase(xe.begin()+todelete);
			ye.erase(ye.begin()+todelete);
			Vxe.erase(Vxe.begin()+todelete);
			Vye.erase(Vye.begin()+todelete);
			Vze.erase(Vze.begin()+todelete);
			lambda_mac.erase(lambda_mac.begin()+todelete);
			}
		for(int i=0;i<MAXNUM;++i)
			{
			lambda_mac[i]*=ratio;
			}
		}
		}
	};
//THIS FUNCTION IS THE SIMPLIEST MERGING FUNCTION. IT REMOVES RANDOMLY PARTICLES FROM THE DISTRIBUTION AND REWEIGHTS ALL OF THE LEFT ONES PROPORTIONALLY. FOR ELECTRON STRUCTUTE
void simple_merge(vector <electron> &ELECTRONS, int MAXNUM, gsl_rng *rr)
	{
	int SIZE=ELECTRONS.size();
	double ratio=(double)SIZE/MAXNUM;
	if(SIZE>MAXNUM)
		{
		int var=SIZE;
		//IF SIZE OF THE CLOUD IS LARGER THAN THE ALLOWED MAXIMUM
		for(int i=0;i<SIZE-MAXNUM;++i)
		//while(xe.size()>MAXNUM)
			{
			int todelete=gsl_rng_uniform_int(rr,var);
			//ELECTRONS.erase(ELECTRONS.begin()+todelete);
			ELECTRONS[todelete]=ELECTRONS.back();
			ELECTRONS.pop_back();
			var--;
			}
		for(int i=0;i<MAXNUM;++i)
			{
			ELECTRONS[i].lambda_mac*=ratio;
			}
		
		}
	};
//THIS FUNCTION SPLITS THE PARTICLES BASED ON THE SYMMETRY OF A SYSTEM (ONLY VALID FOR SYMMETRICAL CASES)
void simple_split(vector <electron> &ELECTRONS, int MINNUM, gsl_rng *rr)
	{
	int SIZE=ELECTRONS.size();
	if(SIZE<MINNUM*0.9)
		{
		//REWEIGHTING COEFFICIENT
		double ratio=(double)SIZE/MINNUM;
		//THE DIFFERENCE BETWEEN A GIVEN NUMBER OF ELECTRONS AND SET MAXIMUM
		int DELTA1=MINNUM-SIZE;
		//LENGTH OF A NUMERIC INTERVAL IN WHICH WE ARE LOOKING FOR SOURCE ELECTRONS
		int DELTA2=SIZE/DELTA1;
		//IF SIZE OF THE CLOUD IS LARGER THAN THE ALLOWED MAXIMUM
		for(int i=0;i<DELTA1;++i)
			{
			int index=i*DELTA2+gsl_rng_uniform_int(rr,DELTA2);
			//GENERATE NUMBER -1 or +1 FOR HORIZONTAL COORDINATE
			int xcoef=gsl_rng_uniform_int(rr,2 )*2-1;
			//INITIAL VALUE IS -1 FOR VERTICAL COORDINATE
			int ycoef=-1;
			//IT CAN BE CHANGED TO +1 IF HORIZONTAL COORDINATE HAS -1
			if(xcoef==-1)
				{
				ycoef=gsl_rng_uniform_int(rr,2 )*2-1;
				}
			//MAKING A NEW ELECTRON
			electron tempo;
			tempo.x=ELECTRONS[index].x*xcoef;
			tempo.y=ELECTRONS[index].y*ycoef;
			tempo.Vx=ELECTRONS[index].Vx*xcoef;
			tempo.Vy=ELECTRONS[index].Vy*ycoef;
			tempo.Vz=xcoef*ELECTRONS[index].Vz;
			tempo.lambda_mac=ELECTRONS[index].lambda_mac;
			ELECTRONS.push_back(tempo);
			}
		SIZE=ELECTRONS.size();
		for(int i=0;i<SIZE;++i)
			{
			ELECTRONS[i].lambda_mac*=ratio;
			}
		
		}
	};






//FUNCTION THAT SIMPLY REMOVES ALL THE ELECTRONS IN RECTANGULAR GEOMETRY
void remove_all_out_rect(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, vector <double> &lambda_mac, double Rx, double Ry)
	{
	int SIZE=xe.size();
	vector <double> xet,yet,Vxet,Vyet,Vzet,lambda_mact;
	/*for(int i=SIZE-1;i>=0;i--)
		{
		if(fabs(xe[i])>Rx || fabs(ye[i])>Ry)
			{
			xe.erase(xe.begin()+i);
			ye.erase(ye.begin()+i);
			Vxe.erase(Vxe.begin()+i);
			Vye.erase(Vye.begin()+i);
			Vze.erase(Vze.begin()+i);
			lambda_mac.erase(lambda_mac.begin()+i);
			}
		SIZE=xe.size();
		}*/
	for(int i=SIZE-1;i>=0;i--)
		{
		if(fabs(xe[i])<Rx && fabs(ye[i])<Ry)
			{
			xet.push_back(xe[i]);
			yet.push_back(ye[i]);
			Vxet.push_back(Vxe[i]);
			Vyet.push_back(Vye[i]);
			Vzet.push_back(Vze[i]);
			lambda_mact.push_back(lambda_mac[i]);
			}
		SIZE=xe.size();
		}
	xe=xet;
	ye=yet;
	Vxe=Vxet;
	Vye=Vyet;
	Vze=Vzet;
	lambda_mac=lambda_mact;
	};











//FUNCTION THAT SIMPLY REMOVES ALL THE ELECTRONS IN RECTANGULAR GEOMETRY (NOT VERY CORRECT BECAUSE DOESN'T INCLUDE THE EFFECT OF A FIELD ON A WAY BACK)
void reflect_all_out_rect(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, double Rx, double Ry, double absorption, gsl_rng *rr)
	{
	int SIZE=xe.size();
	for(int i=0;i<SIZE;++i)
		{
		double xabs=fabs(xe[i]);
		double yabs=fabs(ye[i]);
		if(xabs>Rx && yabs>Ry)
			{
			double dx=xabs-Rx;
			double dy=yabs-Ry;
			xe[i]*=(Rx-dx)/xabs;
			Vxe[i]=-Vxe[i]*absorption;
			ye[i]*=(Ry-dy)/yabs;
			Vye[i]=-Vye[i]*absorption;
			}
		else if(xabs>Rx)
			{
			double dx=xabs-Rx;
			xe[i]*=(Rx-dx)/xabs;
			Vxe[i]=-Vxe[i]*absorption;
			}
		else if(yabs>Ry)
			{
			double dy=yabs-Ry;
			ye[i]*=(Ry-dy)/yabs;
			Vye[i]=-Vye[i]*absorption;
			}
		}
	};












//REMOVE / PRODUCE SECONDARY ELECTRONS
//ENERGY IS V^{2}
//ANOTHER ADDITIONAL PARAMETER IS A PRECALCULATED ARRAY cos_sin TO PRODUCE ELECTRONS WITH cos_thetha ANGULAR DISTRIBUTION 
void wall_effects(vector <electron> &ELECTRONS, 
double Rx, double Ry, 
double dmax, double Emax, double refl0, double Erefl, double sigma_v, 
vector <double> cos_sin, double dt, int geoflag, gsl_rng *rr)
	{
	int SIZE=ELECTRONS.size();
	int cs_size=cos_sin.size()/2;
	//NORMALIZATION
	Erefl*=2.0*CTOM;
	Emax*=2.0*CTOM;
	//VECTORS TO STORE PARTICLES THAT ARE NOT CHANGED / DON'T PARTICIPATE IN THE WALL INTERACTION
	//vector <double> xei,yei,Vxei,Vyei,Vzei,lambda_maci;
	//VECTORS TO STORE PARTICLES THAT ARE OUTSIDE THE PIPE WITH ADDITIONAL PARAMETER ENERGY, NUMBER OF EMITTED ELECTRONS, REFLECTION COEFFICIENT, CROSS POINTS AND NORMALS TO THE SURFACE, dt2 IS A PART OF TIME STEP FROM SURFACE TO THE PARTICLE INSIDE THE SURFACE / EACH ARTICLE HAS A CHANCE TO BE REMOVED
	vector <double> energyo, xc, yc, ex,ey;
	vector <electron> ELECTRONSO;
	//REFLECTION STATUS ABSORBED OR REFLECTED
	vector <int> reflo, nemito;
	//VECTORS TO STORE PARTICLES THAT WERE PRODUCED
	//vector <double> xep,yep,Vxep,Vyep,Vzep,lambda_macp;
	vector <electron> ELECTRONSP;
	//printf("First loop\n");
	//WE CONSTRUCT VECTORS OF PARTICLE OUTSIDE THE GRID AND INSIDE THE GRID
	for(int i=0;i<SIZE;++i)
		{
		//WHILE LOOP IS USED FOR THE CASE OF A CORNER (PROBABLY VERY RARE) WHEN MULTIPLE REFLECTION CAN HAPPEN
		if(outside(geoflag, Rx, Ry, ELECTRONS[i].x, ELECTRONS[i].y))
			{
			//PARTICLES OUTSIDE
			ELECTRONSO.push_back(ELECTRONS[i]);
			//IT IS ASSUMED THAT ENERGY CONSTANTS ARE NORMALIZED
			//energyo.push_back(pow(ELECTRONS[i].Vx,2)+pow(ELECTRONS[i].Vy,2)+pow(ELECTRONS[i].Vz,2));
			//REMOVED THE Z COMPONENT FROM ENERGY CALCULATION
			energyo.push_back(pow(ELECTRONS[i].Vx,2)+pow(ELECTRONS[i].Vy,2)+pow(ELECTRONS[i].Vz,2));
			//CROSS POINTS AND NORM VECTORS
			xc.push_back(0.0);
			yc.push_back(0.0);
			ex.push_back(0.0);
			ey.push_back(0.0);
			//ASSIGNING CROSSING POINTS AND NORMAL VECTORS
			cross_point(geoflag, Rx, Ry, ELECTRONS[i].x, ELECTRONS[i].y, ELECTRONS[i].Vx, ELECTRONS[i].Vy, dt, xc.back(), yc.back(), ex.back(), ey.back());
			//REPLACE PARTICLES
			ELECTRONS[i]=ELECTRONS.back();
			//ERASING LAST ELEMENTS
			ELECTRONS.pop_back();
			SIZE--;
			i--;
			//printf("%e ", energyo.back()*EM/2/EC);
			//CALCULATE IF AN ELECTRON IS REFLECTED
			if(refl(refl0, Erefl, energyo.back())>gsl_rng_uniform(rr))
				{
				reflo.push_back(1);
				//printf("zaebok \n");
				}
			else
				{
				reflo.push_back(0);
				}
			//NUMBER OF ELECTRONS TO EMIT
			nemito.push_back(n_emitted(sey(Emax, dmax, 1.2, energyo.back()),rr));
			//ADD SEARCH FOR LEFT OVER
			}
		}
	int SIZEE=ELECTRONSO.size();
	//NOW WE CONSTRUCT A NEW "OUTPUT" VECTOR OF PARTICLES CONSISTING OF INTERNAL PARTICLES, REFLECTED PARTICLES AND SECONDARIES
	//1. CALCULATE SECONDARIES
	for(int i=0;i<SIZEE;++i)
		{
		//IF THE PARTICLE IS REFLECTED, THEN IT IS ADDED TO THE NEWBORN VECTOR WITH REFLECTED VELOCITIES AND COORDINATES (NOT TAKING INTO ACCOUNT THE EXTERNAL FORCE)
		if(reflo[i]==1)
			{
			//REFLECT AN ELECTRON AND ATTACHING IT TO VECTOR OF SECONDARY PARTICLES
			reflect(ELECTRONSP, ELECTRONSO[i].x, ELECTRONSO[i].y, ELECTRONSO[i].Vx, ELECTRONSO[i].Vy,ELECTRONSO[i].Vz, ELECTRONSO[i].lambda_mac, xc[i], yc[i],ex[i], ey[i], dt, Rx, Ry, geoflag);
			}
		//PRODUCTION OF TRUE SECONDARIES
		for(int j=0;j<nemito[i];++j)
			{
			//CHOSING RANDOMLY THE ANGLE OF EMISSION
			int angle_number=gsl_rng_uniform_int(rr,cs_size);
			//PERPENDICULAR COMPONENT OF THE VELOCITY
			double ep=cos_sin.at(angle_number);
			//PARALLEL TO SURFACE COMPONENT OF VELOCITY
			double eh=cos_sin.at(angle_number+cs_size);
			//GENERATING VELOCITY WITH THE POISSON DISTRIBUTION
			double g_velocity=gen_velocity(sigma_v, rr);
			//HERE WE MAKE A HYBRID MODEL ELECTRONS ARE PRODUCED UNCORRELATED, BUT THE ENERGY OF EACH
			//ELECTRON CAN NOT EXCEED THE IMPACT ENERGY (IN TOTAL STILL POSSIBLE TO EXCEED)
			int checker=0;
			/*while(g_velocity*g_velocity>Vxeo[i]*Vxeo[i]+Vyeo[i]*Vyeo[i]+Vzeo[i]*Vzeo[i])
				{
				checker++;
				if(checker>100 && checker%100==0)
					{
					printf("Alert: More than 100 iterations %e %e\n",g_velocity*g_velocity,ELECTRONSO[i].Vx*ELECTRONSO[i].Vx+ELECTRONSO[i].Vy*ELECTRONSO[i].Vy+ELECTRONSO[i].Vz*ELECTRONSO[i].Vz);
					}
				g_velocity=gen_velocity(sigma_v, rr);
				}*/
			double vgp=g_velocity*ep;
			double vgh=g_velocity*eh;
			electron temp;
			temp.Vx=vgp*ex[i]+vgh*ey[i];
			temp.Vy=vgp*ey[i]-vgh*ex[i];
			//THIS MOMENT SHOULD BE WORKED OUT IN DETAILS TO INCLUDE THE EFFECT OF LOCAL FIELD
			temp.x=xc[i]*0.9999999;
			temp.y=yc[i]*0.9999999;
			temp.Vz=0.0;
			temp.lambda_mac=ELECTRONSO[i].lambda_mac;
			ELECTRONSP.push_back(temp);
			/*if(fabs(Vxep.back())>1e10 || fabs(Vyep.back())>1e10 || fabs(xep.back())>Rx || fabs(yep.back())>Ry)
				{
				printf("%e %e %e %e %e %e\n",xep.back(),yep.back(),ex[i],ey[i],xc[i],yc[i]);
				sleep(0.5);
				}*/
			}
		}
	//MERGING ARRAYS OF INTERNAL AND NEWLY PRODUCED ELECTRONS
	ELECTRONS.insert(ELECTRONS.end(),ELECTRONSP.begin(),ELECTRONSP.end());
	//ASSIGNING ALL THIS TO OUT INPUT ARRAY
	//xe=xei;
	//ye=yei;
	//Vxe=Vxei;
	//Vye=Vyei;
	//Vze=Vzei;
	//lambda_mac=lambda_maci;
	//printf("Left %d\n",xe.size());
	};

//WALL EFFECTS WITH VELOCITY GENERATION HAVING ZERO PROBABILITY FOR ZERO ENERGY
//REMOVE / PRODUCE SECONDARY ELECTRONS
//ENERGY IS V^{2}
//ANOTHER ADDITIONAL PARAMETER IS A PRECALCULATED ARRAY cos_sin TO PRODUCE ELECTRONS WITH cos_thetha ANGULAR DISTRIBUTION 
void wall_effects_peak(vector <electron> &ELECTRONS, 
double Rx, double Ry, 
double dmax, double Emax, double refl0, double Erefl, double sigma_v, double redif,
vector <double> cos_sin, double dt, int geoflag, gsl_rng *rr,
double &energy_deposition,
//COEFFICIENT FOR SECONDARY EMISSION
double sss)
	{
	//CALCULATE PEAK OF TRUE SECONDARIES
	double dsey=dmax-redif;
	int SIZE=ELECTRONS.size();
	int cs_size=cos_sin.size()/2;
	//NORMALIZATION
	Erefl*=2.0*CTOM;
	Emax*=2.0*CTOM;
	double Er=5.0*2.0*CTOM;
	//VECTORS TO STORE PARTICLES THAT ARE NOT CHANGED / DON'T PARTICIPATE IN THE WALL INTERACTION
	//vector <double> xei,yei,Vxei,Vyei,Vzei,lambda_maci;
	//VECTORS TO STORE PARTICLES THAT ARE OUTSIDE THE PIPE WITH ADDITIONAL PARAMETER ENERGY, NUMBER OF EMITTED ELECTRONS, REFLECTION COEFFICIENT, CROSS POINTS AND NORMALS TO THE SURFACE, dt2 IS A PART OF TIME STEP FROM SURFACE TO THE PARTICLE INSIDE THE SURFACE / EACH ARTICLE HAS A CHANCE TO BE REMOVED
	vector <double> energyo, xc, yc, ex,ey;
	vector <electron> ELECTRONSO;
	//REFLECTION STATUS ABSORBED OR REFLECTED
	vector <int> reflo, nemito;
	//STORES NUMBER OF REDIFFUSED ELECTRONS (ACTUALLY 1 OR 0 MULTIPLIED BY THE IMPACT ENERGY)
	vector <double> nemitred;
	//VECTORS TO STORE PARTICLES THAT WERE PRODUCED
	//vector <double> xep,yep,Vxep,Vyep,Vzep,lambda_macp;
	vector <electron> ELECTRONSP;
	energy_deposition=0.0;
	//printf("First loop\n");
	//WE CONSTRUCT VECTORS OF PARTICLE OUTSIDE THE GRID AND INSIDE THE GRID
	for(int i=0;i<SIZE;++i)
		{
		//WHILE LOOP IS USED FOR THE CASE OF A CORNER (PROBABLY VERY RARE) WHEN MULTIPLE REFLECTION CAN HAPPEN
		if(outside(geoflag, Rx, Ry, ELECTRONS[i].x, ELECTRONS[i].y))
			{
			//PARTICLES OUTSIDE
			ELECTRONSO.push_back(ELECTRONS[i]);
			//IT IS ASSUMED THAT ENERGY CONSTANTS ARE NORMALIZED
			//energyo.push_back(pow(ELECTRONS[i].Vx,2)+pow(ELECTRONS[i].Vy,2)+pow(ELECTRONS[i].Vz,2));
			//REMOVED THE Z COMPONENT FROM ENERGY CALCULATION
			energyo.push_back(pow(ELECTRONS[i].Vx,2)+pow(ELECTRONS[i].Vy,2)+pow(ELECTRONS[i].Vz,2));
			//WE ADD THE ENERGY THAT ELECTRON HAS OUTSIDE THE WALL
			energy_deposition+=energyo.back();
			//CROSS POINTS AND NORM VECTORS
			xc.push_back(0.0);
			yc.push_back(0.0);
			ex.push_back(0.0);
			ey.push_back(0.0);
			//ASSIGNING CROSSING POINTS AND NORMAL VECTORS
			cross_point(geoflag, Rx, Ry, ELECTRONS[i].x, ELECTRONS[i].y, ELECTRONS[i].Vx, ELECTRONS[i].Vy, dt, xc.back(), yc.back(), ex.back(), ey.back());
			//REPLACE PARTICLES
			ELECTRONS[i]=ELECTRONS.back();
			//ERASING LAST ELEMENTS
			ELECTRONS.pop_back();
			SIZE--;
			i--;
			//printf("%e ", energyo.back()*EM/2/EC);
			//CALCULATE IF AN ELECTRON IS REFLECTED
			if(refl(refl0, Erefl, energyo.back())>gsl_rng_uniform(rr))
				{
				reflo.push_back(1);
				//printf("zaebok \n");
				}
			else
				{
				reflo.push_back(0);
				}
			//NUMBER OF TRUE SECONDARY ELECTRONS TO EMIT
			double locsey=sey(Emax, dsey, sss, energyo.back());
			nemito.push_back(n_emitted(locsey,rr));
			//NUMBER OF REDIFFUSED ELECTRONS TO EMIT ASSUMING THAT Er=5eV
			nemitred.push_back(energyo.back()*n_emitted(redif*(1.0-exp(-energyo.back()/Er)),rr));
			/*if(redif>0.0 && locsey<=redif)
				{
				nemitred.push_back(energyo.back()*n_emitted(locsey,rr));
				locsey=0.0;
				}
			else if(redif>0.0 && locsey>redif)
				{
				locsey-=redif;
				//NUMBER OF EMITTED REDIFFUSED ELECTRONS

				}
			else
				{
				nemitred.push_back(0.0);
				}*/
			//ADD SEARCH FOR LEFT OVER
			}
		}
	int SIZEE=ELECTRONSO.size();
	//NOW WE CONSTRUCT A NEW "OUTPUT" VECTOR OF PARTICLES CONSISTING OF INTERNAL PARTICLES, REFLECTED PARTICLES AND SECONDARIES
	//1. CALCULATE SECONDARIES
	for(int i=0;i<SIZEE;++i)
		{
		//IF THE PARTICLE IS REFLECTED, THEN IT IS ADDED TO THE NEWBORN VECTOR WITH REFLECTED VELOCITIES AND COORDINATES (NOT TAKING INTO ACCOUNT THE EXTERNAL FORCE)
		if(reflo[i]==1)
			{
			//REFLECT AN ELECTRON AND ATTACHING IT TO VECTOR OF SECONDARY PARTICLES
			reflect(ELECTRONSP, ELECTRONSO[i].x, ELECTRONSO[i].y, ELECTRONSO[i].Vx, ELECTRONSO[i].Vy,ELECTRONSO[i].Vz, ELECTRONSO[i].lambda_mac, xc[i], yc[i],ex[i], ey[i], dt, Rx, Ry, geoflag);
			}
		//PRODUCTION OF TRUE SECONDARIES
		for(int j=0;j<nemito[i];++j)
			{
			//CHOSING RANDOMLY THE ANGLE OF EMISSION
			int angle_number=gsl_rng_uniform_int(rr,cs_size);
			//PERPENDICULAR COMPONENT OF THE VELOCITY
			double ep=cos_sin.at(angle_number);
			//PARALLEL TO SURFACE COMPONENT OF VELOCITY
			double eh=cos_sin.at(angle_number+cs_size);
			//ERATING VELOCITY WITH THE POISSON DISTRIBUTION
			double g_velocity=0.0;
			double absvalue=sqrt(ELECTRONSO[i].Vx*ELECTRONSO[i].Vx+ELECTRONSO[i].Vy*ELECTRONSO[i].Vy+ELECTRONSO[i].Vz*ELECTRONSO[i].Vz);
			//IF VELOCITY OF INCIDENT ELECTRON IS BIGGER THAN sigma_v
			if(absvalue>sigma_v)
			{
			do { g_velocity=gen_velocity_peak(sigma_v, rr);}
			while(g_velocity>absvalue);
			}
			else
			{
			do { g_velocity=gen_velocity_peak(absvalue/2, rr);}
			while(g_velocity>absvalue);
			}
			//EMITTED ELECTRONS TAKE THE ENERGY
			energy_deposition-=g_velocity*g_velocity;
			//HERE WE MAKE A HYBRID MODEL ELECTRONS ARE PRODUCED UNCORRELATED, BUT THE ENERGY OF EACH
			//ELECTRON CAN NOT EXCEED THE IMPACT ENERGY (IN TOTAL STILL POSSIBLE TO EXCEED)
			int checker=0;
			/*while(g_velocity*g_velocity>Vxeo[i]*Vxeo[i]+Vyeo[i]*Vyeo[i]+Vzeo[i]*Vzeo[i])
				{
				checker++;
				if(checker>100 && checker%100==0)
					{
					printf("Alert: More than 100 iterations %e %e\n",g_velocity*g_velocity,ELECTRONSO[i].Vx*ELECTRONSO[i].Vx+ELECTRONSO[i].Vy*ELECTRONSO[i].Vy+ELECTRONSO[i].Vz*ELECTRONSO[i].Vz);
					}
				g_velocity=gen_velocity(sigma_v, rr);
				}*/
			double vgp=g_velocity*ep;
			double vgh=g_velocity*eh;
			electron temp;
			temp.Vx=vgp*ex[i]+vgh*ey[i];
			temp.Vy=vgp*ey[i]-vgh*ex[i];
			//THIS MOMENT SHOULD BE WORKED OUT IN DETAILS TO INCLUDE THE EFFECT OF LOCAL FIELD
			temp.x=xc[i]*0.99999;
			temp.y=yc[i]*0.99999;
			temp.Vz=0.0;
			temp.lambda_mac=ELECTRONSO[i].lambda_mac;
			ELECTRONSP.push_back(temp);
			/*if(fabs(Vxep.back())>1e10 || fabs(Vyep.back())>1e10 || fabs(xep.back())>Rx || fabs(yep.back())>Ry)
				{
				printf("%e %e %e %e %e %e\n",xep.back(),yep.back(),ex[i],ey[i],xc[i],yc[i]);
				sleep(0.5);
				}*/
			}
		//MODELING REDIFFUSEDELECTRONS
		if(nemitred[i]!=0.0)
			{
			//CHOSING RANDOMLY THE ANGLE OF EMISSION
			int angle_number=gsl_rng_uniform_int(rr,cs_size);
			//PERPENDICULAR COMPONENT OF THE VELOCITY
			double ep=cos_sin.at(angle_number);
			//PARALLEL TO SURFACE COMPONENT OF VELOCITY
			double eh=cos_sin.at(angle_number+cs_size);
			double g_velocity=sqrt(gsl_rng_uniform(rr)*nemitred[i]);
			energy_deposition-=g_velocity*g_velocity;
			double vgp=g_velocity*ep;
			double vgh=g_velocity*eh;
			electron temp;
			temp.Vx=vgp*ex[i]+vgh*ey[i];
			temp.Vy=vgp*ey[i]-vgh*ex[i];
			//THIS MOMENT SHOULD BE WORKED OUT IN DETAILS TO INCLUDE THE EFFECT OF LOCAL FIELD
			temp.x=xc[i]*0.99999;
			temp.y=yc[i]*0.99999;
			temp.Vz=0.0;
			temp.lambda_mac=ELECTRONSO[i].lambda_mac;
			ELECTRONSP.push_back(temp);
			}
		}
	//MERGING ARRAYS OF INTERNAL AND NEWLY PRODUCED ELECTRONS
	ELECTRONS.insert(ELECTRONS.end(),ELECTRONSP.begin(),ELECTRONSP.end());
	//ASSIGNING ALL THIS TO OUT INPUT ARRAY
	//xe=xei;
	//ye=yei;
	//Vxe=Vxei;
	//Vye=Vyei;
	//Vze=Vzei;
	//lambda_mac=lambda_maci;
	//printf("Left %d\n",xe.size());
	energy_deposition*=EM*ELECTRONS[0].lambda_mac/EC/2/dt;
	};

void wall_effects(vector <double> &xe, vector <double> &ye, 
vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, 
vector <double> &lambda_mac, 
double Rx, double Ry, 
double dmax, double Emax, double refl0, double Erefl, double sigma_v, 
vector <double> cos_sin, double dt, int geoflag, gsl_rng *rr)
	{
	int SIZE=xe.size();
	int cs_size=cos_sin.size()/2;
	//NORMALIZATION
	Erefl*=2.0*CTOM;
	//VECTORS TO STORE PARTICLES THAT ARE NOT CHANGED / DON'T PARTICIPATE IN THE WALL INTERACTION
	//vector <double> xei,yei,Vxei,Vyei,Vzei,lambda_maci;
	//VECTORS TO STORE PARTICLES THAT ARE OUTSIDE THE PIPE WITH ADDITIONAL PARAMETER ENERGY, NUMBER OF EMITTED ELECTRONS, REFLECTION COEFFICIENT, CROSS POINTS AND NORMALS TO THE SURFACE, dt2 IS A PART OF TIME STEP FROM SURFACE TO THE PARTICLE INSIDE THE SURFACE / EACH ARTICLE HAS A CHANCE TO BE REMOVED
	vector <double> xeo,yeo,Vxeo,Vyeo,Vzeo,lambda_maco, energyo, xc, yc, ex,ey;
	//REFLECTION STATUS ABSORBED OR REFLECTED
	vector <int> reflo, nemito;
	//VECTORS TO STORE PARTICLES THAT WERE PRODUCED
	vector <double> xep,yep,Vxep,Vyep,Vzep,lambda_macp;
	//printf("First loop\n");
	//WE CONSTRUCT VECTORS OF PARTICLE OUTSIDE THE GRID AND INSIDE THE GRID
	for(int i=0;i<SIZE;++i)
		{
		//WHILE LOOP IS USED FOR THE CASE OF A CORNER (PROBABLY VERY RARE) WHEN MULTIPLE REFLECTION CAN HAPPEN
		if(outside(geoflag, Rx, Ry, xe[i], ye[i]))
			{
			//PARTICLES OUTSIDE
			xeo.push_back(xe[i]);
			yeo.push_back(ye[i]);
			Vxeo.push_back(Vxe[i]);
			Vyeo.push_back(Vye[i]);
			Vzeo.push_back(Vze[i]);
			lambda_maco.push_back(lambda_mac[i]);
			//IT IS ASSUMED THAT ENERGY CONSTANTS ARE NORMALIZED
			energyo.push_back(pow(Vxe[i],2)+pow(Vye[i],2)+pow(Vze[i],2));
			//CROSS POINTS AND NORM VECTORS
			xc.push_back(0.0);
			yc.push_back(0.0);
			ex.push_back(0.0);
			ey.push_back(0.0);
			//ASSIGNING CROSSING POINTS AND NORMAL VECTORS
			cross_point(geoflag, Rx, Ry, xe[i], ye[i], Vxe[i], Vye[i], dt, xc.back(), yc.back(), ex.back(), ey.back());
			//REPLACE PARTICLES
			xe[i]=xe.back();
			ye[i]=ye.back();
			Vxe[i]=Vxe.back();
			Vye[i]=Vye.back();
			Vze[i]=Vze.back();
			lambda_mac[i]=lambda_mac.back();
			//ERASING LAST ELEMENTS
			xe.pop_back();
			ye.pop_back();
			Vxe.pop_back();
			Vye.pop_back();
			Vze.pop_back();
			lambda_mac.pop_back();
			SIZE--;
			i--;
			//printf("%e ", energyo.back()*EM/2/EC);
			//CALCULATE IF AN ELECTRON IS REFLECTED
			if(refl(refl0, Erefl, energyo.back())>gsl_rng_uniform(rr))
				{
				reflo.push_back(1);
				//printf("zaebok \n");
				}
			else
				{
				reflo.push_back(0);
				}
			//NUMBER OF ELECTRONS TO EMIT
			nemito.push_back(n_emitted(sey(Emax, dmax, 1.2, energyo.back()),rr));
			//ADD SEARCH FOR LEFT OVER
			}
		}
	int SIZEE=xeo.size();
	//NOW WE CONSTRUCT A NEW "OUTPUT" VECTOR OF PARTICLES CONSISTING OF INTERNAL PARTICLES, REFLECTED PARTICLES AND SECONDARIES
	//1. CALCULATE SECONDARIES
	for(int i=0;i<SIZEE;++i)
		{
		//IF THE PARTICLE IS REFLECTED, THEN IT IS ADDED TO THE NEWBORN VECTOR WITH REFLECTED VELOCITIES AND COORDINATES (NOT TAKING INTO ACCOUNT THE EXTERNAL FORCE)
		if(reflo[i]==1)
			{
			//REFLECT AN ELECTRON AND ATTACHING IT TO VECTOR OF SECONDARY PARTICLES
			reflect(xep,yep,Vxep,Vyep,Vzep,lambda_macp,xeo[i], yeo[i], Vxeo[i], Vyeo[i], lambda_maco[i],xc[i], yc[i],ex[i], ey[i], dt, Rx, Ry, geoflag);
			}
		//PRODUCTION OF TRUE SECONDARIES
		for(int j=0;j<nemito[i];++j)
			{
			//CHOSING RANDOMLY THE ANGLE OF EMISSION
			int angle_number=gsl_rng_uniform_int(rr,cs_size);
			//PERPENDICULAR COMPONENT OF THE VELOCITY
			double ep=cos_sin.at(angle_number);
			//PARALLEL TO SURFACE COMPONENT OF VELOCITY
			double eh=cos_sin.at(angle_number+cs_size);
			//GENERATING VELOCITY WITH THE POISSON DISTRIBUTION
			double g_velocity=gen_velocity(sigma_v, rr);
			//HERE WE MAKE A HYBRID MODEL ELECTRONS ARE PRODUCED UNCORRELATED, BUT THE ENERGY OF EACH
			//ELECTRON CAN NOT EXCEED THE IMPACT ENERGY (IN TOTAL STILL POSSIBLE TO EXCEED)
			int checker=0;
			while(g_velocity*g_velocity>Vxeo[i]*Vxeo[i]+Vyeo[i]*Vyeo[i]+Vzeo[i]*Vzeo[i])
				{
				checker++;
				if(checker>100 && checker%100==0)
					{
					printf("Alert: More than 100 iterations %e %e\n",g_velocity*g_velocity,Vxeo[i]*Vxeo[i]+Vyeo[i]*Vyeo[i]+Vzeo[i]*Vzeo[i]);
					}
				g_velocity=gen_velocity(sigma_v, rr);
				}
			double vgp=g_velocity*ep;
			double vgh=g_velocity*eh;
			Vxep.push_back(vgp*ex[i]+vgh*ey[i]);
			Vyep.push_back(vgp*ey[i]-vgh*ex[i]);
			//THIS MOMENT SHOULD BE WORKED OUT IN DETAILS TO INCLUDE THE EFFECT OF LOCAL FIELD
			xep.push_back(xc[i]*0.99999);
			yep.push_back(yc[i]*0.99999);
			Vzep.push_back(0.0);
			lambda_macp.push_back(lambda_maco[i]);
			if(fabs(Vxep.back())>1e10 || fabs(Vyep.back())>1e10 || fabs(xep.back())>Rx || fabs(yep.back())>Ry)
				{
				printf("%e %e %e %e %e %e\n",xep.back(),yep.back(),ex[i],ey[i],xc[i],yc[i]);
				sleep(0.5);
				}
			}
		}
	//MERGING ARRAYS OF INTERNAL AND NEWLY PRODUCED ELECTRONS
	xe.insert(xe.end(),xep.begin(),xep.end());
	ye.insert(ye.end(),yep.begin(),yep.end());
	Vxe.insert(Vxe.end(),Vxep.begin(),Vxep.end());
	Vye.insert(Vye.end(),Vyep.begin(),Vyep.end());
	Vze.insert(Vze.end(),Vzep.begin(),Vzep.end());
	lambda_mac.insert(lambda_mac.end(),lambda_macp.begin(),lambda_macp.end());
	//ASSIGNING ALL THIS TO OUT INPUT ARRAY
	//xe=xei;
	//ye=yei;
	//Vxe=Vxei;
	//Vye=Vyei;
	//Vze=Vzei;
	//lambda_mac=lambda_maci;
	//printf("Left %d\n",xe.size());
	};
/*
//ANOTHER ADDITIONAL PARAMETER IS A PRECALCULATED ARRAY cos_sin TO PRODUCE ELECTRONS WITH cos_thetha ANGULAR DISTRIBUTION HOWEVER ALL THE PARTICLE DATA IS INSIDE ONE VECTOR OF VECTORS
void wall_effects(vector < vector <double> > &part_param_e, 
double Rx, double Ry, 
double dmax, double Emax, double refl0, double Erefl, double sigma_v, 
vector <double> cos_sin, double dt, int geoflag, gsl_rng *rr)
	{
	int SIZE=xe.size();
	int cs_size=cos_sin.size()/2;
	//NORMALIZATION
	Erefl*=2.0*CTOM;
	Emax*=2.0*CTOM;
	//VECTORS TO STORE PARTICLES THAT ARE OUTSIDE THE PIPE WITH ADDITIONAL PARAMETER ENERGY, NUMBER OF EMITTED ELECTRONS, REFLECTION COEFFICIENT, CROSS POINTS AND NORMALS TO THE SURFACE, dt2 IS A PART OF TIME STEP FROM SURFACE TO THE PARTICLE INSIDE THE SURFACE / EACH ARTICLE HAS A CHANCE TO BE REMOVED
	vector <double> energyo, xc, yc, ex,ey;
	vector < vector <double> > part_param_o,
	//REFLECTION STATUS ABSORBED OR REFLECTED
	vector <int> reflo, nemito;
	//VECTORS TO STORE PARTICLES THAT WERE PRODUCED
	//vector <double> xep,yep,Vxep,Vyep,Vzep,lambda_macp;
	vector < vector <double> > part_param_p,
	//printf("First loop\n");
	//WE CONSTRUCT VECTORS OF PARTICLE OUTSIDE THE GRID AND INSIDE THE GRID
	for(int i=0;i<SIZE;++i)
		{
		//WHILE LOOP IS USED FOR THE CASE OF A CORNER (PROBABLY VERY RARE) WHEN MULTIPLE REFLECTION CAN HAPPEN
		if(outside(geoflag, Rx, Ry, part_param_e[i][0], part_param_e[i][1]))
			{
			//PARTICLES OUTSIDE
			part_param_o.push_back(part_param_e[i]);
			//IT IS ASSUMED THAT ENERGY CONSTANTS ARE NORMALIZED
			energyo.push_back(pow(part_param_e[i][2],2)+pow(part_param_e[i][3],2)+pow(part_param_e[i][4],2));
			//REPLACE PARTICLES
			part_param_e[i]=part_param_e.back();
			//ERASING LAST ELEMENTS
			part_param_e.pop_back();
			SIZE--;
			i--;
			//printf("%e ", energyo.back()*EM/2/EC);
			//CALCULATE IF AN ELECTRON IS REFLECTED
			if(refl(refl0, Erefl, energyo.back())>gsl_rng_uniform(rr))
				{
				reflo.push_back(1);
				//printf("zaebok \n");
				}
			else
				{
				reflo.push_back(0);
				}
			//NUMBER OF ELECTRONS TO EMIT
			nemito.push_back(n_emitted(sey(Emax, dmax, 1.2, energyo.back()),rr));
			//nemito.push_back(n_emitted(sey_const_max(Emax, dmax, 1.2, energyo.back()),rr));
			//CROSS POINTS AND NORM VECTORS
			xc.push_back(0.0);
			yc.push_back(0.0);
			ex.push_back(0.0);
			ey.push_back(0.0);
			//ASSIGNING CROSSING POINTS AND NORMAL VECTORS
			cross_point(geoflag, Rx, Ry, xe[i], ye[i], Vxe[i], Vye[i], dt, xc.back(), yc.back(), ex.back(), ey.back());
			//ADD SEARCH FOR LEFT OVER
			}
		}
	int SIZEE=xeo.size();
	//NOW WE CONSTRUCT A NEW "OUTPUT" VECTOR OF PARTICLES CONSISTING OF INTERNAL PARTICLES, REFLECTED PARTICLES AND SECONDARIES
	//1. CALCULATE SECONDARIES
	for(int i=0;i<SIZEE;++i)
		{
		//IF THE PARTICLE IS REFLECTED, THEN IT IS ADDED TO THE NEWBORN VECTOR WITH REFLECTED VELOCITIES AND COORDINATES (NOT TAKING INTO ACCOUNT THE EXTERNAL FORCE)
		if(reflo[i]==1)
			{
			//REFLECT AN ELECTRON AND ATTACHING IT TO VECTOR OF SECONDARY PARTICLES
			reflect(xep,yep,Vxep,Vyep,Vzep,lambda_macp,xeo[i], yeo[i], Vxeo[i], Vyeo[i], lambda_maco[i],xc[i], yc[i],ex[i], ey[i], dt, Rx, Ry, geoflag);
			}
		//PRODUCTION OF TRUE SECONDARIES
		for(int j=0;j<nemito[i];++j)
			{
			//CHOSING RANDOMLY THE ANGLE OF EMISSION
			int angle_number=gsl_rng_uniform_int(rr,cs_size);
			//PERPENDICULAR COMPONENT OF THE VELOCITY
			double ep=cos_sin.at(angle_number);
			//PARALLEL TO SURFACE COMPONENT OF VELOCITY
			double eh=cos_sin.at(angle_number+cs_size);
			//GENERATING VELOCITY WITH THE POISSON DISTRIBUTION
			double g_velocity=gen_velocity(sigma_v, rr);
			//HERE WE MAKE A HYBRID MODEL ELECTRONS ARE PRODUCED UNCORRELATED, BUT THE ENERGY OF EACH
			//ELECTRON CAN NOT EXCEED THE IMPACT ENERGY (IN TOTAL STILL POSSIBLE TO EXCEED)
			int checker=0;
			while(g_velocity*g_velocity>Vxeo[i]*Vxeo[i]+Vyeo[i]*Vyeo[i]+Vzeo[i]*Vzeo[i])
				{
				checker++;
				if(checker>100 && checker%100==0)
					{
					printf("Alert: More than 100 iterations %e %e\n",g_velocity*g_velocity,Vxeo[i]*Vxeo[i]+Vyeo[i]*Vyeo[i]+Vzeo[i]*Vzeo[i]);
					}
				g_velocity=gen_velocity(sigma_v, rr);
				}
			double vgp=g_velocity*ep;
			double vgh=g_velocity*eh;
			Vxep.push_back(vgp*ex[i]+vgh*ey[i]);
			Vyep.push_back(vgp*ey[i]-vgh*ex[i]);
			//THIS MOMENT SHOULD BE WORKED OUT IN DETAILS TO INCLUDE THE EFFECT OF LOCAL FIELD
			xep.push_back(xc[i]*0.99999);
			yep.push_back(yc[i]*0.99999);
			Vzep.push_back(0.0);
			lambda_macp.push_back(lambda_maco[i]);
			if(fabs(Vxep.back())>1e10 || fabs(Vyep.back())>1e10 || fabs(xep.back())>Rx || fabs(yep.back())>Ry)
				{
				printf("%e %e %e %e %e %e\n",xep.back(),yep.back(),ex[i],ey[i],xc[i],yc[i]);
				sleep(0.5);
				}
			}
		}
	//MERGING ARRAYS OF INTERNAL AND NEWLY PRODUCED ELECTRONS
	xe.insert(xei.end(),xep.begin(),xep.end());
	ye.insert(yei.end(),yep.begin(),yep.end());
	Vxe.insert(Vxei.end(),Vxep.begin(),Vxep.end());
	Vye.insert(Vyei.end(),Vyep.begin(),Vyep.end());
	Vze.insert(Vzei.end(),Vzep.begin(),Vzep.end());
	lambda_mac.insert(lambda_maci.end(),lambda_macp.begin(),lambda_macp.end());
	//ASSIGNING ALL THIS TO OUT INPUT ARRAY
	//xe=xei;
	//ye=yei;
	//Vxe=Vxei;
	//Vye=Vyei;
	//Vze=Vzei;
	//lambda_mac=lambda_maci;
	//printf("Left %d\n",xe.size());
	};
*/

//FUNCTION THAT REFLECTS AN ELECTRON WITH i NUMBER
//UP TO NOW USES SIMPLIFIED ALGORITHMS
void reflect(
//VECTORS TO ADD IN
vector <double> &xep, vector <double> &yep, vector <double> &Vxep, vector <double> &Vyep, vector <double> &Vzep, vector <double> &lambda_macp,
//VECTORS TO TAKE FROM
double xeo, double yeo, double Vxeo, double Vyeo, double lambda_maco,
//CROSSING POINT
double xc, double yc, double ex, double ey, double dt,
//PIPE PARAMETERS
double Rx, double Ry, int geoflag)
	{
	//REFLECTION FOR RACTANGULAR GEOMETRY
	if(geoflag==2)
			{
			//THE EFFECT OF THE FIELD S NEGLECTED AT THIS MOMENT
			double xabs=fabs(xeo);
			double yabs=fabs(yeo);
			//THE CASE WHEN BOTH MAXIMUM PIPE DIMENSIONS ARE EXCEEDED
			if(xabs>Rx && yabs>Ry)
				{
				double dx=xabs-Rx;
				double dy=yabs-Ry;
				xep.push_back(xeo*(Rx-dx)/xabs);
				Vxep.push_back(-Vxeo);
				yep.push_back(yeo*(Ry-dy)/yabs);
				Vyep.push_back(-Vyeo);
				Vzep.push_back(0.0);
				lambda_macp.push_back(lambda_maco);
				}
			//THE CASE WHEN HORIZONTAL PIPE SIZE IS EXCEEDED
			else if(xabs>Rx)
				{
				double dx=xabs-Rx;
				xep.push_back(xeo*(Rx-dx)/xabs);
				Vxep.push_back(-Vxeo);
				yep.push_back(yeo);
				Vyep.push_back(Vyeo);
				Vzep.push_back(0.0);
				lambda_macp.push_back(lambda_maco);
				//printf("%e\n",xep.back());
				}
			//THE CASE WHEN VERTICAL PIPE SIZE IS EXCEEDED
			else if(yabs>Ry)
				{
				double dy=yabs-Ry;
				yep.push_back(yeo*(Ry-dy)/yabs);
				Vyep.push_back(-Vyeo);
				xep.push_back(xeo);
				Vxep.push_back(Vxeo);
				Vzep.push_back(0.0);
				lambda_macp.push_back(lambda_maco);
				}
			}
	//REFLECTION FOR ROUND GEOMETRY
	if(geoflag==1)
		{
		//HOW FAR IS THE ELECTRON OUTSIDE
		double dx=xeo-xc;
		double dy=yeo-yc;
		double prevx=xeo-Vxeo*dt;
		double prevy=yeo-Vyeo*dt;
		//printf("dx=%e dy=%e\n",dx,dy);
		//SCALAR MULTIPLICATION WITH NORM
		double scalar=dx*ex+dy*ey;
		double scalarV=Vxeo*ex+Vyeo*ey;
		xep.push_back(xeo+2.0*ex*scalar);
		yep.push_back(yeo+2.0*ey*scalar);
		Vxep.push_back(Vxeo+2.0*ex*scalarV);
		Vyep.push_back(Vyeo+2.0*ey*scalarV);
		Vzep.push_back(0.0);
		lambda_macp.push_back(lambda_maco);
		if(pow(xep.back(),2)+pow(yep.back(),2)>Rx*Rx)
			{
			printf("Wrong reflection: xeo=%e yeo=%e prevx=%e prevy=%e xc=%e yc=%e -> xe=%e ye=%e\n",xeo,yeo,prevx,prevy,xc,yc,xep.back(),yep.back());
			}
		}
	};

void reflect(
//VECTORS TO ADD IN
vector <electron> &ELECTRONS,
//VECTORS TO TAKE FROM
double xeo, double yeo, double Vxeo, double Vyeo, double Vzeo, double lambda_maco,
//CROSSING POINT
double xc, double yc, double ex, double ey, double dt,
//PIPE PARAMETERS
double Rx, double Ry, int geoflag)
	{
	//REFLECTION FOR RACTANGULAR GEOMETRY
	if(geoflag==2)
			{
			//THE EFFECT OF THE FIELD S NEGLECTED AT THIS MOMENT
			double xabs=fabs(xeo);
			double yabs=fabs(yeo);
			electron temp;
			//THE CASE WHEN BOTH MAXIMUM PIPE DIMENSIONS ARE EXCEEDED
			if(xabs>Rx && yabs>Ry)
				{
				double dx=xabs-Rx;
				double dy=yabs-Ry;
				temp.x=xeo*(Rx-dx)/xabs;
				temp.y=yeo*(Ry-dy)/yabs;
				temp.Vx=-Vxeo;
				temp.Vy=-Vyeo;
				temp.Vz=Vzeo;
				temp.lambda_mac=lambda_maco;
				ELECTRONS.push_back(temp);
				}
			//THE CASE WHEN HORIZONTAL PIPE SIZE IS EXCEEDED
			else if(xabs>Rx)
				{
				double dx=xabs-Rx;
				//xep.push_back(xeo*(Rx-dx)/xabs);
				//Vxep.push_back(-Vxeo);
				//yep.push_back(yeo);
				//Vyep.push_back(Vyeo);
				//Vzep.push_back(0.0);
				//lambda_macp.push_back(lambda_maco);
				temp.x=xeo*(Rx-dx)/xabs;
				temp.y=yeo;
				temp.Vx=-Vxeo;
				temp.Vy=Vyeo;
				temp.Vz=Vzeo;
				temp.lambda_mac=lambda_maco;
				ELECTRONS.push_back(temp);
				//printf("%e\n",xep.back());
				}
			//THE CASE WHEN VERTICAL PIPE SIZE IS EXCEEDED
			else if(yabs>Ry)
				{
				double dy=yabs-Ry;
				//yep.push_back(yeo*(Ry-dy)/yabs);
				//Vyep.push_back(-Vyeo);
				//xep.push_back(xeo);
				//Vxep.push_back(Vxeo);
				//Vzep.push_back(0.0);
				//lambda_macp.push_back(lambda_maco);
				temp.x=xeo;
				temp.y=yeo*(Ry-dy)/yabs;
				temp.Vx=Vxeo;
				temp.Vy=-Vyeo;
				temp.Vz=Vzeo;
				temp.lambda_mac=lambda_maco;
				ELECTRONS.push_back(temp);
				}
			}
	//REFLECTION FOR ROUND GEOMETRY
	if(geoflag==1)
		{
		electron temp;
		//HOW FAR IS THE ELECTRON OUTSIDE
		double dx=xeo-xc;
		double dy=yeo-yc;
		double prevx=xeo-Vxeo*dt;
		double prevy=yeo-Vyeo*dt;
		//printf("dx=%e dy=%e\n",dx,dy);
		//SCALAR MULTIPLICATION WITH NORM
		double scalar=dx*ex+dy*ey;
		double scalarV=Vxeo*ex+Vyeo*ey;
		//xep.push_back(xeo-2.0*ex*scalar);
		//yep.push_back(yeo-2.0*ey*scalar);
		//Vxep.push_back(Vxeo-2.0*ex*scalarV);
		//Vyep.push_back(Vyeo-2.0*ey*scalarV);
		//Vzep.push_back(0.0);
		//lambda_macp.push_back(lambda_maco);
		//REPLACE XEO AND YEO WITH XC AND YC
		temp.x=xc*0.99999-1.0*ex*scalar;
		temp.y=yc*0.99999-1.0*ey*scalar;
		temp.Vx=Vxeo-2.0*ex*scalarV;
		temp.Vy=Vyeo-2.0*ey*scalarV;
		temp.Vz=Vzeo;
		temp.lambda_mac=lambda_maco;
		ELECTRONS.push_back(temp);
		double tempconst=pow(ELECTRONS.back().x/Rx,2)+pow(ELECTRONS.back().y/Ry,2);
		if(tempconst>1.0)
			{
			printf("tempconst %e\n",tempconst);
			printf("summ of normal squares ex*ex+ey*ey=%e\n",ex*ex+ey*ey);
			printf("summ of normalized cross coordinates pow(%e/Rx,2)+pow(%e/Ry,2)=%e\n",xc,yc,pow(xc/Rx,2)+pow(yc/Ry,2));
			printf("coordinates xeo=%e yeo=%e and pow(xeo/Rx,2)+pow(yeo/Ry,2)=%e\n",xeo,yeo,pow(xeo/Rx,2)+pow(yeo/Ry,2));
			tempconst=sqrt(tempconst);
			ELECTRONS.back().x/=tempconst;
			ELECTRONS.back().y/=tempconst;
			//printf("Wrong reflection: xeo=%e yeo=%e prevx=%e prevy=%e xc=%e yc=%e -> xe=%e ye=%e\n",xeo,yeo,prevx,prevy,xc,yc,ELECTRONS.back().x,ELECTRONS.back().y);
			}
		}
	};












//KICK ALL THE ELECTRONS IN THE CLOUD (ONLY CHANGES THE VELOCITY BY adt)
void accelerate_all(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, double *Ex2D, double *Ey2D, int Nx, int Ny, double dx, double dy)
	{
	//CALCULATING NUMBER OF ELECTRON MACROPARTICLES
	int SIZE=xe.size();

	{

	for(int i=0;i<SIZE;++i)
		{
		//CALCULATING THE ACCELERATION ACRING ON THE ELECTRON
		double axdt=get_value(Ex2D,Nx,Ny,dx,dy,xe[i],ye[i]);
		double aydt=get_value(Ey2D,Nx,Ny,dx,dy,xe[i],ye[i]);
		//ACCELERATING THE ELECTRON
		accelerate_electron(axdt,aydt,Vxe[i],Vye[i]);
		}
	}
	}








//SHIFT ALL ELECTRONS (ONLY SHIFTS ELECTRONS BY Vdt)
void shift_all(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, double dt)
	{
	//CALCULATING NUMBER OF ELECTRON MACROPARTICLES
	int SIZE=xe.size();

	{

	for(int i=0;i<SIZE;++i)
		{
		shift_electron(xe[i], ye[i], Vxe[i], Vye[i], dt);
		}
	}
	}








//MOVE ALL ELECTRONS (CHANGES VELOCITY AND COORDINATE SIMULATNEOUSLY)
void move_all(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, double *Ex2D, double *Ey2D, int Nx, int Ny, double dx, double dy, double dt)
	{
	//CALCULATING NUMBER OF ELECTRON MACROPARTICLES
	int SIZE=xe.size();

	{

	for(int i=0;i<SIZE;++i)
		{
		//SHIFTING THE ELECTRON
		shift_electron(xe[i], ye[i], Vxe[i], Vye[i], dt);
		//CALCULATING THE ACCELERATION ACRING ON THE ELECTRON
		double axdt=get_value(Ex2D,Nx,Ny,dx,dy,xe[i],ye[i]);
		double aydt=get_value(Ey2D,Nx,Ny,dx,dy,xe[i],ye[i]);
		//ACCELERATING THE ELECTRON
		accelerate_electron(axdt,aydt,Vxe[i],Vye[i]);
		}
	}
	}









//MOVE ALL ELECTRONS IN MAGNETIC FIELD (CHANGES VELOCITY AND COORDINATE SIMULATNEOUSLY)
void move_all_in_B(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, 
double *Ex2D, double *Ey2D, 
double cosi, double sini, 
double Rx, double Ry, 
int Nx, int Ny, 
double dx, double dy, double dt, 
int geoflag) //geoflag 0 - rect; 1 - circle
	{
	//CALCULATING NUMBER OF ELECTRON MACROPARTICLES
	int SIZE=xe.size();
	double Dx2=pow(dx*(Nx-1)/2.0,2);
	double Dy2=pow(dy*(Nx-1)/2.0,2);
	double Rx2=Rx*Rx;
	for(int i=0;i<SIZE;++i)
		{
		//SHIFTING THE ELECTRON
		shift_electron(xe[i], ye[i], Vxe[i], Vye[i], dt);
		if((geoflag==2 && pow(xe[i],2)<Dx2 && pow(ye[i],2)<Dy2 ) || (geoflag==1 && xe[i]*xe[i]+ye[i]*ye[i]<Rx*Rx ) )
			{
			//CALCULATING THE ACCELERATION ACTING ON THE ELECTRON
			double axdt=get_value(Ex2D,Nx,Ny,dx,dy,xe[i],ye[i]);
			double aydt=get_value(Ey2D,Nx,Ny,dx,dy,xe[i],ye[i]);
			//ACCELERATING THE ELECTRON
			accelerate_electron(axdt/2.0,aydt/2.0,Vxe[i],Vye[i]);
			//ROTATE
			rotate_velocity(Vxe[i],Vze[i],cosi,sini);
			//ACCELERATING
			accelerate_electron(axdt/2.0,aydt/2.0,Vxe[i],Vye[i]);
			}
		
		}
	
	}
//MOVE ALL ELECTRONS IN MAGNETIC FIELD (CHANGES VELOCITY AND COORDINATE SIMULATNEOUSLY) USING ELECTRON STRUCT
void move_all_in_B(vector <electron> &ELECTRONS, 
double *Ex2D, double *Ey2D, 
double cosi, double sini, 
double Rx, double Ry, 
int Nx, int Ny, 
double dx, double dy, double dt, 
int geoflag) //geoflag 0 - rect; 1 - circle
	{
	//CALCULATING NUMBER OF ELECTRON MACROPARTICLES
	int SIZE=ELECTRONS.size();
	double Dx2=pow(dx*(Nx-1)/2.0,2);
	double Dy2=pow(dy*(Nx-1)/2.0,2);
	double Rx2=Rx*Rx;
	for(int i=0;i<SIZE;++i)
		{
		//SHIFTING THE ELECTRON
		shift_electron(ELECTRONS[i].x, ELECTRONS[i].y, ELECTRONS[i].Vx, ELECTRONS[i].Vy, dt);
		if((geoflag==1 && pow(ELECTRONS[i].x,2)<Dx2 && pow(ELECTRONS[i].y,2)<Dy2 ) || (geoflag==2 && pow(ELECTRONS[i].x/Rx,2)+pow(ELECTRONS[i].y/Ry,2)<1.0 ) )
			{
			//CALCULATING THE ACCELERATION ACTING ON THE ELECTRON
			double axdt=get_value(Ex2D,Nx,Ny,dx,dy,ELECTRONS[i].x,ELECTRONS[i].y);
			double aydt=get_value(Ey2D,Nx,Ny,dx,dy,ELECTRONS[i].x,ELECTRONS[i].y);
			//ACCELERATING THE ELECTRON
			accelerate_electron(axdt/2.0,aydt/2.0,ELECTRONS[i].Vx,ELECTRONS[i].Vy);
			//ROTATE
			rotate_velocity(ELECTRONS[i].Vx,ELECTRONS[i].Vz,cosi,sini);
			//ACCELERATING
			accelerate_electron(axdt/2.0,aydt/2.0,ELECTRONS[i].Vx,ELECTRONS[i].Vy);
			}
		
		}
	
	}





//MOVE ALL ELECTRONS IN MAGNETIC FIELD (CHANGES VELOCITY AND COORDINATE SIMULATNEOUSLY) USING THE PARTICLE WEIGHTS FROM INTERPOLATION STEP
void move_all_in_B(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, 
double *Ex2D, double *Ey2D,
//PARAMETERS OF A PARTICLE RELATIVE TO THE GRID
vector <double> &downleft,vector <double> &upleft,vector <double> &downright,vector <double> &upright,vector <int> &low_x,vector <int> &low_y,
double cosi, double sini, 
double Rx, double Ry, 
int Nx, int Ny, 
double dx, double dy, double dt, 
int geoflag) //geoflag 2 - rect 1 - circle
	{
	//CALCULATING NUMBER OF ELECTRON MACROPARTICLES
	int SIZE=xe.size();
	double Dx2=pow(dx*(Nx-1)/2.0,2);
	double Dy2=pow(dx*(Ny-1)/2.0,2);
	double Rx2=Rx*Rx;
	for(int i=0;i<SIZE;++i)
		{
		//SHIFTING THE ELECTRON
		shift_electron(xe[i], ye[i], Vxe[i], Vye[i], dt);
		if((geoflag==2 && pow(xe[i],2)<Dx2 && pow(ye[i],2)<Dy2 ) || (geoflag==1 && xe[i]*xe[i]+ye[i]*ye[i]<Rx*Rx ) )
			{
			//CALCULATING THE ACCELERATION ACTING ON THE ELECTRON
			double axdt=get_value(Ex2D, Nx, Ny, dx, dy, downleft[i], upleft[i], downright[i], upright[i],  low_x[i],  low_y[i]);//get_value(Ex2D,Nx,Ny,dx,dy,xe[i],ye[i]);
			double aydt=get_value(Ey2D, Nx, Ny, dx, dy, downleft[i], upleft[i], downright[i], upright[i],  low_x[i],  low_y[i]);
			//ACCELERATING THE ELECTRON
			accelerate_electron(axdt/2.0,aydt/2.0,Vxe[i],Vye[i]);
			//ROTATE
			rotate_velocity(Vxe[i],Vze[i],cosi,sini);
			//ACCELERATING
			accelerate_electron(axdt/2.0,aydt/2.0,Vxe[i],Vye[i]);
			}
		
		}
	
	}
//MOVE ALL ELECTRONS IN MAGNETIC FIELD (CHANGES VELOCITY AND COORDINATE SIMULATNEOUSLY) USING THE PARTICLE WEIGHTS FROM INTERPOLATION STEP USING ELECTRON STRUCTURE
/*void move_all_in_B(vector <electron> &ELECTRONS, 
double *Ex2D, double *Ey2D,
//PARAMETERS OF A PARTICLE RELATIVE TO THE GRID
vector <weight> WEIGHTS,
double cosi, double sini, 
double Rx, double Ry, 
int Nx, int Ny, 
double dx, double dy, double dt, 
int geoflag) //geoflag 0 - rect 1 - circle
	{
	//CALCULATING NUMBER OF ELECTRON MACROPARTICLES
	int SIZE=ELECTRONS.size();
	double Dx2=pow(dx*(Nx-1)/2.0,2);
	double Dy2=pow(dx*(Ny-1)/2.0,2);
	double Rx2=Rx*Rx;
	vector <int> indices;
	for(int i=0;i<SIZE;++i)
		{
		//SHIFTING THE ELECTRON
		shift_electron(ELECTRONS[i].x, ELECTRONS[i].y, ELECTRONS[i].Vx, ELECTRONS[i].Vy, dt);
		if(
(geoflag==0 && pow(ELECTRONS[i].x,2)<Dx2 && pow(ELECTRONS[i].y,2)<Dy2 ) || 
(geoflag==1 && (pow(ELECTRONS[i].x/Rx,2)+pow(ELECTRONS[i].y/Ry,2))<1.0 ) 
)
			{
			indices.push_back(i);
			}
		}
	int SIZEIN=indices.size();
	int SIZE4=SIZEIN/4;
	int residue=SIZEIN%4;
	for(int i=0;i<SIZE4;i++)
		{
		int BEGIN=indices[4*i];
		int BEGIN1=indices[4*i+1];
		int BEGIN2=indices[4*i+2];
		int BEGIN3=indices[4*i+3];
		//CALCULATING THE ACCELERATION ACTING ON THE ELECTRON 0
		double axdt=get_value(Ex2D, Nx, Ny, dx, dy, WEIGHTS.at(BEGIN).downleft, WEIGHTS.at(BEGIN).upleft, WEIGHTS.at(BEGIN).downright, WEIGHTS.at(BEGIN).upright,  WEIGHTS.at(BEGIN).low_x,  WEIGHTS.at(BEGIN).low_y)/2.0;
		double aydt=get_value(Ey2D, Nx, Ny, dx, dy, WEIGHTS.at(BEGIN).downleft, WEIGHTS.at(BEGIN).upleft, WEIGHTS.at(BEGIN).downright, WEIGHTS.at(BEGIN).upright,  WEIGHTS.at(BEGIN).low_x,  WEIGHTS.at(BEGIN).low_y)/2.0;
		//ACCELERATING THE ELECTRON
		accelerate_electron(axdt,aydt,ELECTRONS.at(BEGIN).Vx,ELECTRONS.at(BEGIN).Vy);
		//ROTATE
		rotate_velocity(ELECTRONS.at(BEGIN).Vx,ELECTRONS.at(BEGIN).Vz,cosi,sini);
		//ACCELERATING
		accelerate_electron(axdt,aydt,ELECTRONS.at(BEGIN).Vx,ELECTRONS.at(BEGIN).Vy);

		//CALCULATING THE ACCELERATION ACTING ON THE ELECTRON 1 
		axdt=get_value(Ex2D, Nx, Ny, dx, dy, WEIGHTS.at(BEGIN1).downleft, WEIGHTS.at(BEGIN1).upleft, WEIGHTS.at(BEGIN1).downright, WEIGHTS.at(BEGIN1).upright,  WEIGHTS.at(BEGIN1).low_x,  WEIGHTS.at(BEGIN1).low_y)/2.0;
		aydt=get_value(Ey2D, Nx, Ny, dx, dy, WEIGHTS.at(BEGIN1).downleft, WEIGHTS.at(BEGIN1).upleft, WEIGHTS.at(BEGIN1).downright, WEIGHTS.at(BEGIN1).upright,  WEIGHTS.at(BEGIN1).low_x,  WEIGHTS.at(BEGIN1).low_y)/2.0;
		//ACCELERATING THE ELECTRON
		accelerate_electron(axdt,aydt,ELECTRONS.at(BEGIN1).Vx,ELECTRONS.at(BEGIN1).Vy);
		//ROTATE
		rotate_velocity(ELECTRONS.at(BEGIN1).Vx,ELECTRONS.at(BEGIN1).Vz,cosi,sini);
		//ACCELERATING
		accelerate_electron(axdt,aydt,ELECTRONS.at(BEGIN1).Vx,ELECTRONS.at(BEGIN1).Vy);

		//CALCULATING THE ACCELERATION ACTING ON THE ELECTRON 2
		axdt=get_value(Ex2D, Nx, Ny, dx, dy, WEIGHTS.at(BEGIN2).downleft, WEIGHTS.at(BEGIN2).upleft, WEIGHTS.at(BEGIN2).downright, WEIGHTS.at(BEGIN2).upright,  WEIGHTS.at(BEGIN2).low_x,  WEIGHTS.at(BEGIN2).low_y)/2.0;
		aydt=get_value(Ey2D, Nx, Ny, dx, dy, WEIGHTS.at(BEGIN2).downleft, WEIGHTS.at(BEGIN2).upleft, WEIGHTS.at(BEGIN2).downright, WEIGHTS.at(BEGIN2).upright,  WEIGHTS.at(BEGIN2).low_x,  WEIGHTS.at(BEGIN2).low_y)/2.0;
		//ACCELERATING THE ELECTRON
		accelerate_electron(axdt,aydt,ELECTRONS.at(BEGIN2).Vx,ELECTRONS.at(BEGIN2).Vy);
		//ROTATE
		rotate_velocity(ELECTRONS.at(BEGIN2).Vx,ELECTRONS.at(BEGIN2).Vz,cosi,sini);
		//ACCELERATING
		accelerate_electron(axdt,aydt,ELECTRONS.at(BEGIN2).Vx,ELECTRONS.at(BEGIN2).Vy);

		//CALCULATING THE ACCELERATION ACTING ON THE ELECTRON 3
		axdt=get_value(Ex2D, Nx, Ny, dx, dy, WEIGHTS.at(BEGIN3).downleft, WEIGHTS.at(BEGIN3).upleft, WEIGHTS.at(BEGIN3).downright, WEIGHTS.at(BEGIN3).upright,  WEIGHTS.at(BEGIN3).low_x,  WEIGHTS.at(BEGIN3).low_y)/2.0;
		aydt=get_value(Ey2D, Nx, Ny, dx, dy, WEIGHTS.at(BEGIN3).downleft, WEIGHTS.at(BEGIN3).upleft, WEIGHTS.at(BEGIN3).downright, WEIGHTS.at(BEGIN3).upright,  WEIGHTS.at(BEGIN3).low_x,  WEIGHTS.at(BEGIN3).low_y)/2.0;
		//ACCELERATING THE ELECTRON
		accelerate_electron(axdt,aydt,ELECTRONS.at(BEGIN3).Vx,ELECTRONS.at(BEGIN3).Vy);
		//ROTATE
		rotate_velocity(ELECTRONS.at(BEGIN3).Vx,ELECTRONS.at(BEGIN3).Vz,cosi,sini);
		//ACCELERATING
		accelerate_electron(axdt,aydt,ELECTRONS.at(BEGIN3).Vx,ELECTRONS.at(BEGIN3).Vy);
		}
	for(int i=4*SIZE4;i<4*SIZE4+residue;i++)
		{
		int ii=indices[i];
			double axdt=get_value(Ex2D, Nx, Ny, dx, dy, WEIGHTS.at(ii).downleft, WEIGHTS.at(ii).upleft, WEIGHTS.at(ii).downright, WEIGHTS.at(ii).upright,  WEIGHTS.at(ii).low_x,  WEIGHTS.at(ii).low_y)/2.0;
			double aydt=get_value(Ey2D, Nx, Ny, dx, dy, WEIGHTS.at(ii).downleft, WEIGHTS.at(ii).upleft, WEIGHTS.at(ii).downright, WEIGHTS.at(ii).upright,  WEIGHTS.at(ii).low_x,  WEIGHTS[i].low_y)/2.0;
			//ACCELERATING THE ELECTRON
			accelerate_electron(axdt,aydt,ELECTRONS.at(ii).Vx,ELECTRONS.at(ii).Vy);
			//ROTATE
			rotate_velocity(ELECTRONS.at(ii).Vx,ELECTRONS.at(ii).Vz,cosi,sini);
			//ACCELERATING
			accelerate_electron(axdt,aydt,ELECTRONS.at(ii).Vx,ELECTRONS.at(ii).Vy);
		}
	
	
	}*/


/*void move_all_in_B(vector <electron> &ELECTRONS, 
double *Ex2D, double *Ey2D,
//PARAMETERS OF A PARTICLE RELATIVE TO THE GRID
vector <weight> WEIGHTS,
double cosi, double sini, 
double Rx, double Ry, 
int Nx, int Ny, 
double dx, double dy, double dt, 
int geoflag) //geoflag 0 - rect 1 - circle
	{
	//CALCULATING NUMBER OF ELECTRON MACROPARTICLES
	int SIZE=ELECTRONS.size();
	double Dx2=pow(dx*(Nx-1)/2.0,2);
	double Dy2=pow(dy*(Ny-1)/2.0,2);
	for(int i=0;i<SIZE;++i)
		{
		//SHIFTING THE ELECTRON
		//shift_electron(ELECTRONS[i].x, ELECTRONS[i].y, ELECTRONS[i].Vx, ELECTRONS[i].Vy, dt);
		ELECTRONS[i].x+=ELECTRONS[i].Vx*dt;
		ELECTRONS[i].y+=ELECTRONS[i].Vy*dt;
		if((geoflag==1 && (pow(ELECTRONS[i].x/Rx,2)+pow(ELECTRONS[i].y/Ry,2))<1.0 ) ||
(geoflag==0 && pow(ELECTRONS[i].x,2)<Dx2 && pow(ELECTRONS[i].y,2)<Dy2 )
 
)			{

			double axdt=get_value(Ex2D, Nx, Ny, dx, dy, WEIGHTS[i].downleft, WEIGHTS[i].upleft, WEIGHTS[i].downright, WEIGHTS[i].upright,  WEIGHTS[i].low_x,  WEIGHTS[i].low_y)/2.0;
			double aydt=get_value(Ey2D, Nx, Ny, dx, dy, WEIGHTS[i].downleft, WEIGHTS[i].upleft, WEIGHTS[i].downright, WEIGHTS[i].upright,  WEIGHTS[i].low_x,  WEIGHTS[i].low_y)/2.0;
			//ACCELERATING THE ELECTRON
			//accelerate_electron(axdt,aydt,ELECTRONS[i].Vx,ELECTRONS[i].Vy);
			ELECTRONS[i].Vx+=axdt;
			ELECTRONS[i].Vy+=aydt;
			//ROTATE
			//rotate_velocity(ELECTRONS[i].Vx,ELECTRONS[i].Vz,cosi,sini);
			double tVx=ELECTRONS[i].Vx*cosi+ELECTRONS[i].Vz*sini;
			double tVz=-ELECTRONS[i].Vx*sini+ELECTRONS[i].Vz*cosi;
			ELECTRONS[i].Vx=tVx;
			ELECTRONS[i].Vz=tVz;
			//ACCELERATING
			//accelerate_electron(axdt,aydt,ELECTRONS[i].Vx,ELECTRONS[i].Vy);
			ELECTRONS[i].Vx+=axdt;
			ELECTRONS[i].Vy+=aydt;

				}

		}
	}*/
void move_all_in_B(vector <electron> &ELECTRONS,
double *Ex2D, double *Ey2D,
//PARAMETERS OF A PARTICLE RELATIVE TO THE GRID
vector <weight> WEIGHTS,
double cosi, double sini, 
double Rx, double Ry, 
int Nx, int Ny, 
double dx, double dy, double dt, 
int geoflag) //geoflag 0 - rect 1 - circle
	{
	//clock_t startm, stopm;
	//double shift_time=0.0;
	//double gettin_values_time=0.0;
	//double kick_rotate_time=0.0;
	//CALCULATING NUMBER OF ELECTRON MACROPARTICLES
	int SIZE=ELECTRONS.size();
	double Dx2=pow(dx*(Nx-1)/2.0,2);
	double Dy2=pow(dy*(Ny-1)/2.0,2);
	for(int i=0;i<SIZE;++i)
		{
		//SHIFTING THE ELECTRON
		//shift_electron(ELECTRONS[i].x, ELECTRONS[i].y, ELECTRONS[i].Vx, ELECTRONS[i].Vy, dt);
		//startm=clock();
		ELECTRONS[i].x+=ELECTRONS[i].Vx*dt;
		ELECTRONS[i].y+=ELECTRONS[i].Vy*dt;
		//stopm=clock();
		//shift_time+=((double)stopm-startm)/CLOCKS_PER_SEC;
		if((geoflag==1 && (ELECTRONS[i].x*ELECTRONS[i].x/Rx/Rx+ELECTRONS[i].y*ELECTRONS[i].y/Ry/Ry)<1.0 ) ||
(geoflag==2 && pow(ELECTRONS[i].x,2)<Dx2 && pow(ELECTRONS[i].y,2)<Dy2 )
 
)			{
			//printf("inside\n");
			//startm=clock();
			double axdt=get_value(Ex2D, Ny, WEIGHTS[i].downleft, WEIGHTS[i].upleft, WEIGHTS[i].downright, WEIGHTS[i].upright,  WEIGHTS[i].low_index)/2.0;
			double aydt=get_value(Ey2D, Ny, WEIGHTS[i].downleft, WEIGHTS[i].upleft, WEIGHTS[i].downright, WEIGHTS[i].upright,  WEIGHTS[i].low_index)/2.0;
			//stopm=clock();
			//gettin_values_time+=((double)stopm-startm)/CLOCKS_PER_SEC;
			//startm=clock();
			//ACCELERATING THE ELECTRON
			//accelerate_electron(axdt,aydt,ELECTRONS[i].Vx,ELECTRONS[i].Vy);
			ELECTRONS[i].Vx+=axdt;
			ELECTRONS[i].Vy+=aydt;
			//ROTATE
			//rotate_velocity(ELECTRONS[i].Vx,ELECTRONS[i].Vz,cosi,sini);
			double tVx=ELECTRONS[i].Vx*cosi+ELECTRONS[i].Vz*sini;
			double tVz=-ELECTRONS[i].Vx*sini+ELECTRONS[i].Vz*cosi;
			ELECTRONS[i].Vx=tVx;
			ELECTRONS[i].Vz=tVz;
			//ACCELERATING
			//accelerate_electron(axdt,aydt,ELECTRONS[i].Vx,ELECTRONS[i].Vy);
			ELECTRONS[i].Vx+=axdt;
			ELECTRONS[i].Vy+=aydt;
			//stopm=clock();
			//kick_rotate_time+=((double)stopm-startm)/CLOCKS_PER_SEC;
			}
		}
		//printf("shift_time=%e\n",shift_time*1000);
		//printf("gettin_values_time=%e\n",gettin_values_time*1000);
		//printf("kick_rotate_time=%e\n\n",kick_rotate_time*1000);
	}

/*void move_all_in_B(vector <electron> &ELECTRONS,
double *Ex2D, double *Ey2D,
//PARAMETERS OF A PARTICLE RELATIVE TO THE GRID
vector <weight> WEIGHTS,
double cosi, double sini, 
double Rx, double Ry, 
int Nx, int Ny, 
double dx, double dy, double dt, 
int geoflag) //geoflag 0 - rect 1 - circle
	{
	//clock_t startm, stopm;
	//double shift_time=0.0;
	//double gettin_values_time=0.0;
	//double kick_rotate_time=0.0;
	//CALCULATING NUMBER OF ELECTRON MACROPARTICLES
	int SIZE=ELECTRONS.size();
	double Dx2=pow(dx*(Nx-1)/2.0,2);
	double Dy2=pow(dy*(Ny-1)/2.0,2);
	for(int i=0;i<SIZE;++i)
		{
		//SHIFTING THE ELECTRON
		//shift_electron(ELECTRONS[i].x, ELECTRONS[i].y, ELECTRONS[i].Vx, ELECTRONS[i].Vy, dt);
		//startm=clock();
		ELECTRONS[i].x+=ELECTRONS[i].Vx*dt;
		ELECTRONS[i].y+=ELECTRONS[i].Vy*dt;
		//stopm=clock();
		//shift_time+=((double)stopm-startm)/CLOCKS_PER_SEC;
		if((geoflag==1 && (ELECTRONS[i].x*ELECTRONS[i].x/Rx/Rx+ELECTRONS[i].y*ELECTRONS[i].y/Ry/Ry)<1.0 ) ||
(geoflag==0 && pow(ELECTRONS[i].x,2)<Dx2 && pow(ELECTRONS[i].y,2)<Dy2 )
 
)			{
			//printf("inside\n");
			//startm=clock();
			double axdt=get_value(Ex2D, Ny, WEIGHTS[i].downleft, WEIGHTS[i].upleft, WEIGHTS[i].downright, WEIGHTS[i].upright,  WEIGHTS[i].low_index)/2.0;
			double aydt=get_value(Ey2D, Ny, WEIGHTS[i].downleft, WEIGHTS[i].upleft, WEIGHTS[i].downright, WEIGHTS[i].upright,  WEIGHTS[i].low_index)/2.0;
			//stopm=clock();
			//gettin_values_time+=((double)stopm-startm)/CLOCKS_PER_SEC;
			//startm=clock();
			//ACCELERATING THE ELECTRON
			//accelerate_electron(axdt,aydt,ELECTRONS[i].Vx,ELECTRONS[i].Vy);
			ELECTRONS[i].Vx+=axdt;
			ELECTRONS[i].Vy+=aydt;
			//ROTATE
			//rotate_velocity(ELECTRONS[i].Vx,ELECTRONS[i].Vz,cosi,sini);
			double tVx=ELECTRONS[i].Vx*cosi+ELECTRONS[i].Vz*sini;
			double tVz=-ELECTRONS[i].Vx*sini+ELECTRONS[i].Vz*cosi;
			ELECTRONS[i].Vx=tVx;
			ELECTRONS[i].Vz=tVz;
			//ACCELERATING
			//accelerate_electron(axdt,aydt,ELECTRONS[i].Vx,ELECTRONS[i].Vy);
			ELECTRONS[i].Vx+=axdt;
			ELECTRONS[i].Vy+=aydt;
			//stopm=clock();
			//kick_rotate_time+=((double)stopm-startm)/CLOCKS_PER_SEC;
			}
		}
		//printf("shift_time=%e\n",shift_time*1000);
		//printf("gettin_values_time=%e\n",gettin_values_time*1000);
		//printf("kick_rotate_time=%e\n\n",kick_rotate_time*1000);
	}*/


//MOVE ALL ELECTRONS IN MAGNETIC FIELD (CHANGES VELOCITY AND COORDINATE SIMULATNEOUSLY)
void move_all_in_rigid_B(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, double *Ey2D, double Rx, double Ry, int Nx, int Ny, double dx, double dy, double dt, int geoflag)
	{
	//CALCULATING NUMBER OF ELECTRON MACROPARTICLES
	int SIZE=ye.size();
	for(int i=0;i<SIZE;++i)
		{
		Vxe[i]=1e-6;
		shift_electron(xe[i], ye[i], Vxe[i], Vye[i], dt);
		if ((geoflag==2 && fabs(xe[i])<dx*(Nx-1)/2.0 && fabs(ye[i])<dy*(Ny-1)/2.0) || (geoflag==1 && pow(xe[i],2)+pow(ye[i],2)<Rx*Rx) )
			{
			//CALCULATING THE ACCELERATION ACTING ON THE ELECTRON
			double axdt=0.0;
			double aydt=get_value(Ey2D,Nx,Ny,dx,dy,xe[i],ye[i]);
			//ACCELERATING THE ELECTRON
			accelerate_electron(axdt/2.0,aydt/2.0,Vxe[i],Vye[i]);
			//ROTATE
			//rotate_velocity(Vxe[i],Vze[i],cosi,sini);
			//ACCELERATING
			accelerate_electron(axdt/2.0,aydt/2.0,Vxe[i],Vye[i]);
			}

		}
	
	}


//MOVE ALL ELECTRONS IN MAGNETIC FIELD (CHANGES VELOCITY AND COORDINATE SIMULATNEOUSLY) USING STRUCTURES FOR ELECTRONS // ELECTROBNS ARE STRUCTURES
void move_all_in_rigid_B(vector <electron> &ELECTRONS, 
double *Ey2D,
vector <weight> WEIGHTS,
double Rx, double Ry, int Nx, int Ny, double dx, double dy, double dt, int geoflag)
{
	//clock_t startm, stopm;
	//double shift_time=0.0;
	//double gettin_values_time=0.0;
	//double kick_rotate_time=0.0;
	//CALCULATING NUMBER OF ELECTRON MACROPARTICLES
	int SIZE=ELECTRONS.size();
	double Dx2=pow(dx*(Nx-1)/2.0,2);
	double Dy2=pow(dy*(Ny-1)/2.0,2);
	for(int i=0;i<SIZE;++i)
		{
		//SHIFTING THE ELECTRON
		//shift_electron(ELECTRONS[i].x, ELECTRONS[i].y, ELECTRONS[i].Vx, ELECTRONS[i].Vy, dt);
		//startm=clock();
		ELECTRONS[i].Vx=0.0;
		ELECTRONS[i].y+=ELECTRONS[i].Vy*dt;
		//stopm=clock();
		//shift_time+=((double)stopm-startm)/CLOCKS_PER_SEC;
		if((geoflag==1 && (ELECTRONS[i].x*ELECTRONS[i].x/Rx/Rx+ELECTRONS[i].y*ELECTRONS[i].y/Ry/Ry)<1.0 ) ||
(geoflag==2 && pow(ELECTRONS[i].x,2)<Dx2 && pow(ELECTRONS[i].y,2)<Dy2 )
 
)			{
			//printf("inside\n");
			//startm=clock();
			double aydt=get_value(Ey2D, Ny, WEIGHTS[i].downleft, WEIGHTS[i].upleft, WEIGHTS[i].downright, WEIGHTS[i].upright,  WEIGHTS[i].low_index)/2.0;
			//stopm=clock();
			//gettin_values_time+=((double)stopm-startm)/CLOCKS_PER_SEC;
			//startm=clock();
			//ACCELERATING THE ELECTRON
			//accelerate_electron(axdt,aydt,ELECTRONS[i].Vx,ELECTRONS[i].Vy);
			ELECTRONS[i].Vy+=aydt;
			//ROTATE
			//rotate_velocity(ELECTRONS[i].Vx,ELECTRONS[i].Vz,cosi,sini);
			//ACCELERATING
			//accelerate_electron(axdt,aydt,ELECTRONS[i].Vx,ELECTRONS[i].Vy);
			ELECTRONS[i].Vy+=aydt;
			//stopm=clock();
			//kick_rotate_time+=((double)stopm-startm)/CLOCKS_PER_SEC;
			}
		}
		//printf("shift_time=%e\n",shift_time*1000);
		//printf("gettin_values_time=%e\n",gettin_values_time*1000);
		//printf("kick_rotate_time=%e\n\n",kick_rotate_time*1000);
	}

void check_x_zero(vector <electron> ELECTRONS)
	{
	int SIZE=ELECTRONS.size();
	for(int i=0;i<SIZE;i++)
		{
		if(ELECTRONS[i].x==0.0)
			{
			printf("Bitch please\n");
			break;
			}
		}
	};




