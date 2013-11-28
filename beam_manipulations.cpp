#include <math.h>
#include <time.h>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include <omp.h>
#include "grid_manipulations.h"
#include "field_manipulations.h"
#include "electron_manipulations.h"
#include "beam_manipulations.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <stdio.h>
#include "/hera/bhs/fpetrov/POSTDOC_CODE/lapack/include/lapackpp/lapackd.h"

#define PI 3.1415926535897932384626433832795 //

using namespace std;

//THE FUNCTION USED FOR RIGID BEAM TO WEIGHT THE FIELD OF SLICES (LOCAL LINE DENSITY RELATIVE TO THE AVERAGE ONE)
vector <double> gaussian_weight(int Nslice,double dz, double sigmaz)
	{
	vector <double> temp_weight;
	double value=0.0;
	double summ=0.0;
	for(int i=0;i<Nslice;++i)
		{
		value=exp(-pow(dz*((double)i-(double)(Nslice)/2),2)/(2.0*pow(sigmaz,2)));
		summ+=value/Nslice;
		temp_weight.push_back(value);
		}
	for(int i=0;i<Nslice;++i)
		{
		temp_weight.at(i)/=summ;
		}
	return temp_weight;
	};
//THE FUNCTION USED FOR RIGID BEAM TO WEIGHT THE FIELD OF SLICES (LOCAL LINE DENSITY RELATIVE TO THE AVERAGE ONE)
vector <double> parabolic_weight(int Nslice,double dz, double totlength)
	{
	vector <double> temp_weight;
	double value=0.0;
	double summ=0.0;
	for(int i=0;i<Nslice;++i)
		{
		value=(dz*i-totlength/2-(double)Nslice/2)*(dz*i+totlength-(double)Nslice/2);
		if(value>0)
		{
		temp_weight.push_back(value);
		summ+=value/Nslice;
		}
		else
		{
		temp_weight.push_back(0.0);
		}
		}
	for(int i=0;i<Nslice;++i)
		{
		temp_weight.at(i)/=summ;
		}
	return temp_weight;
	};



//THIS FUNCTION CREATES A BEAM-PARTICLE DISTRIBUTION CORRESPONDING TO EACH SLICE OF FUNCTION gaussian_weight
//FORGET TO ADD TRANSVERSE MOMENTA
void populate(
//VECTORS USED DIRECTLY IN THE BEAM TRANSFER
vector < vector <double> > &xx, vector < vector <double> > &yy, vector < vector <double> > &Vxx, vector < vector <double> > &Vyy,
//FOCUSING PARAMETERS
double Kx, double Ky,
//OTHER PARAMETERS: TOTAL NUMBER OF SLICES, STEP IN LONGITUDINAL DIRECTION, RMS BUNCH LENGTH, BUNCH TRANSVERSE SIZE, RANDOM NUMBER GENERATOR
int Nslice, double dz, double sigmaz, double Ax, double Ay, gsl_rng *rr)
	{
	//NUMBER OF SLICES FOR WHICH THE BEAM PARTICLES ARE GENERATED //BASICALLY THE LENGTH IS 8 SIGMA
	int slices=4*(int)(2.0*sigmaz/dz);
	//CENTER
	double center=(double)Nslice*dz/2.0;
	//MAXIMUM DISTANCE FROM THE BUNCH MAXIMUM WHERE PARTICLE CAN STILL BE PRESENT
	double Dmax=dz*slices/2.0;
	//STARTING POINT
	double start_point=center-Dmax;
	double end_point=center+Dmax;
	//MINIMUM INDEX
	int start_index=(int)(floor(start_point/dz));
	int end_index=(int)(floor(end_point/dz));
	//FILLING THE XX AND YY VECTORS WITH EMPTY VALUES
	vector <double> temp;
	vector < vector <double> > tempx(slices,temp);
	vector < vector <double> > tempy(slices,temp);
	vector < vector <double> > tempVx(slices,temp);
	vector < vector <double> > tempVy(slices,temp);
	xx=tempx;
	yy=tempy;
	Vxx=tempVx;
	Vyy=tempVy;
	
	//NUMBER OF MACROPARTICLES USED FOR THE BEAM
	int NUM=20000;
	for(int j=0;j<NUM;++j)
		{
		//printf("%d\n",j);
		double coord=dz*Nslice;
		int index=slices;
		while(index>=slices)
		{
		
		while(coord<start_point || coord>end_point)
			{
			coord=gsl_ran_gaussian(rr,sigmaz)+center;
			}
		index=(int)(coord/dz)-start_index;
		}
		//printf("Coord %e %d %d\n",coord,index,tempx.size());
		double rad_distance=sqrt(gsl_rng_uniform(rr));
		double phase=8.0*PI*gsl_rng_uniform(rr);
		double phase2=8.0*PI*gsl_rng_uniform(rr);
		//printf("%d %d\n",index,slices);
		xx.at(index).push_back(Ax*sin(phase)*rad_distance);
		yy.at(index).push_back(Ay*cos(phase)*rad_distance);
		Vxx.at(index).push_back(Ax*sqrt(1.0-rad_distance*rad_distance)*Kx*sin(phase2));
		Vyy.at(index).push_back(Ay*sqrt(1.0-rad_distance*rad_distance)*Ky*cos(phase2));
		
		}
	}


void kick_particle(double dxs, double dys, double &Vx, double &Vy)
	{
	Vx+=dxs;
	Vy+=dys;
	}
//FOR GIVEN 2D TRANSVERSE FIELD KICK THE PARTICLES
void kick_slice(vector <double> &x, vector <double> &y, vector <double> &Vx, vector <double> &Vy, double *Ex2D, double *Ey2D, int Nx, int Ny, double dx, double dy)
	{
	int SIZE=x.size();
	for(int i=0;i<SIZE;i++)
		{
		double dxs=get_value(Ex2D,Nx,Ny,dx,dy,x.at(i),y.at(i));
		double dys=get_value(Ey2D,Nx,Ny,dx,dy,x.at(i),y.at(i));
		//THIS FUNCTION IS TO IMPLEMENT
		kick_particle(dxs,dys,Vx.at(i),Vy.at(i));
		}
	}

//FOR GIVEN 3D TRANSVERSE FIELD KICK THE BEAM PARTICLES
void transverse_beam_kick_3d(
//PARTICLE DATA
vector < vector <double> > &xx, vector < vector <double> > &yy, 
vector < vector <double> > &Vxx, vector < vector <double> > &Vyy,
//FIELD DATA
vector < vector <double> > Ex, vector < vector <double> > Ey,
int Nx, int Ny, double dx, double dy
)
	{
	int SIZE=xx.size();
	for(int i=0;i<SIZE;i++)
		{
		kick_slice(xx.at(i),yy.at(i),Vxx.at(i),Vyy.at(i),&Ex.at(i)[0],&Ey.at(i)[0],Nx,Ny,dx,dy);
		}
	}


//CONSTANT FOCUSING TRANSFER
void cf_transfer(vector < vector <double> > &xx, vector < vector <double> > &yy, 
vector < vector <double> > &Vxx, vector < vector <double> > &Vyy,
double Kx, double Ky,
double cosx, double sinx,
double cosy, double siny,
vector <double> &x, vector <double> &y, vector <double> &Vx, vector <double> &Vy
)
	{
	int BIGSIZE=xx.size();
	if(x.size()==0)
		{
		for(int i=0;i<BIGSIZE;++i)
			{
			int SMALLSIZE=xx.at(i).size();
			for(int j=0;j<SMALLSIZE;++j)
				{
				//SAVING THE PARTICLE POSITION
				x.push_back(xx.at(i).at(j));
				y.push_back(yy.at(i).at(j));
				Vx.push_back(Vxx.at(i).at(j));
				Vy.push_back(Vyy.at(i).at(j));
				}
			}
		}
	for(int i=0;i<BIGSIZE;++i)
		{
		int SMALLSIZE=xx.at(i).size();
		for(int j=0;j<SMALLSIZE;++j)
			{
			double tempx=xx.at(i).at(j);
			double tempy=yy.at(i).at(j);
			double tempVx=Vxx.at(i).at(j);
			double tempVy=Vyy.at(i).at(j);
			//TRANSFERRING
			xx.at(i).at(j) = tempx*cosx+tempVx*sinx/Kx;
			Vxx.at(i).at(j) = -Kx*tempx*sinx+tempVx*cosx;
			yy.at(i).at(j) = tempy*cosy+tempVy*siny/Ky;
			Vyy.at(i).at(j) = -Ky*tempy*siny+tempVy*cosy;
			//SAVING THE PARTICLE POSITION
			x.push_back(xx.at(i).at(j));
			y.push_back(yy.at(i).at(j));
			Vx.push_back(Vxx.at(i).at(j));
			Vy.push_back(Vyy.at(i).at(j));
			}
		}
	}
//CONSTANT FOCUSING TRANSFER
void cf_transfer(
//BEAM PARTICLE PARAMETERS SAVED PER EACH SLICE
vector < vector <double> > &xx, vector < vector <double> > &yy, 
vector < vector <double> > &Vxx, vector < vector <double> > &Vyy,
double Kx, double Ky,
double cosx, double sinx,
double cosy, double siny,
//PARAMETERS OF A BEAM PARTICLES SAVED FOR TUNE SHIFT DIAGRAM CALCULATION
vector <double> &x, vector <double> &y, vector <double> &Vx, vector <double> &Vy,
//SLICE NUMBER
vector <int> &place,
//NUMBER OF KICKS PER ITERATION
int Nkicks,
//NUMBER OF THE KICK UNDER THE EXECUTION
int kick
)
	{
	int BIGSIZE=xx.size();
	if(x.size()==0)
		{
		for(int i=0;i<BIGSIZE;++i)
			{
			int SMALLSIZE=xx.at(i).size();
			for(int j=0;j<SMALLSIZE;++j)
				{
				//SAVING THE PARTICLE POSITION
				x.push_back(xx.at(i).at(j));
				y.push_back(yy.at(i).at(j));
				Vx.push_back(Vxx.at(i).at(j));
				Vy.push_back(Vyy.at(i).at(j));
				place.push_back(i);
				}
			}
		}
	for(int i=0;i<BIGSIZE;++i)
		{
		int SMALLSIZE=xx.at(i).size();
		for(int j=0;j<SMALLSIZE;++j)
			{
			double tempx=xx.at(i).at(j);
			double tempy=yy.at(i).at(j);
			double tempVx=Vxx.at(i).at(j);
			double tempVy=Vyy.at(i).at(j);
			//TRANSFERRING
			xx.at(i).at(j) = tempx*cosx+tempVx*sinx/Kx;
			Vxx.at(i).at(j) = -Kx*tempx*sinx+tempVx*cosx;
			yy.at(i).at(j) = tempy*cosy+tempVy*siny/Ky;
			Vyy.at(i).at(j) = -Ky*tempy*siny+tempVy*cosy;
			//SAVING THE PARTICLE POSITION ONLY ONCE PER FULL REVOLUTION
			if(kick%Nkicks==Nkicks-1)
				{
				x.push_back(xx.at(i).at(j));
				y.push_back(yy.at(i).at(j));
				Vx.push_back(Vxx.at(i).at(j));
				Vy.push_back(Vyy.at(i).at(j));
				}
			}
		}
	}


//TUNESHIFT EVALUATION // iter - NUMBER OF ITERATIONS THAT THE BEAM WAS TRANSFERED THROUGH THE LATTICE // USES LAPACK
void cf_tune_shift_evaluate(
//INPUT PARAMETERS
vector <double> x, vector <double> y, vector <double> Vx, vector <double> Vy, 
double fQx, double fQy,
int iter,
//OUTPUT PARAMETERS
vector <double> &dQx, vector <double> &dQy
)
	{
	if(iter>5)
		{
		//printf("Evaluation of tune spread\n");
		//TOTAL LENGTH OF AN ARRAY
		int tot_length=x.size();
		//LENGTH OF PARTICLE ARRAY
		int one_turn_length=tot_length/(iter+1);
		//MATRICES PARAMETERS  A[m*n]x[n]=b[nrhs*m]
		//THE BIGGES MATRIX RANK
		integer m=iter;
		//NUMBER OF PARAMETERS TO FIND IN TRANSFER MATRIX
		integer n=2;
		//RESULT
		integer nrhs=1;
		//ALLOCATING MATRICES
		double *ax=(double*) malloc(m*n*sizeof(double));
		double *bx=(double*) malloc(m*nrhs*sizeof(double));
		double *ay=(double*) malloc(m*n*sizeof(double));
		double *by=(double*) malloc(m*nrhs*sizeof(double));
		double *aVx=(double*) malloc(m*n*sizeof(double));
		double *bVx=(double*) malloc(m*nrhs*sizeof(double));
		double *aVy=(double*) malloc(m*n*sizeof(double));
		double *bVy=(double*) malloc(m*nrhs*sizeof(double));
		//OTHER PARAMETERS NEEDED FOR LAPACK
		//LWORK -1 MEANS THAT THE FIRST FUNCTION CALL WILL EVALUATE THE NEEDED LWOR
		integer lwork=-1;
		integer lda=m;
		integer ldb=m;
		integer info;
		double wkopt;
		//THIS WILL BE THE ARRAY OF SIZE LWORK
		double *work;
		//WE GET OPTIMAL PARAMETERS FOR THE LEAST SQUARE METHOD
		char lapack_command[]="N";
		dgels_(lapack_command,&m,&n,&nrhs,ax,&lda,bx,&ldb,&wkopt,&lwork,&info);
		//WE CREATE VARIABLES BASED ON THE CALCULATED PARAMETERS
		lwork = (int)wkopt;
		work = (double*)malloc( lwork*sizeof(double) );
		//LEAST SQUARE METHOD FOR EACH PARTICLE
		for(int i=0;i<one_turn_length;i++)
			{
			//POPULATING MATRICES
			for(int j=0;j<m;j++)
				{
				//SOMETHING IS BROCKEN IN THIS PLACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				bx[j]=x.at((j+1)*one_turn_length+i);
				by[j]=y.at((j+1)*one_turn_length+i);
				bVx[j]=Vx.at((j+1)*one_turn_length+i);
				bVy[j]=Vy.at((j+1)*one_turn_length+i);
				//THE COORDINATES IN PHASE SPACE DEPEND ON THE COORDINATES IN THE PREVIOUS TURN. THAT IS WHY BEFOR WE HAD +1
				//COORDINATES OF THE PREVIOUS TURN
				/*ax[j]=x.at(j*one_turn_length+i);
				ax[j+m]=Vx.at(j*one_turn_length+i);
				ay[j]=y.at(j*one_turn_length+i);
				ay[j+m]=Vy.at(j*one_turn_length+i);
				aVx[j]=x.at(j*one_turn_length+i);
				aVx[j+m]=Vx.at(j*one_turn_length+i);
				aVy[j]=y.at(j*one_turn_length+i);
				aVy[j+m]=Vy.at(j*one_turn_length+i);*/
				//MATRIX ELEMENTS FOR X PRODUCTION
				ax[j]=x.at(j*one_turn_length+i);
				ax[j+m]=Vx.at(j*one_turn_length+i);
				//MATRIX ELEMENTS FOR Y PRODUCTION
				ay[j]=y.at(j*one_turn_length+i);
				ay[j+m]=Vy.at(j*one_turn_length+i);
				//MATRIX ELEMENTS FOR VX PRODUCTION
				aVx[j]=x.at(j*one_turn_length+i);
				aVx[j+m]=Vx.at(j*one_turn_length+i);
				//MATRIX ELEMENTS FOR VY PRODUCTION
				aVy[j]=y.at(j*one_turn_length+i);
				aVy[j+m]=Vy.at(j*one_turn_length+i);
				
				}
			dgels_(lapack_command,&m,&n,&nrhs,ax,&lda,bx,&ldb,work,&lwork,&info);
			dgels_(lapack_command,&m,&n,&nrhs,ay,&lda,by,&ldb,work,&lwork,&info);
			dgels_(lapack_command,&m,&n,&nrhs,aVx,&lda,bVx,&ldb,work,&lwork,&info);
			dgels_(lapack_command,&m,&n,&nrhs,aVy,&lda,bVy,&ldb,work,&lwork,&info);
			/*double tempx=acos((2.0*bx[0]*bVx[1]-1.0))/4.0/PI;
			double tempy=acos((2.0*by[0]*bVy[1]-1.0))/4.0/PI;
			double tempx1=acos((1.0+2.0*bx[1]*bVx[0]))/4.0/PI;
			double tempy1=acos((1.0+2.0*by[1]*bVy[0]))/4.0/PI;*/
			double tempx=acos((bx[0]+bVx[1])/2.0)/2.0/PI;
			double tempy=acos((by[0]+bVy[1])/2.0)/2.0/PI;
			/*if(bx[0]<0)
				{
				tempx=0.5-tempx;
				}
			if(by[0]<0)
				{
				tempy=0.5-tempy;
				}*/
			dQx.push_back(tempx-fQx);
			dQy.push_back(tempy-fQy);
			}
		free(ax);
		free(bx);
		free(ay);
		free(by);
		free(aVx);
		free(bVx);
		free(aVy);
		free(bVy);
		free(work);
		}
	else
		{
		printf("Number of turns/iterations is not enough <5");
		}
	};

