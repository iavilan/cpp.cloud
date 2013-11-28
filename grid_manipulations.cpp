#include <math.h>
#include <time.h>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include "electron_manipulations.h"
#include "grid_manipulations.h"

using namespace std;
#define EPS0 8.85e-12
#define PI 3.141592651

//THIS FILE CONTAINS DEFINITIONS OF OPERATIONS THAT CAN BE PERFORMED WITH GRID

//ADDITIONAL FUNCTION TO CALCULATE A DERIVATIVE HAVING 3 POINTS SEPARATED BY dx IN ONE OF THESE POINTS
//THIS FUNCTION IS USED FOR A RECTANGULAR PIPE SHAPE TO CALCULATE ELECTRIC FIELD FROM POTENTIAL ON THE GRID EDGES
double derivative(double V1, double V2, double V3, double dL, int num)
	{
	//num=-1 - left point; num=0 - middle point; num=1 - right point
	return ((V1-V2*2+V3)*num+(V3-V1)/2.0)/dL;
	}

//GET VALUE FROM POINT x,y
double get_value(double *rho, int Nx, int Ny, double dx, double dy, double x, double y)
	{
	if(fabs(x)<=dx*(Nx-1)/2 && fabs(y)<=dy*(Ny-1)/2)
		{
		//OBTAINING COORDINATES OF PARTICLE IN GRID SYSTEM
		double nx=(x+dx/2.0*(Nx-1))/dx;
		double ny=(y+dy/2.0*(Ny-1))/dy;
		//THE CLOSEST KNOTS WITH SMALLER COORDINATES
		int low_x=(int)(nx);
		int low_y=(int)(ny);
		//IF ON THE UPPER EDGE OF THE GRID (RECTANGULAR GEOMETRY)
		if(low_x==Nx)
			{
			low_x=low_x-1;
			}
		if(low_y==Ny)
			{
			low_y=low_y-1;
			}
		//THE CLOSEST KNOTS WITH HIGHER COORDINATES
		int high_x=low_x+1;
		int high_y=low_y+1;
		//DISTANCE TO THE LOWER KNOTS
		double dnx=nx-low_x;
		double dny=ny-low_y;
		//DISTANCE THE UPPER KNOTS
		double dnxh=1.0-dnx;
		double dnyh=1.0-dny;
		//RETURN WEIGHED AVERAGE OF THE NEARBY KNOTS
		return rho[low_y+Ny*low_x]*dnxh*dnyh+rho[low_y+Ny*high_x]*dnx*dnyh+rho[high_y+Ny*low_x]*dnxh*dny+rho[high_y+Ny*high_x]*dnx*dny;
		}
	else
		{
		//printf("Trying to get value outside the grid");
		return 0.0;
		}
	};

//GET VALUE FROM POINT x,y FOR GIVEN WEIGHT
double get_value(double *rho, int Nx, int Ny, double dx, double dy, double dl, double ul, double dr, double ur, int low_x, int low_y)
	{
	//if(fabs(x)<=dx*(Nx-1)/2 && fabs(y)<=dy*(Ny-1)/2)
	//	{
		//IF ON THE UPPER EDGE OF THE GRID (RECTANGULAR GEOMETRY)
		if(low_x==Nx)
			{
			low_x=low_x-1;
			}
		if(low_y==Ny)
			{
			low_y=low_y-1;
			}
		//THE CLOSEST KNOTS WITH HIGHER COORDINATES
		int high_x=low_x+1;
		int high_y=low_y+1;
		//RETURN WEIGHED AVERAGE OF THE NEARBY KNOTS
		return (rho[low_y+Ny*low_x]*dl+rho[low_y+Ny*high_x]*dr+rho[high_y+Ny*low_x]*ul+rho[high_y+Ny*high_x]*ur)*dx*dy;
	//	}
	//else
	//	{
	//	//printf("Trying to get value outside the grid");
	//	return 0.0;
	//	}
	};

//GET VALUE FROM POINT x,y FOR GIVEN WEIGHT WITH SIMPLIFIED INDEX
double get_value(double *rho, int Ny, double dl, double ul, double dr, double ur, int low_index)
	{
	//if(fabs(x)<=dx*(Nx-1)/2 && fabs(y)<=dy*(Ny-1)/2)
	//	{
		//IF ON THE UPPER EDGE OF THE GRID (RECTANGULAR GEOMETRY)
		//if(low_x==Nx)
		//	{
		//	low_x=low_x-1;
		//	}
		//if(low_y==Ny)
		//	{
		//	low_y=low_y-1;
		//	}
		//THE CLOSEST KNOTS WITH HIGHER COORDINATES
		//int high_x=low_x+1;
		//int high_y=low_y+1;
		//RETURN WEIGHED AVERAGE OF THE NEARBY KNOTS
		return (rho[low_index]*dl+rho[low_index+Ny]*dr+rho[low_index+1]*ul+rho[low_index+1+Ny]*ur);
	//	}
	//else
	//	{
	//	//printf("Trying to get value outside the grid");
	//	return 0.0;
	//	}
	};


//SET VALUE OF PARTICLE WITH x,y AND lambda CHARGE LINE DENSITY
void add_value(double *rho, int Nx, int Ny, double dx, double dy, double x, double y, double lambda)
	{

	//OBTAINING COORDINATES OF PARTICLE IN GRID SYSTEM
	double nx=(x+dx/2.0*(Nx-1))/dx;
	double ny=(y+dy/2.0*(Ny-1))/dy;
	if(fabs(x)<(Nx-1)*dx/2 && fabs(y)<(Ny-1)*dy/2)
		{
		//IDEALY SHOULD CHECK WETHER A PARTICLE IS INSIDE A GRID

		//THE CLOSEST KNOTS WITH SMALLER COORDINATES
		int low_x=(int)(nx);
		int low_y=(int)(ny);
		//THE CLOSEST KNOTS WITH HIGHER COORDINATES
		int high_x=low_x+1;
		int high_y=low_y+1;
		//DISTANCE TO THE LOWER KNOTS
		double dnx=nx-low_x;
		double dny=ny-low_y;
		//DISTANCE THE UPPER KNOTS
		double dnxh=1.0-dnx;
		double dnyh=1.0-dny;
		//PRODUCT OF dx AND dy
		double dxdy=dx*dy;
		rho[low_y+Ny*low_x]+=lambda*dnxh*dnyh/dxdy;
		rho[high_y+Ny*high_x]+=lambda*dnx*dny/dxdy;
		rho[high_y+Ny*low_x]+=lambda*dnxh*dny/dxdy;
		rho[low_y+Ny*high_x]+=lambda*dnx*dnyh/dxdy;
		}
	else
		{
		1+1;
		//printf("Trying to add value out of the grid %e %e\n",x,y);
		}
	};

//SET VALUE OF PARTICLE WITH x,y AND lambda CHARGE LINE DENSITY AND SAVING THE WEIGHT
void add_value(double *rho, int Nx, int Ny, double dx, double dy, double x, double y, double lambda, double &dl, double &ul, double &dr, double &ur, int &low_x, int &low_y)
	{

	//OBTAINING COORDINATES OF PARTICLE IN GRID SYSTEM
	double nx=(x+dx/2.0*(Nx-1))/dx;
	double ny=(y+dy/2.0*(Ny-1))/dy;
	if(fabs(x)<(Nx-1)*dx/2 && fabs(y)<(Ny-1)*dy/2)
		{
		//IDEALY SHOULD CHECK WETHER A PARTICLE IS INSIDE A GRID

		//THE CLOSEST KNOTS WITH SMALLER COORDINATES
		low_x=(int)(nx);
		low_y=(int)(ny);
		//THE CLOSEST KNOTS WITH HIGHER COORDINATES
		int high_x=low_x+1;
		int high_y=low_y+1;
		//DISTANCE TO THE LOWER KNOTS
		double dnx=nx-low_x;
		double dny=ny-low_y;
		//DISTANCE THE UPPER KNOTS
		double dnxh=1.0-dnx;
		double dnyh=1.0-dny;
		//PRODUCT OF dx AND dy
		double dxdy=dx*dy;
		dl=dnxh*dnyh/dxdy;
		ul=dnxh*dny/dxdy;
		dr=dnx*dnyh/dxdy;
		ur=dnx*dny/dxdy;
		rho[low_y+Ny*low_x]+=lambda*dl;
		rho[high_y+Ny*high_x]+=lambda*ur;
		rho[high_y+Ny*low_x]+=lambda*ul;
		rho[low_y+Ny*high_x]+=lambda*dr;
		}
	};

//SET VALUE OF PARTICLE WITH x,y AND lambda CHARGE LINE DENSITY AND SAVING THE WEIGHT SIMPLIFIED INDEXING
void add_value(double *rho, int Nx, int Ny, double dx, double dy, 
double x, double y, double lambda, 
double &dl, double &ul, double &dr, double &ur, int &low_index)
	{

	//OBTAINING COORDINATES OF PARTICLE IN GRID SYSTEM
	double nx=(x+dx/2.0*(Nx-1))/dx;
	double ny=(y+dy/2.0*(Ny-1))/dy;
	if(nx>=0 && ny>=0 && nx<Nx && ny<Ny)
		{
		//THE CLOSEST KNOTS WITH SMALLER COORDINATES
		int low_x=(int)(nx);
		int low_y=(int)(ny);
		low_index=low_x*Ny+low_y;
		//THE CLOSEST KNOTS WITH HIGHER COORDINATES
		//int high_x=low_x+1;
		//int high_y=low_y+1;
		//DISTANCE TO THE LOWER KNOTS
		double dnx=nx-low_x;
		double dny=ny-low_y;
		//DISTANCE THE UPPER KNOTS
		double dnxh=1.0-dnx;
		double dnyh=1.0-dny;
		//PRODUCT OF dx AND dy
		dl=dnxh*dnyh;
		ul=dnxh*dny;
		dr=dnx*dnyh;
		ur=dnx*dny;
		rho[low_index]+=lambda*dl;
		rho[low_index+Ny+1]+=lambda*ur;
		rho[low_index+1]+=lambda*ul;
		rho[low_index+Ny]+=lambda*dr;
		}
	};

//SUM OF THE GRIDS
void sum_grids(double *rho1,double *rho2, int NxNy)
	{
	for(int i=0; i<NxNy; ++i)
		{
		//printf("Going well?\n");
		rho1[i]=rho1[i]+rho2[i];//rho2[i];
		}
	}

//MULTIPLY THE GRID BY NUMBER
double *grid_by_number(double *rho, double number, int NxNy)
	{
	double *temp=(double*) malloc(NxNy*sizeof(double));
	for(int i=0; i<NxNy; ++i)
		{
		temp[i]=rho[i]*number;
		}
	return temp;
	}

//SET ELLIPSE
void set_ellipse(double *rho,int Nx, int Ny, double dx, double dy, double Ax, double Ay, double lambda0)
	{
	double rho0=lambda0/Ax/Ay/PI;
	double ix=Ax/dx;
	double iy=Ay/dy;
	for(int i=0;i<Nx;++i)
		{
		for(int j=0;j<Ny;++j)
			{
			double xx=((double)i-(double)Nx/2.0)/ix;
			double yy=((double)j-(double)Ny/2.0)/iy;;
			double r=xx*xx+yy*yy;
			if(r<1.0)
				{
				rho[j*Nx+i]=rho0;
				}
			else
				{
				rho[j*Nx+i]=0.0;
				}
			}
		}	
	};

//SETTING GRID VALUES TO ZERO
void set_zero(double *rho, int NxNy)
	{
	for(int i=0;i<NxNy;++i)
		{
		rho[i]=0.0;
		}
	}

//SETTING GRID VALUES TO ZERO
void set_zero(long double *rho, int NxNy)
	{
	for(int i=0;i<NxNy;++i)
		{
		rho[i]=0.0;
		}
	}


//SIMPLE LINEAR INTERPOLATION OF ELECTRON VECTORS ON THE GRID
void simple_interpolate(vector <double> xe, vector <double> ye, vector <double> lambda_mac, double *rho, int Nx, int Ny, double dx, double dy)
	{
	printf("Nummer 1\n");
	//SIZE OF ELECTRONS ARRAY
	int SIZE=xe.size();
	//SETTING GRID VALUES TO ZERO
	set_zero(rho,Nx*Ny);
	int residue=SIZE%4;
	SIZE=SIZE-residue;
	int size=SIZE/4;
	//INTERPOLATING
	for(int i=0;i<size;++i)
		{
		int BEGIN=4*i;
		int BEGIN1=4*i+1;
		int BEGIN2=4*i+2;
		int BEGIN3=4*i+3;
		add_value(rho, Nx, Ny, dx, dy, xe[BEGIN], ye[BEGIN], lambda_mac[BEGIN]);
		add_value(rho, Nx, Ny, dx, dy, xe[BEGIN1], ye[BEGIN1], lambda_mac[BEGIN1]);
		add_value(rho, Nx, Ny, dx, dy, xe[BEGIN2], ye[BEGIN2], lambda_mac[BEGIN2]);
		add_value(rho, Nx, Ny, dx, dy, xe[BEGIN3], ye[BEGIN3], lambda_mac[BEGIN3]);
		}
	for(int i=SIZE;i<SIZE+residue;i++)
		{
		add_value(rho, Nx, Ny, dx, dy, xe[i], ye[i], lambda_mac[i]);
		}
	}
//SIMPLE LINEAR INTERPOLATION OF ELECTRON VECTORS ON THE GRID USING ELECTRONS STRUCTURE
void simple_interpolate(vector <electron> ELECTRONS, double *rho, int Nx, int Ny, double dx, double dy)
	{
	printf("Nummer 2\n");
	//SIZE OF ELECTRONS ARRAY
	int SIZE=ELECTRONS.size();
	//SETTING GRID VALUES TO ZERO
	set_zero(rho,Nx*Ny);
	int residue=SIZE%4;
	SIZE=SIZE-residue;
	int size=SIZE/4;
	//INTERPOLATING
	for(int i=0;i<size;++i)
		{
		int BEGIN=4*i;
		int BEGIN1=4*i+1;
		int BEGIN2=4*i+2;
		int BEGIN3=4*i+3;
		add_value(rho, Nx, Ny, dx, dy, ELECTRONS[BEGIN].x, ELECTRONS[BEGIN].y, ELECTRONS[BEGIN].lambda_mac);
		add_value(rho, Nx, Ny, dx, dy, ELECTRONS[BEGIN1].x, ELECTRONS[BEGIN1].y, ELECTRONS[BEGIN1].lambda_mac);
		add_value(rho, Nx, Ny, dx, dy, ELECTRONS[BEGIN2].x, ELECTRONS[BEGIN2].y, ELECTRONS[BEGIN2].lambda_mac);
		add_value(rho, Nx, Ny, dx, dy, ELECTRONS[BEGIN3].x, ELECTRONS[BEGIN3].y, ELECTRONS[BEGIN3].lambda_mac);
		}
	for(int i=SIZE;i<SIZE+residue;i++)
		{
		add_value(rho, Nx, Ny, dx, dy, ELECTRONS[i].x, ELECTRONS[i].y, ELECTRONS[i].lambda_mac);
		}
	}

//SIMPLE LINEAR INTERPOLATION OF ELECTRON VECTORS ON THE GRID STORING DATA FOR THE KICK
void simple_interpolate(vector <double> xe, vector <double> ye, vector <double> lambda_mac, double *rho, int Nx, int Ny, double dx, double dy,vector <double> &downleft,vector <double> &upleft,vector <double> &downright,vector <double> &upright,vector <int> &low_x,vector <int> &low_y)
	{
	printf("Nummer 3\n");
	//SIZE OF ELECTRONS ARRAY
	int SIZE=xe.size();
	downleft.assign(SIZE,0);
	upleft.assign(SIZE,0);
	downright.assign(SIZE,0);
	upright.assign(SIZE,0);
	low_x.assign(SIZE,0);
	low_y.assign(SIZE,0);
	//SETTING GRID VALUES TO ZERO
	set_zero(rho,Nx*Ny);
	int residue=SIZE%4;
	SIZE=SIZE-residue;
	int size=SIZE/4;
	//INTERPOLATING
	for(int i=0;i<size;++i)
		{
		int BEGIN=4*i;
		int BEGIN1=4*i+1;
		int BEGIN2=4*i+2;
		int BEGIN3=4*i+3;
		add_value(rho, Nx, Ny, dx, dy, xe[BEGIN], ye[BEGIN], lambda_mac[BEGIN],downleft[BEGIN],upleft[BEGIN],downright[BEGIN],upright[BEGIN],low_x[BEGIN],low_y[BEGIN]);
		add_value(rho, Nx, Ny, dx, dy, xe[BEGIN1], ye[BEGIN1], lambda_mac[BEGIN1],downleft[BEGIN1],upleft[BEGIN1],downright[BEGIN1],upright[BEGIN1],low_x[BEGIN1],low_y[BEGIN1]);
		add_value(rho, Nx, Ny, dx, dy, xe[BEGIN2], ye[BEGIN2], lambda_mac[BEGIN2],downleft[BEGIN2],upleft[BEGIN2],downright[BEGIN2],upright[BEGIN2],low_x[BEGIN2],low_y[BEGIN2]);
		add_value(rho, Nx, Ny, dx, dy, xe[BEGIN3], ye[BEGIN3], lambda_mac[BEGIN3],downleft[BEGIN3],upleft[BEGIN3],downright[BEGIN3],upright[BEGIN3],low_x[BEGIN3],low_y[BEGIN3]);
		}
	for(int i=SIZE;i<SIZE+residue;++i)
		{
		add_value(rho, Nx, Ny, dx, dy, xe[i], ye[i], lambda_mac[i],downleft[i],upleft[i],downright[i],upright[i],low_x[i],low_y[i]);
		}
	}
//SIMPLE LINEAR INTERPOLATION OF ELECTRON VECTORS ON THE GRID STORING DATA FOR THE KICK USING ELECTRON STRUCTURE
void simple_interpolate(vector <electron> &ELECTRONS, double *rho, int Nx, int Ny, double dx, double dy,
//vector <double> &downleft,vector <double> &upleft,vector <double> &downright,vector <double> &upright,vector <int> &low_x,vector <int> &low_y
vector <weight> &WEIGHTS
)
	{
	//printf("Simply interpolating\n");
	//SIZE OF ELECTRONS ARRAY
	int SIZE=ELECTRONS.size();
	weight temp;
	WEIGHTS.assign(SIZE,temp);
	//upleft.assign(SIZE,0);
	//downright.assign(SIZE,0);
	//upright.assign(SIZE,0);
	//low_x.assign(SIZE,0);
	//low_y.assign(SIZE,0);
	//SETTING GRID VALUES TO ZERO
	set_zero(rho,Nx*Ny);
	int residue=SIZE%4;
	SIZE=SIZE-residue;
	int size=SIZE/4;
	//INTERPOLATING
	for(int i=0;i<size;++i)
		{
		int BEGIN=4*i;
		int BEGIN1=4*i+1;
		int BEGIN2=4*i+2;
		int BEGIN3=4*i+3;
		add_value(rho, Nx, Ny, dx, dy, ELECTRONS[BEGIN].x, ELECTRONS[BEGIN].y, ELECTRONS[BEGIN].lambda_mac, WEIGHTS[BEGIN].downleft, WEIGHTS[BEGIN].upleft, WEIGHTS[BEGIN].downright, WEIGHTS[BEGIN].upright, WEIGHTS[BEGIN].low_index);
		add_value(rho, Nx, Ny, dx, dy, ELECTRONS[BEGIN1].x, ELECTRONS[BEGIN1].y, ELECTRONS[BEGIN1].lambda_mac, WEIGHTS[BEGIN1].downleft, WEIGHTS[BEGIN1].upleft, WEIGHTS[BEGIN1].downright, WEIGHTS[BEGIN1].upright, WEIGHTS[BEGIN1].low_index);
		add_value(rho, Nx, Ny, dx, dy, ELECTRONS[BEGIN2].x, ELECTRONS[BEGIN2].y, ELECTRONS[BEGIN2].lambda_mac, WEIGHTS[BEGIN2].downleft, WEIGHTS[BEGIN2].upleft, WEIGHTS[BEGIN2].downright, WEIGHTS[BEGIN2].upright, WEIGHTS[BEGIN2].low_index);
		add_value(rho, Nx, Ny, dx, dy, ELECTRONS[BEGIN3].x, ELECTRONS[BEGIN3].y, ELECTRONS[BEGIN3].lambda_mac, WEIGHTS[BEGIN3].downleft, WEIGHTS[BEGIN3].upleft, WEIGHTS[BEGIN3].downright, WEIGHTS[BEGIN3].upright, WEIGHTS[BEGIN3].low_index);
		}
	for(int i=SIZE;i<SIZE+residue;++i)
		{
		add_value(rho, Nx, Ny, dx, dy, ELECTRONS[i].x, ELECTRONS[i].y, ELECTRONS[i].lambda_mac,WEIGHTS[i].downleft, WEIGHTS[i].upleft, WEIGHTS[i].downright, WEIGHTS[i].upright, WEIGHTS[i].low_index);
		}
	}


//INTERPOLATE ON A RADIALLY SYMMETRIC GRID
void radial_interpolate(vector <double> xe, vector <double> ye, vector <double> lambda_mac, double *rho_r, double Rr, int Nr, double dr)
	{
	//SIZE OF ELECTRONS ARRAY
	int SIZE=xe.size();
	//SETTING GRID VALUES TO ZERO
	set_zero(rho_r,Nr);
	//CONSTANTS FOR NORMALIZATION
	double consta=PI*dr*dr;
	//INTERPOLATING
	for(int i=0;i<SIZE;++i)
		{
		double rr=sqrt(xe[i]*xe[i]+ye[i]*ye[i]);
		int ind=(int)(rr/dr);
		if(ind>=Nr)
			{
			printf("Pepets!\n");
			}
		else
			{
			rho_r[ind]+=lambda_mac[i];
			}
		}
	for(int i=0;i<Nr;++i)
		{
		rho_r[i]=rho_r[i]/(consta*(1.0+2.0*i));
		}
	}
//INTERPOLATE ON A RADIALLY SYMMETRIC GRID USING ELECTRON STRUCTURE
void radial_interpolate(vector <electron> ELECTRONS, double *rho_r, double Rr, int Nr, double dr)
	{
	//SIZE OF ELECTRONS ARRAY
	int SIZE=ELECTRONS.size();
	//SETTING GRID VALUES TO ZERO
	set_zero(rho_r,Nr);
	//CONSTANTS FOR NORMALIZATION
	double consta=PI*dr*dr;
	//INTERPOLATING
	for(int i=0;i<SIZE;++i)
		{
		double rr=sqrt(ELECTRONS[i].x*ELECTRONS[i].x+ELECTRONS[i].y*ELECTRONS[i].y);
		int ind=(int)(rr/dr);
		if(ind>=Nr)
			{
			//printf("Pepets!\n");
			}
		else
			{
			rho_r[ind]+=ELECTRONS[i].lambda_mac;
			}
		}
	for(int i=0;i<Nr;++i)
		{
		rho_r[i]=rho_r[i]/(consta*(1.0+2.0*i));
		}
	}

//SAVE RADIAL DISTRIBUTION
void radial_save(FILE *fd, double *rho_r, double Rr, int Nr, double dr)
	{
	for(int i=0;i<Nr;++i)
		{
		fprintf(fd,"%e ",rho_r[i]);
		}
	fprintf(fd,"\n");
	}

//SAVE DISTRIBUTION IN AS A FUNCTION OF X IN j-th SLICE, j<Ny;
void x_save(FILE *fd, double *rho_e, int Nx, int Ny, int j)
	{
	for(int i=0;i<Nx;++i)
		{
		fprintf(fd,"%e ",rho_e[j+Ny*i]);
		}
	fprintf(fd,"\n");
	}
//SAVE DISTRIBUTION IN AS A FUNCTION OF Y IN j-th SLICE, j<Nx;
void y_save(FILE *fd, double *rho_e, int Nx, int Ny, int j)
	{
	for(int i=0;i<Ny;++i)
		{
		fprintf(fd,"%e ",rho_e[i+Ny*j]);
		}
	fprintf(fd,"\n");
	}




