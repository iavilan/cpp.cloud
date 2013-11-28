#include "slu_ddefs.h"
#include <math.h>
#include <time.h>
#include <vector>
#include <cstring>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include "supermatrix.h"
#include "grid_manipulations.h"
#include "field_manipulations.h"
#include "/hera/bhs/fpetrov/POSTDOC_CODE/lapack/include/lapackpp/lapackd.h"

#define EPS0 8.85e-12
#define PI 3.1415926535897932384626433832795 //

using namespace std;
//THIS FILE CONTAINS DEFINITIONS OF OPERATIONS FOR THE FIELD EVALUATION

//THIS STRUCTURE IS PASSED TO THE SOLVER 





//GREEN FUNCTION FOR ROUND BOUNDARY
//x0,y0 - source point
//x,y - calculation point
//R - radius of the pipe
//GREEN FUNCTION FOR ROUND BOUNDARY IS A PROBLEM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FOR SOLVER
double green_round(double x0, double y0, double x, double y, double a0)
	{
	double r02=x0*x0+y0*y0;
	double r0=sqrt(r02);
	double a02=a0*a0;
	double R=sqrt((x0-x)*(x0-x)+(y0-y)*(y0-y));
	double Rbar = (r02*x-a02*x0)*(r02*x-a02*x0)+(r02*y-a02*y0)*(r02*y-a02*y0);
	if(Rbar==0.0 || R==0.0 || r0==0.0)
		{
		printf("Apocalypto Rbar=%e R=%e r0=%e x0=%e y0=%e x=%e y=%e\n",Rbar,R,r0,x0,y0,x,y);
		return 0.0;
		}
	//return (-log((R*R*a02*r02)/Rbar )) /(4.0*PI);
	//REPLACE FOR A WHILE WITH OPEN BOUNDARY
	return -log(R)/(2.0*PI);
	}
//THE SLOWEST SOLVER USING GREEN FUNCTION TO SOLVE THE POISSON EQUATION
void calc_potential_green(double *rhoe,int Nx,int Ny,double dx, double dy, double R)
	{
	int NxNy=Nx*Ny;
	double dxdy=dx*dy;
	double *rho_temp=(double*)malloc(NxNy*sizeof(double));
	set_zero(rho_temp,NxNy);
	double shiftx=dx*(Nx-1)/2.0;
	double shifty=dy*(Ny-1)/2.0;
	double Rcorner=sqrt(shiftx*shiftx+shifty*shifty);
	for(int i=0;i<Nx;++i)
		{
		for(int j=0;j<Ny;++j)
			{
			//SOURCE POINT
			double x0=dx*i-shiftx;
			double y0=dy*j-shifty;
			//double r0=x0
			if(rhoe[Nx*i+j]!=0.0)
				{
				for(int ii=0;ii<Nx;++ii)
					{
					for(int jj=0;jj<Ny;++jj)
						{
						//TARGET POINT
						double x=dx*ii-shiftx;
						double y=dy*jj-shifty;
						if((x!=x0 || y!=y0))
							{
							//POTENTIAL IN A TARGET POINT
							rho_temp[Nx*jj+ii]+=rhoe[Nx*i+j]*green_round(x0, y0, x, y, R)*dxdy/EPS0;
							}
						}
					}
				}
			}
		}
	memcpy(rhoe,rho_temp,NxNy*sizeof(double));
	free(rho_temp);
	};

//GENERATE 2D GREEN FUNCTION FOR THE GRID
void gen_green_round(double *greens, int *border_x, int *border_y ,double Radius, int Nx, int Ny, double dx, double dy)
	{
	//SOME TEMPORARY GRID PARAMETERS
	double dxdy=dx*dy;
	int NxNy=Nx*Ny;
	//COORDINATE SHIFTS
	double shiftx=dx*(double)(Nx-1)/2.0;
	double shifty=dy*(double)(Ny-1)/2.0;
	double Rcorner=sqrt(shiftx*shiftx+shifty*shifty);
	//ARRAYS OF BORDER POINTS
	int nx=Nx-1;
	int ny=Ny-1;
	int nxny=nx*ny;
	int border_length=2*Nx+2*Ny-4;
	set_zero(greens,border_length*NxNy);
	for(int i=0;i<Nx;++i)
		{
		border_x[i]=i;
		border_y[i]=0;
		border_x[i+Nx]=i;
		border_y[i+Nx]=Nx-1;
		//printf("%d ",(*border_x)[i]);
		}
	for(int i=1;i<Ny-1;++i)
		{
		border_y[2*Nx+i-1]=i;
		border_x[2*Nx+i-1]=0;
		border_y[2*Nx+Ny-3+i]=i;
		border_x[2*Nx+Ny-3+i]=Nx-1;
		}
	double Radius2=pow(Radius,2);
	for(int i=0;i<border_length;++i)
		{
		//printf("%d %d\n",border_x[i],border_y[i]);
		//CALCULATING SQUARED DISTANCE FROM CENTER TO THE BORDER POINT
		double xb=dx*border_x[i]-shiftx;
		//double xb2=pow(xb,2);
		double yb=dy*border_y[i]-shifty;
		//double yb2=pow(yb,2);
		//if(fabs(xb)!=shiftx-dx && fabs(yb)!=shifty-dy)
		//printf("%d ",i);
		double coef=0.0;
		if((border_x[i]==0 || border_x[i]==Nx-1) && (border_y[i]==0 || border_y[i]==Ny-1))
			{
			coef=4.0;
			}
		else if(border_x[i]==0 || border_x[i]==Nx-1)
			{
			coef=2.0*dy/dx;
			}
		else
			{
			coef=2.0*dx/dy;
			}
		for(int j=1;j<Nx-1;++j)
			{
			for(int k=1;k<Ny-1;++k)
				{
				double x0=dx*j-shiftx;
				double y0=dy*k-shifty;
				//WE OBTAIN THE DENSITY SO WE DONT NEED EPS0
				double tempo=green_round(x0, y0, xb, yb, Radius);

				greens[i*NxNy+k*Nx+j]=coef*tempo;
				//}
				if(greens[i*NxNy+k*Nx+j]!=greens[i*NxNy+k*Nx+j])
					{
					printf("Green function go crazy!\n");
					}
				}
			}
		
		}
	};

void set_boundary_condition(double *rho, double *greens, int *border_x, int *border_y, int Nx, int Ny, int geoflag)
	{
	if(geoflag==1)
		{
		int nx=Nx-1;
		int ny=Ny-1;
		int borderlength=2*Nx+2*Ny-4;
		int NxNy=Nx*Ny;
		int nxny=nx*ny;
		
		//WE NEED TO GENERATE POTENTIAL NOT USING THE BOUNDARY POINTS
					for(int i=0;i<borderlength;++i)
						{
					double tempo=0.0;
		for(int j=1;j<Nx-1;++j)
			{
			for(int k=1;k<Ny-1;++k)
				{
				if(rho[Nx*k+j]!=0.0)
					{
						double coef=0.0;
						//printf("%d",i);
						if (greens[NxNy*i+Nx*k+j]!=0.0)
							{
							tempo+=(greens[NxNy*i+Nx*k+j]*rho[Nx*k+j]);
							//rho[Nx*border_y[i]+border_x[i]]+=(rho[Nx*k+j]*greens[NxNy*i+Nx*k+j]);
							}

						}
					}
				}
			rho[Nx*border_y[i]+border_x[i]]=tempo;
			}
		}
	};


//GETS POTENTIAL FROM A GIVEN DENSITY, GRID DIMENSIONS, AND GRID SIZE
//COEFFICIENTS FOR POISSON SOLVER OF THE FORM: cosnx[j]=cos(PI*((double)j+1)/(Nx));
//FFTW PLANS FOR SINE TRANSFORM OF THE FORM
//	fplan=fftw_plan_r2r_2d(Ny,Nx,rho_e,inv_rho_e,FFTW_RODFT00,FFTW_RODFT00,FFTW_ESTIMATE);
//	bplan=fftw_plan_r2r_2d(Ny,Nx,inv_rho_e,rho_e,FFTW_RODFT00,FFTW_RODFT00,FFTW_ESTIMATE);
//NX-1, NY-1 NUMBER OF GRID CELLS
//dx=2*Rx/(Nx-1)
//*rho STORES DENSITY, *rho_inv IS USED TO STORE FOURIER TRANSFORM OF DENSITY
void calc_potential(double *rho, double *rho_inv, double *cosnx, double *cosny, int Nx, int Ny, double dx, double dy, fftw_plan forward, fftw_plan backward)
	{
	//EXECUTING FFTW, TRANSFORM IS WRITTEN TO rho_inv
	fftw_execute(forward);
	//CALCULATING SOME INTERMEDIATE VARIABLES
	double dxdx=dx*dx;
	double dydy=dy*dy;
	double Norma=1.0/(Nx+1)/(Ny+1)/4.0/EPS0;
	//GOING THROUGH THE TRANSFORM
	for(int i=0;i<Nx;++i)
		{
		for(int j=0;j<Ny;++j)
			{
			//GRID DIMENSIONS ARE NOT EQUAL IN X AND Y
			//if( (i!=0 && j!=0) && (i!=0 && j!=Ny-1) && (i!=Nx-1 && j!=Ny-1) && (i!=Nx-1 && j!=0))
				
				rho_inv[i+Nx*j]=-Norma*rho_inv[i+Nx*j]/((2.0-2.0*cosnx[i])/(dxdx)+(2.0-2.0*cosny[j])/(dydy));
				if(!(rho_inv[i+Nx*j]==0) && !(rho_inv[i+Nx*j]!=0))
					{
					printf("%e\n",rho_inv[i+Nx*j]);
					}
				
			}
		}
	//EXECUTING INVERSE SINE FFT TRANSFORM, WRITTEN TO rho
	fftw_execute(backward);
	};
//THIS FUNCTION USES LU DECOMPOSITION TO SOLVE THE ELECTRIC POTENTIAL IN RECTANGULAR GRID
void calc_potential_lu(double *rho, int Nx, int Ny, double dx, double dy)
	{
	//CREATE A MATRIX
	integer size=(Nx*Ny);
	doublereal *matr=(doublereal*)malloc(Nx*Nx*Ny*Ny*sizeof(doublereal));
	printf("Matrix generation\n");
	for(int j=0;j<size;++j)
		{
		printf("%d\n",j);
		for(int i=0;i<size;++i)
			{
			if(i==j)
				{
				matr[i*size+j]=4.0;
				}
			else if(i==j-1 || i==j+1 || i==j+Nx || i==j-Nx)
				{
				matr[i*size+j]=-1.0;
				}
			else
				{
				matr[i*size+j]=0.0;
				}
			}
		}
	printf("Matrix decomposition\n");
	integer m=size;
	integer n=size;
	integer lda=size;
	integer info;
	integer *ipiv=(integer*)malloc((Nx*Nx)*sizeof(integer));
	integer lwork=size*size;
	doublereal *work = (doublereal*)malloc(Nx*Nx*sizeof(doublereal));
	dgetrf_(&size, &size, matr, &lda, ipiv, &info);
	printf("Matrix solution\n");
	dgetri_(&size,rho,&size,ipiv,work,&lwork,&info);
	FILE *fd=fopen("poisson.dat","w");
	printf("Output to file\n");
	for(int i=0;i<size;++i)
		{
		printf("%d\n",i);
		for(int j=0;j<size;++j)
			{
			fprintf(fd,"%e ",(double)matr[i*size+j]);
			//if(matr[i*(Ny-2)+j]!=0.0)
			//	{
			//	printf("%e ",matr[i*(Ny-2)+j]);
			//	}
			}
		fprintf(fd,"\n");
		}
	fclose(fd);
	free(ipiv);
	free(matr);
	free(work);
	};

//PREPARATION OF LU DECOMPOSITION FOR A GIVEN 
void preparation_of_matrix(double Rx, double Ry,  long int Nx, long int Ny, double dx, double dy,
double *&A_lu, long int *&IPIV, 
//RELATIVE CUT CELL EDGE LENGTH
vector <double> &alphax, vector <double> &alphay,
//INDICES OF POINT INSIDE THE PHYSICAL DOMAIN
vector <int> &ix, vector <int> &iy,
//POINTS NEAR THE BORDER
vector <int> &bx, vector <int> &by,
//ANGLE BETWEEN THE X AND Y AXIS AND THE WALL NEAR THE CUT CELL
vector <double> &tanx, vector <double> &tany)
	{
	int ixmin=Nx*Ny;
	int iymin=Nx*Ny;
	int ixmax=0;
	int iymax=0;
	double Dx=dx*(Nx-1)/2.0;
	double Dy=dy*(Ny-1)/2.0;
	vector <int> incolumn;
	vector <double> alphatempx;
	vector <double> alphatempy;
	//SAVING THE INDICES OF OF THE POINTS INCIDE THE ELLIPSE
	//COLUMN
	for(int j=0;j<Nx;++j)
		{
		int num=0;
		//ROW
		for(int i=0;i<Ny;++i)
			{
			double xx=(dx*j-Dx);
			double yy=(dy*i-Dy);
			double x=fabs(xx);
			double y=fabs(yy);
			double r2=(x*x)/(Rx*Rx)+(y*y)/(Ry*Ry);
			if(r2<1.0)
				{
				ix.push_back(j);
				if(j<ixmin)
					{
					ixmin=j;
					}
				//printf("%d\n",j);
				iy.push_back(i);
				num++;
				//SEARCHING INSIDE-AND-NEAR-THE-BOUBNDARY INDICES
				if( (pow((x+dx)/Rx,2)+pow((y+dy)/Ry,2)>1.0) && (pow(x/Rx,2)+pow(y/Ry,2)<=1.0) )
					{
					//FINDING THE CUT CELL EDGE LENGTH
					if(pow((x+dx)/Rx,2)+pow((y)/Ry,2)>1.0)
						{
						alphax.push_back( (Rx*sqrt(1.0-y*y/Ry/Ry)-x)/dx);
						alphatempx.push_back(alphax.back());
						tanx.push_back(-Rx*y/Ry/Ry/sqrt(1.0-y*y/Ry/Ry));
						bx.push_back(j);
						by.push_back(i);
						//alphax.push_back(1.0);
						}
					else
						{
						alphax.push_back(1.0);
						alphatempx.push_back(alphax.back());
						bx.push_back(j);
						by.push_back(i);
						}
					if(pow((x)/Rx,2)+pow((y+dy)/Ry,2)>1.0)
						{
						alphay.push_back( (Ry*sqrt(1.0-x*x/Rx/Rx)-y)/dy);
						tanx.push_back(-Ry*x/Rx/Rx/sqrt(1.0-x*x/Rx/Rx));
						alphatempy.push_back(alphay.back());
						if(alphatempy.size()>bx.size())
						{
						bx.push_back(j);
						by.push_back(i);
						}
						//alphay.push_back(1.0);
						}
					else
						{
						alphay.push_back(1.0);
						alphatempy.push_back(alphay.back());
						if(alphatempy.size()>bx.size())
						{
						bx.push_back(j);
						by.push_back(i);
						}
						}

					//FIRST WE TRY LINEARIZED CUT CELL
					//WE START WITH POINTS INSIDE THE PIPE												
					}
				else
					{
					alphax.push_back(1.0);
					alphay.push_back(1.0);
					}
				}
			}
		if(num>0)
		{
		incolumn.push_back(num);
		}
		}
	//HETTING TOTAL SIZE OF THE POINTS
	long int totsize=ix.size();
	long int bandwidth=3*Ny+1;
	//A_lu=new double[totsize*bandwidth];//(double*)malloc(totsize*bandwidth*sizeof(double));
	A_lu=(double*) malloc (totsize*bandwidth * sizeof(double));
	//MAKING ZEROS EVERYWHERE
	set_zero(A_lu,totsize*bandwidth);
	//printf("%e\n",A_lu[totsize*bandwidth-1]);
	//ADDING DIAGONALS
	for(int i=0;i<totsize;++i)
		{
		//MAIN DIAGONAL
		A_lu[bandwidth-1-Ny+i*bandwidth]=EPS0*(2.0/dx/dx/alphax[i]+2.0/dy/dy/alphay[i]);//2.0;
		//LOWER SMALL -1 DIAGONAL
		if(i>0 && ix[i]==ix.at(i-1))
			{
			A_lu[bandwidth-2-Ny+i*bandwidth]=-EPS0*1.0*2.0/(1.0+alphay[i])/dy/dy;//2.0;
			}
		//UPPER SMALL -1 DIAGONAL
		if(i<totsize-1 && ix[i]==ix.at(i+1))
			{
			A_lu[bandwidth-Ny+i*bandwidth]=-EPS0*1.0*2.0/(1.0+alphay[i])/dy/dy;//2.0;
			}
		int index=ix[i]-ixmin;
		//printf("index %d ",index);
		//LOWER BIG DIAGONAL
		if(index>0)
			{
			int delta=(incolumn.at(index)+incolumn.at(index-1))/2;
			//printf("delta %d i-delta %d ",delta,i-delta);
			if(i-delta>=0 && iy[i]==iy.at(i-delta))
				{
				//printf("accepted %d %d %d",i-delta,iy[i],iy.at(i-delta));
				A_lu[bandwidth-1-Ny-delta+i*bandwidth]=-EPS0*1.0*2.0/(1.0+alphax[i])/dx/dx;//2.0;
				}
			}
		//printf("\n");
		//UPPER BIG DIAGONAL
		if(index<incolumn.size()-1)
			{
			int delta=(incolumn.at(index)+incolumn.at(index+1))/2;
			if((i+delta)<totsize && iy[i]==iy.at(i+delta))
				{
				A_lu[bandwidth-1-Ny+delta+i*bandwidth]=-EPS0*1.0*2.0/(1.0+alphax[i])/dx/dx;//2.0;
				}
			}
		}
	IPIV = new long int[totsize];//(long int*)malloc(totsize*sizeof(long int));
	long int INFO=-1;
	dgbtrf_(&totsize,&totsize,&Ny,&Ny,A_lu,&bandwidth,IPIV,&INFO);
	alphax=alphatempx;
	alphay=alphatempy;
	}

//PREPARATION OF MATRIX FOR POISSON SOLVER
/*void preparation_of_matrix_superlu(double Rx, double Ry,  long int Nx, long int Ny, double dx, double dy,
struct superlu_params &Par,
 vector <double> &alphax, vector <double> &alphay,
vector <int> &ix, vector <int> &iy,
vector <int> &bx, vector <int> &by,
//ANGLE BETWEEN THE X AND Y AXIS AND THE WALL NEAR THE CUT CELL
vector <double> &tanx, vector <double> &tany)
	{
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Enter main()");
#endif
     Defaults 
    Par.lwork = 0;
    Par.nrhs  = 1;
    Par.equil = YES;	
    Par.u     = 1.0;
    Par.trans = NOTRANS;

    set_default_options(&Par.options);

    Par.options.Equil = Par.equil;
    Par.options.DiagPivotThresh = Par.u;
    Par.options.Trans = Par.trans;


    if ( Par.lwork > 0 ) {
	Par.work = SUPERLU_MALLOC(Par.lwork);
	if ( !Par.work ) {
	    ABORT("DLINSOLX: cannot allocate work[]");
	}
    }
printf("Creating arrays\n");
vector <double> Aa;
vector <int> Arows;
vector <int> Acolumn;
vector <int> incolumn;
//TEMPORARY VARIABLE FOR THE CUT-CELL INFORMATION
vector <double> alphatempx;
vector <double> alphatempy;
double Dx=dx*(Nx-1)/2.0;
double Dy=dy*(Ny-1)/2.0;
//HORIZONTAL INDEXING
Acolumn.push_back(0);

for(int i=0;i<Nx;++i)
	{
	int num=0;
	for(int j=0;j<Ny;++j)
		{
		double xx=(dx*j-Dx);
		double yy=(dy*i-Dy);
		double x=fabs(xx);
		double y=fabs(yy);
		//TAKING ONLY THE POINT INSIDE THE PIPE
		if(pow(x/Rx,2)+pow(y/Ry,2)<1.0)
			{
			//ADDING INDICES TO TO VECTORS
			ix.push_back(i);
			iy.push_back(j);
			num++;
			//IF SHIFT ALONG X OR ALONG Y PER ONE CELL GETS US OUTSIDE THE PIPE
			if( (pow((x+dx)/Rx,2)+pow((y)/Ry,2)>=1.0) || (pow((x)/Rx,2)+pow((y+dy)/Ry,2)>=1.0) )
					{
					//FINDING THE HORIZONTAL CUT CELL EDGE LENGTH
					if(pow((x+dx)/Rx,2)+pow((y)/Ry,2)>=1.0)
						{
						alphax.push_back( (Rx*sqrt(1.0-y*y/Ry/Ry)-x)/dx);
						//HERE THE  yy IS USED BECAUSE THE SIGN MATTERS
						tanx.push_back(-Rx*yy/Ry/Ry/sqrt(1.0-y*y/Ry/Ry));
						}
					else
						{
						alphax.push_back(1.0);
						}
					//FINDING THE VERTTICAL CUT CELL EDGE LENGTH
					if(pow((x)/Rx,2)+pow((y+dy)/Ry,2)>=1.0)
						{
						alphay.push_back( (Ry*sqrt(1.0-x*x/Rx/Rx)-y)/dy);
						//HERE THE  xx IS USED BECAUSE THE SIGN MATTERS
						tany.push_back(-Ry*xx/Ry/Ry/sqrt(1.0-x*x/Rx/Rx));
						}
					else
						{
						alphay.push_back(1.0);
						}

					if(tanx.size()>tany.size())
						{
						//JUST SOME VALUE THAT WILL NOT BE USED
						tany.push_back(1.0);
						}
					else if(tany.size()>tanx.size())
						{
						//JUST SOME VALUE THAT WILL NOT BE USED
						tanx.push_back(1.0);
						}
					
					alphatempx.push_back(alphax.back());
					alphatempy.push_back(alphay.back());
					bx.push_back(i);
					by.push_back(j);							
					}
				else
					{
					alphax.push_back(1.0);
					alphay.push_back(1.0);
					}
			}
		}
	if(num>0)
		{
		incolumn.push_back(num);
		}
	}

int totsize=ix.size();
int ixmin=ix.at(0);
double dxdx=dx*dx;
double dydy=dy*dy;
//CONSTRUCTING MATRIX X
for(int i=0;i<totsize;++i)
	{
	int index=ix[i]-ixmin;
	//UPPER DIAGONAL
	if(index>0)
		{
		int deltam=(incolumn.at(index)+incolumn.at(index-1))/2;
		if(i-deltam>=0 && iy[i]==iy.at(i-deltam))
			{
			Aa.push_back(-EPS0*2.0/dxdx/(1.0+alphax[i]));
			Arows.push_back(i-deltam);
			}
		}
	//UP DIAGONAL
	if(i>0 && ix[i]==ix.at(i-1))
		{
		//printf("Before up\n");
		Aa.push_back(-EPS0*2.0/dydy/(1.0+alphay[i]));
		Arows.push_back(i-1);
		}
	//MAIN DIAGONAL
		Aa.push_back(2.0*EPS0*(1.0/dxdx/alphax[i]+1.0/dydy/alphay[i]));
		Arows.push_back(i);
	//LOW DIAGONAL
	if(i<totsize-1 && ix[i]==ix.at(i+1))
		{
		Aa.push_back(-EPS0*2.0/dydy/(1.0+alphay[i]));
		Arows.push_back(i+1);
		}
	//LOWER DIAGONAL
	if(index<incolumn.size()-1)
		{
		int deltap=(incolumn.at(index)+incolumn.at(index+1))/2;
		if(i+deltap<totsize && iy[i]==iy.at(i+deltap))
			{
			Aa.push_back(-EPS0*2.0/dxdx/(1.0+alphax[i]));
			Arows.push_back(i+deltap);
			}
		}
	Acolumn.push_back(Aa.size());
	}
	printf("Rewriting vectors to dynamic arrays\n");
	double *a=new double[Aa.size()];
	memcpy(a,&Aa[0],Aa.size()*sizeof(double));
	int *asub=new int[Arows.size()];
	memcpy(asub,&Arows[0],Arows.size()*sizeof(int));
	int *xa=new int[Acolumn.size()];
	memcpy(xa,&Acolumn[0],Acolumn.size()*sizeof(int));

Par.nnz=Aa.size();
Par.rhsb=new double[totsize];
Par.rhsx=new double[totsize];
for(int j=0;j<totsize;++j)
	{
	Par.rhsb[j]=1.0;
	Par.rhsx[j]=1.0;
	}
    Par.m=totsize;
    Par.n=totsize;
	printf("Creating matrices in SuperLU format\n");
    dCreate_CompCol_Matrix(&(Par.A), Par.m, Par.n, Par.nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&(Par.B), Par.m, Par.nrhs, Par.rhsb, Par.m, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&(Par.X), Par.m, Par.nrhs, Par.rhsx, Par.m, SLU_DN, SLU_D, SLU_GE);
    Par.xact = doubleMalloc(Par.n * Par.nrhs);
    Par.ldx = Par.n;
    //dGenXtrue(n, nrhs, xact, ldx);
    //dFillRHS(trans, nrhs, xact, ldx, &A, &B);

	printf("Other SuperLU variables\n");
    if ( !(Par.etree = intMalloc(Par.n)) ) ABORT("Malloc fails for etree[].");
    if ( !(Par.perm_r = intMalloc(Par.m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(Par.perm_c = intMalloc(Par.n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(Par.R = (double *) SUPERLU_MALLOC(Par.A.nrow * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(Par.C = (double *) SUPERLU_MALLOC(Par.A.ncol * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for C[].");
    if ( !(Par.ferr = (double *) SUPERLU_MALLOC(Par.nrhs * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(Par.berr = (double *) SUPERLU_MALLOC(Par.nrhs * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for berr[].");
   //ONLY PERFORM THE LU DECOMPOSITION 
    Par.B.ncol = Par.nrhs;  //Indicate not to solve the system 
	printf("Finally LU decomposition\n");
    StatInit(&Par.stat);
    dgssvx(&(Par.options), &(Par.A), Par.perm_c, Par.perm_r, Par.etree, Par.equed, Par.R, Par.C,
           &(Par.L), &(Par.U), Par.work, Par.lwork, &(Par.B), &(Par.X), &(Par.rpg), &(Par.rcond), Par.ferr, Par.berr,
           &(Par.mem_usage), &(Par.stat), &(Par.info));

    StatFree(&(Par.stat));
    Par.options.Fact = FACTORED; // Indicate the factored form of A is supplied.
    Par.B.ncol = Par.nrhs;  // Set the number of right-hand side
	printf("Freeing stuff\n");
	delete [] xa;
	delete [] asub;
	delete [] a;
	printf("ix_size %d alpha_size %d\n",ix.size(),alphax.size());
	alphax=alphatempx;
	alphay=alphatempy;
	printf("New alphax_size %d\n",alphax.size());
	printf("Size of tangent vectors %d %d\n",tanx.size(),tany.size());
	}*/

void preparation_of_matrix_superlu(double Rx, double Ry,  long int Nx, long int Ny, double dx, double dy,
struct superlu_params &Par,
 vector <double> &alphax, vector <double> &alphay,
vector <int> &iindex,
vector <int> &bindex,
//ANGLE BETWEEN THE X AND Y AXIS AND THE WALL NEAR THE CUT CELL
vector <double> &tanx, vector <double> &tany)
	{
double dxdy=dx*dy;
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Enter main()");
#endif
//     Defaults 
    Par.lwork = 0;
    Par.nrhs  = 1;
    Par.equil = YES;
    //Par.IterRefine = ::DOUBLE;	
    Par.u     = 0.1;
    Par.trans = NOTRANS;
    set_default_options(&Par.options);

    Par.options.Equil = Par.equil;
    Par.options.DiagPivotThresh = Par.u;
   // Par.options.IterRefine = SLU_DOUBLE;
    Par.options.Trans = Par.trans;
    Par.options.ColPerm=MMD_AT_PLUS_A;
    Par.options.SymmetricMode = YES;


    if ( Par.lwork > 0 ) {
	Par.work = SUPERLU_MALLOC(Par.lwork);
	if ( !Par.work ) {
	    ABORT("DLINSOLX: cannot allocate work[]");
	}
    }
printf("Creating arrays\n");
vector <double> Aa;
vector <int> Arows;
vector <int> Acolumn;
vector <int> incolumn;
vector <int> ix;
vector <int> iy;
//TEMPORARY VARIABLE FOR THE CUT-CELL INFORMATION
vector <double> alphatempx;
vector <double> alphatempy;
double Dx=dx*(Nx-1)/2.0;
double Dy=dy*(Ny-1)/2.0;
//HORIZONTAL INDEXING
Acolumn.push_back(0);

for(int i=0;i<Nx;++i)
	{
	int num=0;
	for(int j=0;j<Ny;++j)
		{
		double xx=(dx*i-Dx);
		double yy=(dy*j-Dy);
		double x=fabs(xx);
		double y=fabs(yy);
		//TAKING ONLY THE POINT INSIDE THE PIPE
		if(pow(x/Rx,2)+pow(y/Ry,2)<=1.0)
			{
			//ADDING INDICES TO TO VECTORS
			ix.push_back(i);
			iy.push_back(j);
			iindex.push_back(i*Ny+j);
			num++;
			//IF SHIFT ALONG X OR ALONG Y PER ONE CELL GETS US OUTSIDE THE PIPE
			if( (pow((x+dx)/Rx,2)+pow((y)/Ry,2)>1.0) || (pow((x)/Rx,2)+pow((y+dy)/Ry,2)>1.0) )
					{
					//FINDING THE HORIZONTAL CUT CELL EDGE LENGTH
					if(pow((x+dx)/Rx,2)+pow((y)/Ry,2)>=1.0)
						{
						double alphatemp=(Rx*sqrt(1.0-y*y/Ry/Ry)-x)/dx;
						/*if(alphatemp==1.0)
							{
							alphatemp=0.0;
							}*/
						alphax.push_back( alphatemp);
						//HERE THE  yy IS USED BECAUSE THE SIGN MATTERS
						tanx.push_back(-Rx*y/Ry/Ry/sqrt(1.0-y*y/Ry/Ry));
						//alphax.push_back(1.0);
						}
					else
						{
						alphax.push_back(1.0);
						}
					//FINDING THE VERTTICAL CUT CELL EDGE LENGTH
					if(pow((x)/Rx,2)+pow((y+dy)/Ry,2)>=1.0)
						{
						double alphatemp=(Ry*sqrt(1.0-x*x/Rx/Rx)-y)/dy;
						/*if(alphatemp==1.0)
							{
							alphatemp=0.0;
							}*/
						alphay.push_back( alphatemp);
						//alphay.push_back(1.0);
						//HERE THE  xx IS USED BECAUSE THE SIGN MATTERS
						tany.push_back(-Ry*x/Rx/Rx/sqrt(1.0-x*x/Rx/Rx));
						}
					else
						{
						alphay.push_back(1.0);
						}

					if(tanx.size()>tany.size())
						{
						//JUST SOME VALUE THAT WILL NOT BE USED
						tany.push_back(1.0);
						}
					else if(tany.size()>tanx.size())
						{
						//JUST SOME VALUE THAT WILL NOT BE USED
						tanx.push_back(1.0);
						}
					
					alphatempx.push_back(alphax.back());
					alphatempy.push_back(alphay.back());
					//bx.push_back(i);
					//by.push_back(j);	
					bindex.push_back(i*Ny+j);						
					}
				else
					{
					alphax.push_back(1.0);
					alphay.push_back(1.0);
					}
			}
		}
	if(num>0)
		{
		incolumn.push_back(num);
		}
	}

int totsize=ix.size();
int ixmin=ix.at(0);
double dxdx=dx*dx;
double dydy=dy*dy;
//CONSTRUCTING MATRIX X
//vector <double> matrix(totsize*totsize,0.0);
for(int i=0;i<totsize;++i)
	{
	int index=ix[i]-ixmin;
	//UPPER DIAGONAL
	if(index>0)
		{
		int deltam=(incolumn.at(index)+incolumn.at(index-1))/2;
		//IF SHIFT BY deltam PRESERVES THE ROW NUMBER
		if(i-deltam>=0 && iy[i]==iy.at(i-deltam))
			{
			Aa.push_back(-EPS0*dxdy/dxdx/*/(1.0+alphax[i])*/);
			Arows.push_back(i-deltam);
			//matrix.at(i*totsize+i-deltam)=Aa.back()/EPS0;
			}
		}
	//UP DIAGONAL//IF SHIFT BY ONE ELEMENT PRESERVES THE COLUMN
	if(i>0 && ix[i]==ix.at(i-1))
		{
		//printf("Before up\n");
		Aa.push_back(-EPS0*dxdy/dydy/*/(1.0+alphay[i])*/);
		Arows.push_back(i-1);
		//matrix.at(i*totsize+i-1)=Aa.back()/EPS0/dydy;
		}
	//MAIN DIAGONAL
		//Aa.push_back(2.0*EPS0*(1.0/dxdx/alphax[i]+1.0/dydy/alphay[i]));
		//Aa.push_back(EPS0*(2.0/dxdx/alphax[i]+2.0/dydy/alphay[i]) );
		Aa.push_back(EPS0*dxdy*((1.0+alphax[i])/alphax[i]/dxdx+(1.0+alphay[i])/alphay[i]/dydy) );
		Arows.push_back(i);
		//matrix.at(i*totsize+i)=Aa.back()/EPS0;
	//LOW DIAGONAL//IF SHIFT BY ONE ELEMENT PRESERVES THE COLUMN
	if(i+1<totsize && ix[i]==ix.at(i+1))
		{
		Aa.push_back(-EPS0*dxdy/dydy/*/(1.0+alphay[i])*/);
		Arows.push_back(i+1);
		//matrix.at(i*totsize+i+1)=Aa.back()/EPS0;
		}
	//LOWER DIAGONAL
	if(index+1<incolumn.size())
		{
		int deltap=(incolumn.at(index)+incolumn.at(index+1))/2;
		//IF SHIFT BY deltap PRESERVES THE ROW NUMBER
		if(i+deltap<totsize && iy[i]==iy.at(i+deltap))
			{
			Aa.push_back(-EPS0*dxdy/dxdx/*/(1.0+alphax[i])*/);
			Arows.push_back(i+deltap);
			//matrix.at(i*totsize+i+deltap)=Aa.back()/EPS0;
			}
		}
	Acolumn.push_back(Aa.size());
	}
	/*if(totsize<=10)
	{
	FILE *fd=fopen("matrice.dat","w");
	for(int i=0;i<totsize;i++)
		{
		for(int j=0;j<totsize;j++)
			{
			fprintf(fd,"%e ",matrix[i*totsize+j]);
			}
		fprintf(fd, "\n");
		}
	fclose(fd);
	}*/
	printf("Rewriting vectors to dynamic arrays\n");
	double *a=new double[Aa.size()];
	memcpy(a,&Aa[0],Aa.size()*sizeof(double));
	int *asub=new int[Arows.size()];
	memcpy(asub,&Arows[0],Arows.size()*sizeof(int));
	int *xa=new int[Acolumn.size()];
	memcpy(xa,&Acolumn[0],Acolumn.size()*sizeof(int));

Par.nnz=Aa.size();
Par.rhsb=new double[totsize];
Par.rhsx=new double[totsize];
for(int j=0;j<totsize;++j)
	{
	Par.rhsb[j]=1.0;
	Par.rhsx[j]=1.0;
	}
    Par.m=totsize;
    Par.n=totsize;
	printf("Creating matrices in SuperLU format\n");
    dCreate_CompCol_Matrix(&(Par.A), Par.m, Par.n, Par.nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&(Par.B), Par.m, Par.nrhs, Par.rhsb, Par.m, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&(Par.X), Par.m, Par.nrhs, Par.rhsx, Par.m, SLU_DN, SLU_D, SLU_GE);
    Par.xact = doubleMalloc(Par.n * Par.nrhs);
    Par.ldx = Par.n;
    //dGenXtrue(n, nrhs, xact, ldx);
    //dFillRHS(trans, nrhs, xact, ldx, &A, &B);

	printf("Other SuperLU variables\n");
    if ( !(Par.etree = intMalloc(Par.n)) ) ABORT("Malloc fails for etree[].");
    if ( !(Par.perm_r = intMalloc(Par.m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(Par.perm_c = intMalloc(Par.n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(Par.R = (double *) SUPERLU_MALLOC(Par.A.nrow * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(Par.C = (double *) SUPERLU_MALLOC(Par.A.ncol * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for C[].");
    if ( !(Par.ferr = (double *) SUPERLU_MALLOC(Par.nrhs * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(Par.berr = (double *) SUPERLU_MALLOC(Par.nrhs * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for berr[].");
   //ONLY PERFORM THE LU DECOMPOSITION 
    Par.B.ncol = Par.nrhs;  //Indicate not to solve the system 
	printf("Finally LU decomposition\n");
    StatInit(&Par.stat);
    dgssvx(&(Par.options), &(Par.A), Par.perm_c, Par.perm_r, Par.etree, Par.equed, Par.R, Par.C,
           &(Par.L), &(Par.U), Par.work, Par.lwork, &(Par.B), &(Par.X), &(Par.rpg), &(Par.rcond), Par.ferr, Par.berr,
           &(Par.mem_usage), &(Par.stat), &(Par.info));

    StatFree(&(Par.stat));
    Par.options.Fact = FACTORED; // Indicate the factored form of A is supplied.
    Par.B.ncol = Par.nrhs;  // Set the number of right-hand side
	printf("Freeing stuff\n");
	delete [] xa;
	delete [] asub;
	delete [] a;
	printf("ix_size %d alpha_size %d\n",ix.size(),alphax.size());
	alphax=alphatempx;
	alphay=alphatempy;
	printf("New alphax_size %d\n",alphax.size());
	printf("Size of tangent vectors %d %d\n",tanx.size(),tany.size());
	}


//FOR RECTANGULAR GEOMETRY
void preparation_of_matrix_superlu_rectangular(long int Nx, long int Ny, double dx, double dy,
struct superlu_params &Par,
 vector <double> &alphax, vector <double> &alphay,
vector <int> &iindex,
vector <int> &bindex,
//ANGLE BETWEEN THE X AND Y AXIS AND THE WALL NEAR THE CUT CELL
vector <double> &tanx, vector <double> &tany)
	{
double dxdy=dx*dy;
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Enter main()");
#endif
//     Defaults 
    Par.lwork = 0;
    Par.nrhs  = 1;
    Par.equil = YES;
    //Par.IterRefine = ::DOUBLE;	
    Par.u     = 0.1;
    Par.trans = NOTRANS;
    set_default_options(&Par.options);

    Par.options.Equil = Par.equil;
    Par.options.DiagPivotThresh = Par.u;
   // Par.options.IterRefine = SLU_DOUBLE;
    Par.options.Trans = Par.trans;
    Par.options.ColPerm=MMD_AT_PLUS_A;
    Par.options.SymmetricMode = YES;


    if ( Par.lwork > 0 ) {
	Par.work = SUPERLU_MALLOC(Par.lwork);
	if ( !Par.work ) {
	    ABORT("DLINSOLX: cannot allocate work[]");
	}
    }
printf("Creating arrays\n");
vector <double> Aa;
vector <int> Arows;
vector <int> Acolumn;
vector <int> incolumn;
vector <int> ix;
vector <int> iy;
//TEMPORARY VARIABLE FOR THE CUT-CELL INFORMATION
vector <double> alphatempx;
vector <double> alphatempy;
double Dx=dx*(Nx-1)/2.0;
double Dy=dy*(Ny-1)/2.0;
//HORIZONTAL INDEXING
Acolumn.push_back(0);

for(int i=1;i<Nx-1;++i)
	{
	int num=0;
	for(int j=1;j<Ny-1;++j)
		{
		double xx=(dx*i-Dx);
		double yy=(dy*j-Dy);
		double x=fabs(xx);
		double y=fabs(yy);
		//TAKING ONLY THE POINT INSIDE THE PIPE

			//ADDING INDICES TO TO VECTORS
			ix.push_back(i);
			iy.push_back(j);
			iindex.push_back(i*Ny+j);
			num++;
			//IF SHIFT ALONG X OR ALONG Y PER ONE CELL GETS US OUTSIDE THE PIPE
			{
			alphax.push_back(1.0);
			alphay.push_back(1.0);
			}
			
		}
	if(num>0)
		{
		incolumn.push_back(num);
		}
	}

int totsize=ix.size();
int ixmin=ix.at(0);
double dxdx=dx*dx;
double dydy=dy*dy;
//CONSTRUCTING MATRIX X
//vector <double> matrix(totsize*totsize,0.0);
for(int i=0;i<totsize;++i)
	{
	int index=ix[i]-ixmin;
	//UPPER DIAGONAL
	if(index>0)
		{
		int deltam=(incolumn.at(index)+incolumn.at(index-1))/2;
		//IF SHIFT BY deltam PRESERVES THE ROW NUMBER
		if(i-deltam>=0 && iy[i]==iy.at(i-deltam))
			{
			Aa.push_back(-EPS0*dxdy/dxdx/*/(1.0+alphax[i])*/);
			Arows.push_back(i-deltam);
			//matrix.at(i*totsize+i-deltam)=Aa.back()/EPS0;
			}
		}
	//UP DIAGONAL//IF SHIFT BY ONE ELEMENT PRESERVES THE COLUMN
	if(i>0 && ix[i]==ix.at(i-1))
		{
		//printf("Before up\n");
		Aa.push_back(-EPS0*dxdy/dydy/*/(1.0+alphay[i])*/);
		Arows.push_back(i-1);
		//matrix.at(i*totsize+i-1)=Aa.back()/EPS0/dydy;
		}
	//MAIN DIAGONAL
		//Aa.push_back(2.0*EPS0*(1.0/dxdx/alphax[i]+1.0/dydy/alphay[i]));
		//Aa.push_back(EPS0*(2.0/dxdx/alphax[i]+2.0/dydy/alphay[i]) );
		Aa.push_back(EPS0*dxdy*((1.0+alphax[i])/alphax[i]/dxdx+(1.0+alphay[i])/alphay[i]/dydy) );
		Arows.push_back(i);
		//matrix.at(i*totsize+i)=Aa.back()/EPS0;
	//LOW DIAGONAL//IF SHIFT BY ONE ELEMENT PRESERVES THE COLUMN
	if(i+1<totsize && ix[i]==ix.at(i+1))
		{
		Aa.push_back(-EPS0*dxdy/dydy/*/(1.0+alphay[i])*/);
		Arows.push_back(i+1);
		//matrix.at(i*totsize+i+1)=Aa.back()/EPS0;
		}
	//LOWER DIAGONAL
	if(index+1<incolumn.size())
		{
		int deltap=(incolumn.at(index)+incolumn.at(index+1))/2;
		//IF SHIFT BY deltap PRESERVES THE ROW NUMBER
		if(i+deltap<totsize && iy[i]==iy.at(i+deltap))
			{
			Aa.push_back(-EPS0*dxdy/dxdx/*/(1.0+alphax[i])*/);
			Arows.push_back(i+deltap);
			//matrix.at(i*totsize+i+deltap)=Aa.back()/EPS0;
			}
		}
	Acolumn.push_back(Aa.size());
	}
	/*if(totsize<=10)
	{
	FILE *fd=fopen("matrice.dat","w");
	for(int i=0;i<totsize;i++)
		{
		for(int j=0;j<totsize;j++)
			{
			fprintf(fd,"%e ",matrix[i*totsize+j]);
			}
		fprintf(fd, "\n");
		}
	fclose(fd);
	}*/
	printf("Rewriting vectors to dynamic arrays\n");
	double *a=new double[Aa.size()];
	memcpy(a,&Aa[0],Aa.size()*sizeof(double));
	int *asub=new int[Arows.size()];
	memcpy(asub,&Arows[0],Arows.size()*sizeof(int));
	int *xa=new int[Acolumn.size()];
	memcpy(xa,&Acolumn[0],Acolumn.size()*sizeof(int));

Par.nnz=Aa.size();
Par.rhsb=new double[totsize];
Par.rhsx=new double[totsize];
for(int j=0;j<totsize;++j)
	{
	Par.rhsb[j]=1.0;
	Par.rhsx[j]=1.0;
	}
    Par.m=totsize;
    Par.n=totsize;
	printf("Creating matrices in SuperLU format\n");
    dCreate_CompCol_Matrix(&(Par.A), Par.m, Par.n, Par.nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&(Par.B), Par.m, Par.nrhs, Par.rhsb, Par.m, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&(Par.X), Par.m, Par.nrhs, Par.rhsx, Par.m, SLU_DN, SLU_D, SLU_GE);
    Par.xact = doubleMalloc(Par.n * Par.nrhs);
    Par.ldx = Par.n;
    //dGenXtrue(n, nrhs, xact, ldx);
    //dFillRHS(trans, nrhs, xact, ldx, &A, &B);

	printf("Other SuperLU variables\n");
    if ( !(Par.etree = intMalloc(Par.n)) ) ABORT("Malloc fails for etree[].");
    if ( !(Par.perm_r = intMalloc(Par.m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(Par.perm_c = intMalloc(Par.n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(Par.R = (double *) SUPERLU_MALLOC(Par.A.nrow * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(Par.C = (double *) SUPERLU_MALLOC(Par.A.ncol * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for C[].");
    if ( !(Par.ferr = (double *) SUPERLU_MALLOC(Par.nrhs * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(Par.berr = (double *) SUPERLU_MALLOC(Par.nrhs * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for berr[].");
   //ONLY PERFORM THE LU DECOMPOSITION 
    Par.B.ncol = Par.nrhs;  //Indicate not to solve the system 
	printf("Finally LU decomposition\n");
    StatInit(&Par.stat);
    dgssvx(&(Par.options), &(Par.A), Par.perm_c, Par.perm_r, Par.etree, Par.equed, Par.R, Par.C,
           &(Par.L), &(Par.U), Par.work, Par.lwork, &(Par.B), &(Par.X), &(Par.rpg), &(Par.rcond), Par.ferr, Par.berr,
           &(Par.mem_usage), &(Par.stat), &(Par.info));

    StatFree(&(Par.stat));
    Par.options.Fact = FACTORED; // Indicate the factored form of A is supplied.
    Par.B.ncol = Par.nrhs;  // Set the number of right-hand side
	printf("Freeing stuff\n");
	delete [] xa;
	delete [] asub;
	delete [] a;
	printf("ix_size %d alpha_size %d\n",ix.size(),alphax.size());
	alphax=alphatempx;
	alphay=alphatempy;
	printf("New alphax_size %d\n",alphax.size());
	printf("Size of tangent vectors %d %d\n",tanx.size(),tany.size());
	}

void calc_potential_superlu(double *rho, struct superlu_params *Par, long int Nx, long int Ny, double dx, double dy,
vector <int> iindex)
	{
	int totsize=iindex.size();
	for(int i=0;i<totsize;++i)
		{
		Par->rhsb[i]=rho[iindex[i]];
		Par->rhsx[i]=rho[iindex[i]];
		}
	dCreate_Dense_Matrix(&Par->B, Par->m, Par->nrhs, Par->rhsb, Par->m, SLU_DN, SLU_D, SLU_GE);
	dCreate_Dense_Matrix(&Par->X, Par->m, Par->nrhs, Par->rhsx, Par->m, SLU_DN, SLU_D, SLU_GE);
    StatInit(&(Par->stat));
    dgssvx(&(Par->options), &(Par->A), Par->perm_c, Par->perm_r, Par->etree, Par->equed, Par->R, Par->C,
           &(Par->L), &(Par->U), Par->work, Par->lwork, &(Par->B), &(Par->X), &(Par->rpg), &(Par->rcond), Par->ferr, Par->berr,
           &(Par->mem_usage), &(Par->stat), &(Par->info));
    StatFree(&(Par->stat));
	double *data=(double*) ((DNformat*) Par->X.Store)->nzval;
	set_zero(rho,Nx*Ny);
	for(int i=0;i<totsize;++i)
		{
		rho[iindex[i]]=data[i];
		}
	//free(data);
	}


//THIS FUNCTION USES LU DECOMPOSITION TO SOLVE THE ELECTRIC POTENTIAL IN RECTANGULAR GRID
//A_banded IS AN LU DECOMPOSED MATRIX REPRESENTATION FROM LAPACK IN BANDED FORM
void calc_potential_cut_cell(double *rho, double *A_banded, long int *IPIV, long int Nx, long int Ny, double dx, double dy,
vector <int> ix, vector <int> iy)
	{
	//NUMBER OF INTERNAL POINTS
	long int totsize=ix.size();
	//NUMBER OF ELEMENTS IN COLUMN OF A_BAND
	long int bandwidth=2*Ny+1+Ny;
	//CREATING RIGHT HAND SIDE
	double *rhotemp=(double*)malloc(totsize*sizeof(double));
	for(int i=0;i<totsize;++i)
		{
		rhotemp[i]=rho[iy[i]+Ny*ix[i]];//*area_cut[i];
		}


	long int nrhs=1;
	long int INFO=-1;
	dgbtrs_("N",&totsize,&Ny,&Ny,&nrhs,A_banded,&bandwidth,IPIV,rhotemp,&totsize,&INFO);
	for(int i=0;i<totsize;++i)
		{
		rho[ix[i]*Ny+iy[i]]=rhotemp[i];
		}
	free(rhotemp);
	};


void calc_potential_1D(double *rhox, double *rhox_inv, double *cosnx, int Nx, double dx, fftw_plan forward, fftw_plan backward)
	{
	//EXECUTING FFTW, TRANSFORM IS WRITTEN TO rho_inv
	fftw_execute(forward);
	//CALCULATING SOME INTERMEDIATE VARIABLES
	double dxdx=dx*dx;
	double Norma=1.0/(Nx+1)/2.0/EPS0;
	//GOING THROUGH THE TRANSFORM
	for(int i=0;i<Nx;++i)
		{		
		rhox_inv[i]=-Norma*rhox_inv[i]/((2.0-2.0*cosnx[i])/(dxdx));
		}
	//EXECUTING INVERSE SINE FFT TRANSFORM, WRITTEN TO rho
	fftw_execute(backward);
	};

//ELECTRIC FIELD SOLUTION. GRADIENT OF CALCULATED POTENTIAL
//ctomdt IS 1.0 IF YOU WANT TO CALCULATE THE FIELD, IF DIRECTLY THE ACCELERATION THEN ELECTRON CHARGE TO MASS
void calc_field(double *phi, double *Ex, double *Ey, int Nx, int Ny, double dx, double dy, double ctomdt)
	{
	int flagx=0;
	int flagy=0;
	int globind=0;
	for(int i=0;i<Nx;++i)
		{
		for(int j=0;j<Ny;++j)
			{
			//GLOBAL INDEX
			globind=i+Nx*j;
			//DERIVATIVE FLAG FOR X DIRECTION
			if(i==0)
				{
				flagx=-1;
				}
			else if(i==Nx-1)
				{
				flagx=1;
				}	
			else
				{
				flagx=0;
				}
			//DERIVATIVE FLAG FOR Y DIRECTION
			if(j==0)
				{
				flagy=-1;
				}
			else if(j==Ny-1)
				{
				flagy=1;
				}	
			else
				{
				flagy=0;
				}
			//TAKING DERIVATIVE OF POTENTIAL
			//Ex[globind]=ctomdt*derivative(phi[globind-1-flagx],phi[globind-flagx],phi[globind+1-flagx],dx,flagx);
			//Ey[globind]=ctomdt*derivative(phi[globind-(1+flagy)*Nx],phi[globind-flagy*Nx],phi[globind+(1-flagy)*Nx],dy,flagy);
			Ex[globind]=ctomdt*derivative(phi[globind-(1+flagx)*Ny],phi[globind-flagx*Ny],phi[globind+(1-flagx)*Ny],dx,flagx);
			Ey[globind]=ctomdt*derivative(phi[globind-1-flagy],phi[globind-flagy],phi[globind+1-flagy],dy,flagy);
			}
		}
	};
//ELECTRIC FIELD SOLUTION. GRADIENT OF CALCULATED POTENTIAL. TAKES CUT CELL INFORMATION TO CALCULATE THE CORRECT FIELD AT THE BOUNDARY AND EXTRAPOLATE IT ONE CELL AWAY FROM THE BOUNDARY

void calc_field_cut_cell(
//FIELD ARRAYS
double *phi, double *Ex, double *Ey,
//VECTOR OF INTERNAL GRID COORDINATES
vector <int> iindex,
//VECTOR OF POINTS NEAR THE BOUNDARY
vector <int> bindex,
//VECTORS OF RELATIVE CUT-CELL EDGES
vector <double> alphax, vector <double> alphay,
vector <double> tanx, vector <double> tany,
//SIZE AND STEP SIZE OF THE GRID
 int Nx, int Ny, double dx, double dy, 
double ctomdt)
	{
	double Dx=dx*(Nx-1)/2;
	double Dy=dy*(Ny-1)/2;
	int flagx=0;
	int flagy=0;
	int globind=0;
	int totsize=iindex.size();
	double dxdx=dx*dx;
	double dydy=dy*dy;
	double constx=ctomdt/2.0/dx;
	double consty=ctomdt/2.0/dy;
	for(int i=0;i<totsize;++i)
		{
		globind=iindex[i];
		//INSIDE OF PIPE
		if(phi[globind+Ny]!=0.0 && phi[globind-Ny]!=0.0)
			{
			Ex[globind]=constx*(phi[globind+Ny]-phi[globind-Ny]);
			}
		if(phi[globind+1]!=0.0 && phi[globind-1]!=0.0)
			{
			Ey[globind]=consty*(phi[globind+1]-phi[globind-1]);
			}
		}

	int bordersize=bindex.size();
	//FIRST TIME WE CALCULATE THE SIMPLEST DERIVATIVES
	for(int i=0;i<bordersize;++i)
		{
		globind=bindex[i];
		int by=globind%Ny;
		int bx=globind/Ny;
		//AT THE RIGHT SIDE OF A PIPE
		if(phi[globind+Ny]==0.0)
			{
			double E1=(phi[globind]-phi[globind-Ny])/dx;
			double E2=-phi[globind]/alphax[i]/dx;
			if(alphax[i]==1.0)
				{			
				printf("%e %e %e\n",alphax[i],(dx*bx-Dx)/0.0609, (dy*by-Dy)/0.02425);
				sleep(1);
				}
			Ex[globind]=ctomdt*(E1+E2*alphax[i])/(1.0+alphax[i]);
			//FIELD AT THE BOUNDARY
			//double E0x=2.0*E2-Ex[globind];
			Ex[globind+Ny]=2.0*Ex[globind]-Ex[globind-Ny];
			}
		//AT THE LEFT SIDE OF A PIPE
		if(phi[globind-Ny]==0.0)
			{
			double E1=(phi[globind+Ny]-phi[globind])/dx;
			double E2=phi[globind]/alphax[i]/dx;
			Ex[globind]=ctomdt*(E1+E2*alphax[i])/(1.0+alphax[i]);
			Ex[globind-Ny]=2.0*Ex[globind]-Ex[globind+Ny];
			}
		//AT THE UPPER SIDE OF A PIPE
		if(phi[globind+1]==0.0)
			{
			double E1=(phi[globind]-phi[globind-1])/dy;
			double E2=-phi[globind]/alphay[i]/dy;
			Ey[globind]=ctomdt*(E1+E2*alphay[i])/(1.0+alphay[i]);
			Ey[globind+1]=2.0*Ey[globind]-Ey[globind-1];
			}
		//AT THE LOWER SIDE OF A PIPE
		if(phi[globind-1]==0.0)
			{
			double E1=(phi[globind+1]-phi[globind])/dy;
			double E2=phi[globind]/alphay[i]/dy;
			Ey[globind]=ctomdt*(E1+E2*alphay[i])/(1.0+alphay[i]);
			Ey[globind-1]=2.0*Ey[globind]-Ey[globind+1];
			}
		}
	//SECOND TIME WE INTERPOLATE TRANSVERSE DERIVATIVES (NOT USING ANGLES AT A TIME)
	for(int i=0;i<bordersize;++i)
		{
		globind=bindex[i];
		//int by=globind%Ny;
		//int bx=globind/Ny;
		if(phi[globind+Ny]==0.0 && Ey[globind+Ny]==0.0)
			{
			//WE INTERPOLATE OUTSIDE THE GRID
			Ey[globind+Ny]=2.0*Ey[globind]-Ey[globind-Ny];
			}
		if(phi[globind+1]==0.0 && Ex[globind+1]==0.0)
			{

			Ex[globind+1]=2.0*Ex[globind]-Ex[globind-1];
			}
		if(phi[globind-Ny]==0.0 && Ey[globind-Ny]==0.0)
			{
			//WE CALCULATE THE Ex ON THE BORDER
			Ey[globind-Ny]=2.0*Ey[globind]-Ey[globind+Ny];
			}
		if(phi[globind-1]==0.0 && Ex[globind-1]==0.0)
			{

			Ex[globind-1]=2.0*Ex[globind]-Ex[globind+1];
			}
		}

	};
/*void interpolate_stuff(//FIELD ARRAYS
double *Ex, double *Ey,
vector <int> bx, vector<int> by,
//VECTORS OF RELATIVE CUT-CELL EDGES
vector <double> alphax, vector <double> alphay,
//SIZE AND STEP SIZE OF THE GRID
 int Nx, int Ny)
	{
	int bordersize=bx.size();
	for(int i=0;i<bordersize;i++)
		{
		
		}
	}*/


//SOLVER OF A POTENTIAL IN PURE RECTANGULAR 3D GEOMETRY
void calc_potential_3d(double *rho, double *rho_inv, double *cosnx, double *cosny, double *cosnz, int Nx, int Ny, int Nz, double dx, double dy, double dz, fftw_plan forward, fftw_plan backward)
	{
	fftw_execute(forward);
	double dxdx=dx*dx;
	double dydy=dy*dy;
	double dzdz=dz*dz;
	int slice=Nx*Ny;
	double Norma=1.0/EPS0/(Nx+1)/(Ny+1)/(Nz+1)/8.0;
	for(int i=0;i<Nx;++i)
		{
		for(int j=0;j<Ny;++j)
			{
			for(int k=0;k<Nz;++k)
				{
				//GRID DIMENSIONS ARE NOT EQUAL IN X AND Y
				rho_inv[i+Nx*j+slice*k]=-Norma*rho_inv[i+Nx*j+slice*k]/((2.0-2.0*cosnx[i])/(dxdx)+(2.0-2.0*cosny[j])/(dydy)+(2.0-2.0*cosnz[k])/(dzdz));
				}
			}
		}
	fftw_execute(backward);
	};
//ELECTRIC FIELD SOLUTION
void calc_field_3d(double *phi, double *Ez, int Nx, int Ny, int Nz, double dx, double dy, double dz, double ctomdt)
	{
	int flagx=0;
	int flagy=0;
	int flagz=0;
	int globind=0;
	int NxNy=Nx*Ny;
	for(int i=0;i<Nx;++i)
		{
		for(int j=0;j<Ny;++j)
			{
			for(int k=0;k<Nz;++k)
				{
				//GLOBAL INDEX
				globind=i+Nx*j+NxNy*k;
				//DERIVATIVE FLAG FOR Y DIRECTION
				if(k==0)
					{
					flagz=-1;
					}
				else if(k==Nz-1)
					{
					flagz=1;
					}	
				else
					{
					flagz=0;
					}
				//TAKING DERIVATIVE OF POTENTIAL
				Ez[globind]=ctomdt*derivative(phi[globind-(1+flagz)*NxNy],phi[globind-flagz*NxNy],phi[globind+(1-flagz)*NxNy],dz,flagz);
				}
			}
		}
	};
//CALCULATE AVERAGE FIELD ACTING ON THE BEAM
void calc_Ez_aver(double *wake_z, double *Ez, double *rho_i, int Nx, int Ny, int Nz, double dx, double dy, double dz)
	{
	int NxNy=Nx*Ny;
	int globind=0;
	double sum=0.0;
	//THIS VALUE IS USED FOR NORMALIZATION
	for(int j=0;j<Ny;++j)
		{
		for(int k=0;k<Nx;++k)
			{
			sum=sum+rho_i[Nx*j+k];
			}
		}
	//CALCULATING FIELD
	for(int i=0;i<Nz;++i)
		{
		wake_z[i]=0.0;
		for(int j=0;j<Ny;++j)
			{
			for(int k=0;k<Nx;++k)
				{
				globind=k+Nx*j+NxNy*i;
				if(rho_i[Nx*j+k]!=0.0)
					{
					wake_z[i]+=rho_i[Nx*j+k]*Ez[globind];
					}
				}
			}
		wake_z[i]=wake_z[i]/sum;
		}
	
	}




