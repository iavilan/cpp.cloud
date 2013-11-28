//THIS FILE CONTAINS DEFINITIONS OF OPERATIONS THAT CAN BE PERFORMED WITH GRID
struct electron; //DEFINED IN ELECTRONS_MANIPULATIONS.H
struct weight 
{
//GRID WEIGHTS OF PARTICLES
double downleft;
double upleft;
double downright;
double upright;
//INDICES OF A LOWR GRID POINT CLOSEST TO A PARTICLE
int low_index;
};
using namespace std;
//DERIVATIVE USING THREE POINTS
double derivative(double V1, double V2, double V3, double dL, int num);
//GET VALUE FROM POINT x,y
double get_value(double *rho, int Nx, int Ny, double dx, double dy, double x, double y);
//GET VALUE FROM POINT x,y FOR GIVEN WEIGHT
double get_value(double *rho, int Nx, int Ny, double dx, double dy, double dl, double ul, double dr, double du, int low_x, int low_y);
double get_value(double *rho, int Ny, double dl, double ul, double dr, double ur, int low_index);
//SET VALUE OF PARTICLE WITH x,y AND lambda CHARGE LINE DENSITY
void add_value(double *rho, int Nx, int Ny, double dx, double dy, double x, double y, double lambda);
//SET VALUE OF PARTICLE WITH x,y AND lambda CHARGE LINE DENSITY AND SAVING THE WEIGHT
void add_value(double *rho, int Nx, int Ny, double dx, double dy, double lambda, double &dl, double &ul, double &dr, double &ur, int &low_x, int &low_y);
//SET VALUE OF PARTICLE WITH x,y AND lambda CHARGE LINE DENSITY AND SAVING THE WEIGHT
void add_value(double *rho, int Nx, int Ny, double dx, double dy, double x, double y, double lambda, double &dl, double &ul, double &dr, double &ur, int &low_index);
//SUM OF THE GRIDS
void sum_grids(double *rho1,double *rho, int NxNy);
//MULTIPLY THE GRID BY NUMBER
double *grid_by_number(double *rho, double number, int NxNy);
//SET ELLIPSE
void set_ellipse(double *rho, int Nx, int Ny, double dx, double dy, double Ax, double Ay, double lambda0);
//SETTING GRID VALUES TO ZERO
void set_zero(double *rho, int NxNy);
//SETTING GRID VALUES TO ZERO
void set_zero(long double *rho, int NxNy);


//SIMPLE LINEAR INTERPOLATION OF ELECTRON VECTORS ON THE GRID
void simple_interpolate(vector <double> xe, vector <double> ye, vector <double> lambda_mac, 
	double *rho, int Nx, int Ny, double dx, double dy);
//SIMPLE LINEAR INTERPOLATION OF ELECTRON VECTORS ON THE GRID USING VECTOR OF ELECTRON STRUCTURES
void simple_interpolate(vector <electron> ELECTRONS, 
	double *rho, int Nx, int Ny, double dx, double dy);
//SIMPLE LINEAR INTERPOLATION OF ELECTRON VECTORS ON THE GRID STORING DATA FOR THE KICK
void simple_interpolate(vector <double> xe, vector <double> ye, vector <double> lambda_mac, 
	double *rho, int Nx, int Ny, double dx, double dy,
	vector <double> &downleft,vector <double> &upleft,vector <double> &downright,vector <double> &upright,vector <int> &low_x,vector <int> &low_y);

void simple_interpolate(vector <electron> &ELECTRONS, 
	double *rho, int Nx, int Ny, double dx, double dy,
	vector <weight> &WEIGHTS);

//RADIALLY SYMMETRIC INTERPOLATION
void radial_interpolate(vector <double> xe, vector <double> ye, vector <double> lambda_mac, double *rho_r, double Rr, int Nr, double dr);
//RADIALLY SYMMETRIC INTERPOLATION USING ELECTRON STRUCTURE
void radial_interpolate(vector <electron> ELECTRONS, double *rho_r, double Rr, int Nr, double dr);


//SAVING RADIALLY SYMMETRIC DATA
void radial_save(FILE *fd, double *rho_r, double Rr, int Nr, double dr);

//SAVE DISTRIBUTION IN AS A FUNCTION OF X IN j-th SLICE, j<Ny;
void x_save(FILE *fd, double *rho_e, int Nx, int Ny, int j);

//SAVE DISTRIBUTION IN AS A FUNCTION OF Y IN j-th SLICE, j<Nx;
void y_save(FILE *fd, double *rho_e, int Nx, int Ny, int j);


