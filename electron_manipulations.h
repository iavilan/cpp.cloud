
struct electron {
//PHYSICAL PARAMETERS
double x;
double y;
double Vx;
double Vy;
double Vz;
double lambda_mac;
//GRID PARAMETERS
//GRID WEIGHTS OF PARTICLES
/*double downleft;
double upleft;
double downright;
double upright;
//INDICES OF A LOWR GRID POINT CLOSEST TO A PARTICLE
int low_y;
int low_x;*/
};

struct weight;

//THIS FILE CONTAINS DEFINITIONS OF FUNCTIONS THAT ARE USED FOR ELECTRON MANIPULATIONS
using namespace std;
//GETTING THE POSITION OF A POINT WHERE ELECTRON CROSSES THE WALL AND NORM TO THE SURFACE
//xc,yc - coordinates of a cross point
///ex, ey - is the unit norm vector
void cross_point(int geoflag, double Rx, double Ry, double x, double y, double Vx, double Vy, double dt, double &xc, double &yc, double &ex, double &ey);
//GET NUMBER OF ELECTRONS TO PRODUCE
int number_to_produce(double impactsey);

//KICK ELECTRON VELOCITY
void accelerate_electron(double *force_dt, double &Vx, double &Vy, double dt);


//MOVE ELECTRON WITH A GIVEN VELOCITY AND TIME STEP
void shift_electron(double &x, double &y, double Vx, double Vy, double dt);


//GENERATE VELOCITY OF AN ELECTRON (TAKES SEY PARAMETERS AND RANDOM NUMBWER GENERATOR)
double gen_velocity(double sigmav, gsl_rng *rr);
//GENERATE VELOCITY OF AN ELECTRON (WITH ZERO PROBABILITY OF ZERO VELOCITY)
double gen_velocity_peak(double sigV, gsl_rng *rr);
//GENERATE ANGLE OF AN EMITTED ELECTRON (COS THETHA DISTRIBUTION)
double gen_angle(gsl_rng *rr);


//CREATE INITIAL UNIFORM ELLIPTICAL DISTRIBUTION
void uniform_cloud_ellipse(double Rx, double Ry, vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, vector <double> &lambda_mac, double lambda_mac0, int NUM, gsl_rng *rr);
void uniform_cloud_ellipse(double Rx, double Ry, vector <electron> &ELECTRONS, double lambda_mac0, int NUM, gsl_rng *rr);
//ELLIPTICAL DISTRIBUTION
void gaussian_cloud_ellipse(double Rx, double Ry, double DX, vector <electron> &ELECTRONS, double lambda_mac0, int NUM, gsl_rng *rr);

//CREATE AN ELECTRON CLOUD RING
void uniform_cloud_ring(double Rmin, double Rmax, vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, vector <double> &lambda_mac, double lambda_mac0, int NUM, gsl_rng *rr);
void uniform_cloud_ring(double Rmin, double Rmax, vector <electron> &ELECTRONS, double lambda_mac0, int NUM, gsl_rng *rr);

//CREATE RECTANGULAR DISTRIBUTION
void uniform_cloud_rect(double Rx, double Ry, vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, vector <double> &lambda_mac, double lambda_mac0, int NUM, gsl_rng *rr);
void uniform_cloud_rect(double Rx, double Ry, vector <electron> &ELECTRONS, double lambda_mac0, int NUM, gsl_rng *rr);


//MERGE(DELETE SOME) PARTICLES WHEN THEIR TOTAL NUMBER EXCEEDS THE LIMIT Nmax
void simple_merge(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, vector <double> &lambda_mac, int MAXINUM, gsl_rng *rr);
void simple_merge(vector <electron> &ELECTRONS, int MAXINUM, gsl_rng *rr);
//SIMPLE SPLITTING OF PARTICLES IN CASE OF SYMMETRY
void simple_split(vector <electron> &ELECTRONS, int MINNUM, gsl_rng *rr);

//FUNCTION THAT SIMPLY REMOVES ALL THE ELECTRONS 
void remove_all_out_rect(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, vector <double> &weights, double Rx, double Ry);
void remove_all_out_rect(vector <electron> &ELECTRONS, double Rx, double Ry);

//FUNCTION THAT REFLECTS  ALL THE ELECTRONS
void reflect_all_out_rect(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, double Rx, double Ry, double absorption);
void reflect_all_out_rect(vector <electron> &ELECTRONS, double Rx, double Ry, double absorption);



//FUNCTION THAT REFLECTS AN ELECTRON WITH i NUMBER
void reflect(
	//VECTORS TO ADD IN
	vector <double> &xep, vector <double> &yep, vector <double> &Vxep, vector <double> &Vyep, vector <double> &Vzep, vector 	<double> &lambda_macp,
	//VECTORS TO TAKE FROM
	double xeo, double yeo, double Vxeo, double Vyeo, double lambda_maco,
	//CROSSING POINT
double xc, double yc, double ex, double ey, double dt,
	//PIPE PARAMETERS
	double Rx, double Ry, int geoflag);
void reflect(
	//VECTORS TO ADD IN
	vector <electron> &ELECTRONS,
	//VECTORS TO TAKE FROM
	double xeo, double yeo, double Vxeo, double Vyeo, double Vzeo, double lambda_maco,
	//CROSSING POINT
	double xc, double yc, double ex, double ey, double dt,
	//PIPE PARAMETERS
	double Rx, double Ry, int geoflag);

//IMPLEMENTATION OF WALL EFFECTS IN RECTANGULAR GEOMETRY
void wall_effects(
	vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, vector <double> &lambda_mac, 
	double Rx, double Ry, double dmax, double Emax, double refl0, double Erefl, double sigma_v, 
	vector <double> cos_sin, double dt, 
	int geoflag, gsl_rng *rr);

void wall_effects(
	vector <electron> &ELECTRONS, 
	double Rx, double Ry, double dmax, double Emax, double refl0, double Erefl, double sigma_v, 
	vector <double> cos_sin, double dt, 
	int geoflag, gsl_rng *rr);
//WALL EFFECTS WITH PEAK VELOCITY
void wall_effects_peak(vector <electron> &ELECTRONS, 
double Rx, double Ry, 
double dmax, double Emax, double refl0, double Erefl, double sigma_v,  double redif,
vector <double> cos_sin, double dt, int geoflag, gsl_rng *rr,
double &energy_deposition,
double sss);


//KICK ALL THE ELECTRONS IN THE CLOUD
void accelerate_all(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, 
	double *Ex2D, double *Ey2D, 
	int Nx, int Ny, 
	double dx, double dy);
void accelerate_all(vector <electron> &ELECTRONS, 
	double *Ex2D, double *Ey2D, 
	int Nx, int Ny, 
	double dx, double dy);

//SHIFT ALL ELECTRONS
void shift_all(vector <double> &xe, vector <double> &ye, 
	vector <double> &Vxe, vector <double> &Vye, 
	double dt);
void shift_all(vector <electron> &ELECTRONS, 
	double dt);

//MOVE ALL THE ELECTRONS (ACCELERATES AND SHIFTS)
void move_all(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, double *Ex2D, double *Ey2D, 	int Nx, int Ny, double dx, double dy,double dt);
void move_all(vector <electron> &ELECTRONS, double *Ex2D, double *Ey2D, int Nx, int Ny, double dx, double dy,double dt);


//MOVE ALL THE ELECTRONS IN MAGNETIC FIELD (ACCELERATES / ROTATES / ACCELERATES / SHIFTS)
void move_all_in_B(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, 		double *Ex2D, double *Ey2D, 
	double cosi, double sini, 
	double Rx, double Ry, 
	int Nx, int Ny, 
	double dx, double dy,double dt, 
	int geoflag);
void move_all_in_B(vector <electron> &ELECTRONS, 	
	double *Ex2D, double *Ey2D, 
	double cosi, double sini, 
	double Rx, double Ry, 
	int Nx, int Ny, 
	double dx, double dy,double dt, 
	int geoflag);


//MOVE ALL ELECTRONS IN MAGNETIC FIELD (CHANGES VELOCITY AND COORDINATE SIMULATNEOUSLY) USING THE PARTICLE WEIGHTS FROM INTERPOLATION STEP
void move_all_in_B(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, vector <double> &Vze, 
	double *Ex2D, double *Ey2D,
	//PARAMETERS OF A PARTICLE RELATIVE TO THE GRID
	vector <double> &downleft,vector <double> &upleft,vector <double> &downright,vector <double> &upright,vector <int> &low_x,vector <int> &low_y,
	double cosi, double sini, 
	double Rx, double Ry, 
	int Nx, int Ny, 
	double dx, double dy, double dt, 
	int geoflag); //geoflag 0 - rect 1 - circle
void move_all_in_B(vector <electron> &ELECTRONS, 
	double *Ex2D, double *Ey2D,
	//PARAMETERS OF A PARTICLE RELATIVE TO THE GRID
	//vector <double> downleft,vector <double> upleft,vector <double> downright,vector <double> upright,vector <int> low_x,vector <int> low_y,
	vector <weight> WEIGHTS,
	double cosi, double sini, 
	double Rx, double Ry, 
	int Nx, int Ny, 
	double dx, double dy, double dt, 
	int geoflag); //geoflag 0 - rect 1 - circle


//ASSUME THAT ELECTRONS ARE MOVING ONLY VERTICALLY
void move_all_in_rigid_B(vector <double> &xe, vector <double> &ye, vector <double> &Vxe, vector <double> &Vye, double *Ey2D, double Rx, double Ry, int Nx, int Ny, double dx, double dy, double dt, int geoflag);
//MOVE ALL ELECTRONS IN MAGNETIC FIELD (CHANGES VELOCITY AND COORDINATE SIMULATNEOUSLY) USING STRUCTURES FOR ELECTRONS // ELECTROBNS ARE STRUCTURES
void move_all_in_rigid_B(vector <electron> &ELECTRONS, 
double *Ey2D,
vector <weight> WEIGHTS,
double Rx, double Ry, int Nx, int Ny, double dx, double dy, double dt, int geoflag);

void check_x_zero(vector <electron> ELECTRONS);

