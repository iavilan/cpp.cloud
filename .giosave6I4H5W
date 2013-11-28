//THIS FILE CONTAINS FUNCTIONS USED FOR DIFFERENT BEAM OPERATIONS

using namespace std;

//THE FUNCTION USED FOR RIGID BEAM TO WEIGHT THE FIELD OF SLICES (LOCAL LINE DENSITY RELATIVE TO THE AVERAGE ONE)
vector <double> gaussian_weight(int Nslice,double dz, double sigmaz);

vector <double> parabolic_weight(int Nslice,double dz, double total);
//THIS FUNCTION CREATES A BEAM-PAERTICLE DISTRIBUTION CORRESPONDING TO EACH SLICE OF FUNCTION gaussian_weight
void populate(vector < vector <double> > &xx, vector < vector <double> > &yy, 
vector < vector <double> > &xxp, vector < vector <double> > &yyp,
double Kx, double Ky,
int Nslice, double dz, double sigmaz, double Ax, double Ay, gsl_rng *rr);
//CONSTANT FOCUSING TRANSFER
void cf_transfer(vector < vector <double> > &xx, vector < vector <double> > &yy, 
vector < vector <double> > &Vxx, vector < vector <double> > &Vyy,
double Kx, double Ky,
double cosx, double sinx,
double cosy, double siny,
vector <double> &x, vector <double> &y, vector <double> &Vx, vector <double> &Vy);
//KICK PARTICLE
void kick_particle(double dxs, double dys,double &Vx, double &Vy);
//FOR GIVEN 2D TRANSVERSE FIELD KICK THE PARTICLES
void kick_slice(vector <double> &x, vector <double> &y, vector <double> &Vx, vector <double> &Vy, 
double *Ex2D, double *Ey2D, 
int Nx, int Ny, double dx, double dy);
//TRANSVERSE KICK OF A BEAM BY THE PRECALCULATED FIELD
//FOR GIVEN 3D TRANSVERSE FIELD KICK THE BEAM PARTICLES
void transverse_beam_kick_3d(
//PARTICLE DATA
vector < vector <double> > &xx, vector < vector <double> > &yy, 
vector < vector <double> > &Vxx, vector < vector <double> > &Vyy,
//FIELD DATA
vector < vector <double> > Ex, vector < vector <double> > Ey,
int Nx, int Ny, double dx, double dy
);
//TUNESHIFT EVALUATION // iter - NUMBER OF ITERATIONS THAT THE BEAM WAS TRANSFERED THROUGH THE LATTICE // USES LAPACK
void cf_tune_shift_evaluate(
//INPUT PARAMETERS
vector <double> x, vector <double> y, vector <double> Vx, vector <double> Vy, 
double fQx, double fQy,
int iter,
//OUTPUT PARAMETERS
vector <double> &dQx, vector <double> &dQy
);

//TUNESHIFT EVALUATUION AND SAVING THE SLICE NUMBER
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
);
