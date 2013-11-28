#include "slu_ddefs.h"
//THIS FILE CONTAINS DEFINITIONS OF OPERATIONS FOR THE FIELD EVALUATION
using namespace std;

//STRUCTURE FOR PASSING TO SUPERLU SOLVER
struct superlu_params {
    SuperMatrix A, L, U, B, X;
    yes_no_t       equil;
    double   *rhs;
    int      *perm_r; // row permutations from partial pivoting
    int      *perm_c; // column permutation vector
    double         *R, *C;
    int            *etree;
    void           *work;
    int            lwork, ldx;
    double         *rhsb, *rhsx, *xact;
    double         *ferr, *berr;
    double         u, rpg, rcond;
    int      nrhs, info, m, n, nnz, permc_spec;
    mem_usage_t    mem_usage;
    char           equed[1];
    superlu_options_t options;
    SuperLUStat_t stat;
    trans_t        trans;
};

//SETS BOUNDARY CONDITIONS
void set_boundary(double *rho, double Rx, double Ry, int Nx, int Ny, double dx, double dy);

//CALCULATING BOUNDARY VALUES FOR GIVEN GREEN FUNCTION
void set_boundary_condition(double *rho, double *greens, int *border_x, int *border_y, int Nx, int Ny, int geoflag);

//GETS POTENTIAL FROM A GIVEN DENSITY, GRID DIMENSIONS, AND GRID SIZE
void calc_potential(double *rho, double *rho_inv, double *cosnx, double *cosny, int Nx, int Ny, double dx, double dy, fftw_plan forward, fftw_plan backward);

//THIS FUNCTION USES LU DECOMPOSITION TO SOLVE THE ELECTRIC POTENTIAL FOR RECTANGULAR BOUNDARY
void calc_potential_lu(double *rho, int Nx, int Ny, double dx, double dy);

//PREPARATION OF LU DECOMPOSITION FOR A GIVEN 
void preparation_of_matrix(double Rx, double Ry,  long int Nx, long int Ny, double dx, double dy,
double *&A_lu, long int *&IPIV, 
//CUT CELL RELATIVE LENGTH
vector <double> &alphax, vector <double> &alphay,
//INDICES OF A RECTANGULAR GRID INSIDE THE DOMAIN
vector <int> &ix, vector <int> &iy,
//INDICES INSIDE THE DOMAIN NEAR THE BORDER (IF ONE CELL TO THE LEFT/RIGTH/UP/DOWN IS ALREADY OUTSIDE)
vector <int> &bx, vector <int> &by,
vector <double> &tanx, vector <double> &tany);

void modify(struct superlu_params *Par);


//PREPARATION OF MATRIX FOR ELLIPTICAL GEOMETRY
void preparation_of_matrix_superlu(double Rx, double Ry,  long int Nx, long int Ny, double dx, double dy,
struct superlu_params &Par,vector <double> &alphax, vector <double> &alphay,
vector <int> &iindex,
vector <int> &bindex,
vector <double> &tanx, vector <double> &tany);
//PREPARATION OF MATRIX FORRECTANGULAR GEOMETRY
void preparation_of_matrix_superlu_rectangular(long int Nx, long int Ny, double dx, double dy,
struct superlu_params &Par,
 vector <double> &alphax, vector <double> &alphay,
vector <int> &iindex,
vector <int> &bindex,
//ANGLE BETWEEN THE X AND Y AXIS AND THE WALL NEAR THE CUT CELL
vector <double> &tanx, vector <double> &tany);
/*THIS FUNCTION USES LU DECOMPOSITION TO SOLVE THE ELECTRIC POTENTIAL FOR ANY BOUNDARY GIVEN:
rho - density NxNy matrix
A_banded - a matrix in a banded form with dimensions 3*Ny+1 representing the LU decomposition (lapack)
IPIV - is calculated during the LU decomposition
ix,iy - are internal indices
*/
void calc_potential_cut_cell(double *rho, double *A_banded, long int *IPIV, long int Nx, long int Ny, double dx, double dy,
vector <int> ix, vector <int> iy);
//CALCULATING POTENTIAL USING THE SUPERLU ALGORITHMS
//void calc_potential_superlu(double *rho, struct superlu_params *Par, long int Nx, long int Ny, double dx, double dy,
//vector <int> ix, vector <int> iy);
void calc_potential_superlu(double *rho, struct superlu_params *Par, long int Nx, long int Ny, double dx, double dy,
vector <int> iindex);

//CALCULATE POTENTIAL IN 1D
void calc_potential_1D(double *rhox, double *rhox_inv, double *cosnx, int Nx, double dx, fftw_plan forward, fftw_plan backward);

//ELECTRIC FIELD SOLUTION (MAYBE DIRECTLY ACCELERATION MULTIPLIED BY TIMESTEP?)
void calc_field(double *phi, double *Ex, double *Ey, int Nx, int Ny, double dx, double dy,double ctomdt);

//CALCULATING FIELD WITH CUT-CELL KNOWLEDGE
void calc_field_cut_cell(
//FIELD ARRAYS
double *phi, double *Ex, double *Ey,
//VECTOR OF INTERNAL GRID COORDINATES
vector <int> iindex,
vector <int> bindex,
//VECTORS OF RELATIVE CUT-CELL EDGES
vector <double> alphax, vector <double> alphay,
vector <double> tanx, vector <double> tany,
//SIZE AND STEP SIZE OF THE GRID
 int Nx, int Ny, double dx, double dy, 
double ctomdt);

//SOLVER OF A POTENTIAL IN PURE RECTANGULAR 3D GEOMETRY
void calc_potential_3d(double *rho, double *rho_inv, double *cosnx, double *cosny, double *cosnz, int Nx, int Ny, int Nz, double dx, double dy, double dz, fftw_plan forward, fftw_plan backward);

//CALCULATING GRADIENT PART IN THE LONGITUDINAL PLANE
void calc_field_3d(double *phi, double *Ez, int Nx, int Ny, int Nz, double dx, double dy, double dz, double ctomdt);

//GENERATE GREEN FUNCTION FOR ROUND GEOMETRY AT ALL BOUNDARY POINTS
void gen_green_round(double *greens, int *border_x, int *border_y,double R, int Nx, int Ny, double dx, double dy);

//GREEN FUNCTION FOR ROUND BOUNDARY
double green_round(double x0, double y0, double x, double y, double R);

//SOLVE USING GREE
void calc_potential_green(double *rhoe,int Nx,int Ny,double dx, double dy, double R);

//AVERAGE FIELD ACTING ON EACH SLICE OF THE BEAM
void calc_Ez_aver(double *wake_z, double *Ez, double *rho_i, int Nx, int Ny, int Nz, double dx, double dy, double dz);


