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

//CALCULATING POTENTIAL USING THE SUPERLU ALGORITHMS
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




