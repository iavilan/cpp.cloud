//THIS FILE CONTAINS DEFINITIONS OF FUNCTIONS DESCRIBING THE PHYSICAL PROPERTIES OF THE WALL
using namespace std;

//SECONDARY EMISSION FOR A GIVEN IMPACT ENERGY
double sey(double Emax, double dmax, double s, double Eimp);
//REDIFFUSED ELECTRONS
double rediffused(double Er, double rmax, double Eimp);
//AFTER MAXIMUM THIS IS CONSTANT
double sey_const_max(double Emax, double dmax, double s, double Eimp);
//REFLECTION COEFFICIENT FOR A GIVEN IMPACT ENERGY
double refl(double refl0, double E0, double Eimp);
//FUNCTION GENERATES A VECTOR WITH PARAMETERS OF UNIT VECTORS NEEDED TO EMIT ELECTRONS FROM THE WALL
//COSINE DISTRIBUTION IS CREATED
void set_angles(vector <double> &cos_sin, int Num);

//NUMBER OF EMITTED ELECTRONS FOR A SIMPLE EMISSION MODEL
//EXAMPLE: IF impact_sey=2.7, THEN IT IS ASSUMED THAT 2 ELECTRONS ARE EMITTED FOR SURE AND A THIRD ELECTRON IS EMITTED WITH THE PROBABILITY OF 0.7
double n_emitted(double impact_sey, gsl_rng *rr);
//FUNCTION GIVING TIME
string taketime();
