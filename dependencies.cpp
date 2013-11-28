#include <cmath>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fftw3.h>

using namespace std;

//SECONDARY EMISSION FOR A GIVEN IMPACT ENERGY
double sey(double Emax, double dmax, double s, double Eimp)
	{
	double x=Eimp/Emax;

	return dmax*s*x/(s-1.0+pow(x,s));

	};
//REDIFFUSED ELECTRONS
double rediffused(double Er, double rmax, double Eimp)
	{
	double x=Eimp/Er;
	if(x<5)
	{
	return rmax*(1.0-exp(-x));
	}
	else
	{
	return rmax;
	}
	};
//SECONDARY EMISSION CONSTANT AFTER MAX
double sey_const_max(double Emax, double dmax, double s, double Eimp)
	{
	double x=Eimp/Emax;
	if(x<1)
		{
		return dmax*s*x/(s-1.0+pow(x,s));
		}
	else
		{
		return dmax;
		}
	};

//REFLECTION COEFFICIENT FOR A GIVEN IMPACT ENERGY
double refl(double refl0, double E0, double Eimp)
	{
	double first=sqrt(Eimp+E0);
	double second=sqrt(Eimp);
	return refl0*pow((first-second)/(first+second),2);
	};
//EMISSION ANGLE WITH COSINE DISTRIBUTION
void set_angles(vector <double> &cos_sin, int Num)
	{
	cos_sin.assign(2*Num,0.0);
	for(int i=0;i<Num;i++)
		{
		double phi=asin(2.0*((double)i+0.5)/Num-1.0);
		cos_sin.at(i)=cos(phi);
		cos_sin.at(i+Num)=sin(phi);
		}
	}
//NUMBER OF ELECTRONS GENERATED IN SIMPLE EMISSION MODEL
int n_emitted(double impact_sey, gsl_rng *rr)
	{
	double ne=floor(impact_sey);
	if(gsl_rng_uniform(rr)<impact_sey-ne)
		{
		ne+=1.0;
		}
	return (int)ne;
	}

//FUNCTION GIVING TIME
string taketime()
	{
	char ctime[25];
	string s;
	time_t curtime;
	struct tm *loctime;
  /* Get the current time. */
	curtime = time (NULL);
  /* Convert it to local time representation. */
	loctime = localtime (&curtime);
	strftime (ctime, 25, "%F-%H%M%S", loctime);
	return s.append(ctime);
	}
