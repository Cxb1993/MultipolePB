#pragma once

//
#include <complex>

typedef std::complex<double> complex_t;



//
// Bessel functions
//

double mathModSphBesselI(int n, double x);
double mathModSphBesselK(int n, double x);
double mathAdaptModSphBesselI(int n, double x);
double mathAdaptModSphBesselK(int n, double x);

//
// Legrendre polynomials
//

double mathAssocLegendreP(int l, int m, double x);

double mathCondonShortleyPhase(int m);


//
// Spherical Harmonics 
//
complex_t mathSphHarmonic(int n, int m, double theta, double phi);


//
// Gamma func.
//
double mathGamma(double x);


//
//
//
double mathFactorial(int n);
double mathDoubleFactorial(int n);



//
// Constants
//




