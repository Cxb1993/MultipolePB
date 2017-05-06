
#include <cmath>

// we use implementations from BOOST
#include <boost/math/special_functions.hpp>

//
#include "mathfunc.h"


////////////////////////////////////////////////////////////////
// Bessel function
////////////////////////////////////////////////////////////////



double mathModSphBesselI(int n, double x) {
	double val = boost::math::cyl_bessel_i((double) n+0.5, x);
	val = sqrt(M_PI_2/x) * val;

	return val;
}

double mathModSphBesselK(int n, double x) {
	double val = boost::math::cyl_bessel_k((double) n+0.5, x);
	val = sqrt(M_PI_2/x) * val;

	return val;
}

double mathAdaptModSphBesselI(int n, double x) {
	double coef = mathDoubleFactorial(2*n+1) / pow(x,n);
	double in = mathModSphBesselI(n,x);
	double ihat = coef * in;
	return ihat;
}

double mathAdaptModSphBesselK(int n, double x) {
	double coef = exp(x) * pow(x,n+1) / mathDoubleFactorial(2*n-1);
	double kn = mathModSphBesselK(n,x);
	double khat = 2.0/M_PI * coef * kn;
	return khat;
}


////////////////////////////////////////////////////////////////
// Legendre poly.
////////////////////////////////////////////////////////////////


double mathAssocLegendreP(int l, int m, double x) {
	double p = boost::math::legendre_p(l,m,x);
	return p;
}

double mathCondonShortleyPhase(int m) {
	int am = abs(m);
	if (am%2 == 0) {
		return 1;
	} else {
		return -1;
	}
}


////////////////////////////////////////////////////////////////
// Spherical Harmonics
////////////////////////////////////////////////////////////////

complex_t mathSphHarmonic(int n, int m, double theta, double phi) {
	complex_t ynm;
	if (m >= 0) {
		ynm = boost::math::spherical_harmonic(n, m, theta, phi);
	} else {
		ynm = boost::math::spherical_harmonic(n, -m, theta, phi);
		ynm = std::conj(ynm);
	}

	if (1) { // further include phase factor
		ynm *= mathCondonShortleyPhase(m);
	}

	if (1) { 
		// cancel full normalization factor
		// so it falls back to the pseodo-regularized Schmidt SH
		double coef = sqrt(((double) 2*n+1)/(M_PI*4));
		ynm /= coef;
	}

	return ynm;
}


////////////////////////////////////////////////////////////////
// Gamma
////////////////////////////////////////////////////////////////
double mathGamma(double x) {
	double val = boost::math::tgamma(x);
	return val;
}

double mathFactorial(int n) {
	double val = boost::math::factorial<double>(n);
	return val;
}


double mathDoubleFactorial(int n) {
	double val = 1.0;
	for (int i=1; i<=n; i+=2) {
		val *= (double) i;
	}
	return val;
}

