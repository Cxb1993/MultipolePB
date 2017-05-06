
#include <cstdio>
#include <cstdlib>
#include <cmath>


#include <iostream>

#include "mppb.h"
#include "mathfunc.h"

//
#include "twobody.h"

//


static void test_specfunc() {
	double p = mathAssocLegendreP(1, 1, 0.5);
	std::cout << p << std::endl;

	double i = mathModSphBesselI(2, 0.5);
	std::cout << i << std::endl;

	double k = mathModSphBesselK(2, 0.5);
	std::cout << k << std::endl;

	double g = mathGamma(4.0);
	std::cout << g << std::endl;

	complex_t ya = mathSphHarmonic(2, 1, M_PI/3, M_PI/5);
	complex_t yb = mathSphHarmonic(2, -1, M_PI/3, M_PI/5);
	std::cout << ya << std::endl;
	std::cout << yb << std::endl;
}

static void test_twobody(int argc, char *argv[]) {

	if (argc != 1 + 3) {
		std::cerr << "twobody h kappa cutoff" << std::endl;
		exit(1);
	}

	const double a1 = 1.0;
	const double a2 = 1.0;
	const double phi1 = 1.0;
	const double phi2 = 1.0;

	const double h = atof(argv[1]);
	const double kappa = atof(argv[2]);
	const int cutoff = atoi(argv[3]);

	TwoBodyProblem prob;

	solveTwoBodyConstPotent(prob, 
		a1, a2, h, phi1, phi2, kappa, cutoff);

	std::cout << "a = " << std::endl;
	std::cout << prob.acoef << std::endl;
	
	std::cout << "aK = " << std::endl;
	std::cout << prob.kacoef << std::endl;
}



int main(int argc, char *argv[]) {

	std::cout << "hoge" << std::endl;

	//test_specfunc();

	//test_twobody(argc, argv);

	extern void hoge();
	hoge();

	return 0;
}


