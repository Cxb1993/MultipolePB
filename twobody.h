#pragma once

//#include "mppb.h"
#include "mathvecmat.h"



struct TwoBodyProblem
{
	double a1, a2, h;
	double psi1, psi2;
	double kappa;

	//
	int ncutoff;

	// 
	VectorXd acoef, bcoef;
	// coef multiplied by Bessel func
	VectorXd kacoef, kbcoef;

};

void solveTwoBodyConstPotent(
	TwoBodyProblem &prob,
	double a1, double a2, double h, 
	double phi1, double phi2, 
	double kappa, int ncutoff);



