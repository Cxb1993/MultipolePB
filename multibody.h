#pragma once

#include <vector>
#include "mathvecmat.h"



struct MultiSphereProblem
{
	// 
	// physical parameter
	//

	// kappa
	double invDebyeLength;

	// number
	int numSphere;

	// sphere center
	std::vector<Vector3d> position;

	// sphere size
	std::vector<double> radius;

	// zeta-potential as BC
	std::vector<double> zetaPotent;



	//
	// discretized
	//

	int ncutoff;

	MatrixXcd coef;

};


