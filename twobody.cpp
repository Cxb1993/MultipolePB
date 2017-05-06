
#include <utility>

#include "mathfunc.h"
#include "twobody.h"

static double func_a(int nu, int n, int m) {
	double a1 = mathGamma(n-nu+0.5) * mathGamma(m-nu+0.5) * mathGamma(nu+0.5) / mathGamma(m+n-nu+1.5);
	double a2 = mathFactorial(n+m-nu) / mathFactorial(n-nu) / mathFactorial(m-nu) / mathFactorial(nu);
	double a = ((double) n + m - 2*nu + 0.5) / M_PI * a1 * a2;
	return a;
}

static double func_b(int n, int m, double kr) {
	const int numax = std::min(n,m);

	double bnm = 0;
	for (int nu=0; nu<=numax; nu++) {
		double anm = func_a(nu, n, m);
		double knm = mathModSphBesselK(n+m-2*nu,kr);
		
		bnm += anm * knm;
	}

	return bnm;
}


void solveTwoBodyConstPotent(
	TwoBodyProblem &prob,
	double a1, double a2, double h, 
	double phi1, double phi2, 
	double kappa, int ncutoff)
{
	//
	prob.a1 = a1;
	prob.a2 = a2;
	prob.h = h;
	prob.psi1 = phi1;
	prob.psi2 = phi2;
	prob.kappa = kappa;
	prob.ncutoff = ncutoff;

	prob.acoef.resize(ncutoff);
	prob.bcoef.resize(ncutoff);
	prob.kacoef.resize(ncutoff);
	prob.kbcoef.resize(ncutoff);


	//
	const double a1k = a1 * kappa;
	const double a2k = a2 * kappa;

	const double r = a1 + a2 + h;
	const double rk = r * kappa;

	// precalc Bessel
	VectorXd Ka1k(ncutoff);
	VectorXd Ia1k(ncutoff);
	VectorXd Ka2k(ncutoff);
	VectorXd Ia2k(ncutoff);
	for (int n=0; n<ncutoff; n++) {
		Ka1k(n) = mathModSphBesselK(n, a1k);
		Ia1k(n) = mathModSphBesselI(n, a1k);
		Ka2k(n) = mathModSphBesselK(n, a2k);
		Ia2k(n) = mathModSphBesselI(n, a2k);
	}

	// B
	MatrixXd Bmat(ncutoff,ncutoff);
	for (int n=0; n<ncutoff; n++) {
		for (int m=0; m<ncutoff; m++) {
			double bnm = func_b(n, m, rk);
			Bmat(n,m) = bnm;
		}
	}

	// L,M
	MatrixXd Lmat(ncutoff,ncutoff);
	MatrixXd Mmat(ncutoff,ncutoff);
	for (int j=0; j<ncutoff; j++) {
		for (int n=0; n<ncutoff; n++) {
			double Bnj = Bmat(n,j);
			
			double Ljn = ((double) 2*j+1) * Bnj * Ia1k(j) / Ka2k(n);
			double Mjn = ((double) 2*j+1) * Bnj * Ia2k(j) / Ka1k(n);

			Lmat(j,n) = Ljn;
			Mmat(j,n) = Mjn;
		}
	}

	// identity
	MatrixXd Imat = MatrixXd::Identity(ncutoff,ncutoff);

	//
	MatrixXd mat(ncutoff*2,ncutoff*2);
	mat << Imat, Lmat, Mmat, Imat;

	//
	VectorXd rhs1(ncutoff);
	VectorXd rhs2(ncutoff);
	rhs1.setZero();
	rhs1(0) = phi1;
	rhs2.setZero();
	rhs2(0) = phi2;

	VectorXd rhs(ncutoff*2);
	rhs << rhs1, rhs2;


	//
	VectorXd sol(ncutoff*2);
	sol = mat.colPivHouseholderQr().solve(rhs);

	prob.kacoef = sol.block(0,0,ncutoff,1);
	prob.kbcoef = sol.block(ncutoff,0,ncutoff,1);

	prob.acoef = prob.kacoef.cwiseQuotient(Ka1k);
	prob.bcoef = prob.kbcoef.cwiseQuotient(Ka2k);

	//for (int i=0; i<ncutoff; i++) {
	//	prob.acoef(i) /= Ka1k(i);
	//}
	//for (int i=0; i<ncutoff; i++) {
	//	prob.bcoef(i) /= Ka2k(i);
	//}

}

