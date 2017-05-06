
#include <cassert>

#include <iostream>
#include <utility>
#include <algorithm>
#include <vector>

#include "mathfunc.h"
#include "mathvecmat.h"

#include "multibody.h"

//
namespace 
{
	const complex_t ImUnit(0.0, 1.0);
}


inline int sh_numpole(int nmax) {
	int nmax1 = nmax + 1;
	int npole = nmax1 * nmax1;
	return npole;
}

inline int sh_sub2ind(int n, int m) {
	int ind = n + m + n*n;
	return ind;
}

inline void sh_cart2sph(
	double x, double y, double z,
	double &r, double &theta, double &phi)
{
	r = sqrt(x*x + y*y + z*z);
	theta = acos(z/r);
	phi = atan2(y,x);
	if (phi < 0) {
		phi = phi + M_PI*2;
	}
}

inline Vector3d sh_cart2sph(const Vector3d &vin) {
	Vector3d vout;
	sh_cart2sph(vin(0),vin(1),vin(2), vout(0),vout(1),vout(2));
	return vout;
}

struct ReExpansion
{
	//
	double kappa;

	// coord transform
	double dist;
	double theta;
	double phi;
	double chi;
	Matrix3d Q;

	//
	int nmax;
	int npole;

	MatrixXcd R;
	MatrixXd S;
};

double reexpRegularFunc(int n, double r, double kappa) {
	double ihat = mathAdaptModSphBesselI(n,kappa*r);
	return pow(r,n) * ihat;
}

double reexpSingularFunc(int n, double r, double kappa) {
	double khat = mathAdaptModSphBesselK(n,kappa*r);
	return exp(-kappa*r) / pow(r,n+1) * khat;
}

// F  = K(r,kappa)
// dF = dK/dr
void reexpRecurSingular(int nmax, double r, double kappa,
	double f[], double dfdr[]) 
{
	assert(nmax >= 0);
	assert(r > 0);
	assert(kappa >= 0);

	// 
	const double x = kappa * r;
	const double r2 = r * r;
	const double kappa2 = kappa * kappa;
	const double enx = exp(-x);


	// recurrent variables
	double fm = 0;
	double fn = 0;
	double fp = 0;

	{
		// initial value for n=0,1
		fm = 0;
		fn = enx / r;
		fp = enx / r2 * (x+1);

		if (f) {
			f[0] = fn;
		}

		if (dfdr) {
			// derivative n=0 is -K1
			dfdr[0] = -fp;
		}
	}
	
	for (int n=1; n<=nmax; n++) {

		fm = fn;
		fn = fp;

		// next function
		fp = kappa2 * fm + (double) (2*n+1) / r * fn;

		if (f != NULL) {
			f[n] = fn;
		}

		if (dfdr != NULL) {
			double dfn = -(double) n * kappa2 * fm - (double) (n+1) * fp;
			dfn /= (2*n+1);
			dfdr[n] = dfn;
		}
	}

	if (1) {
		// add (2n-1)!! denominator
		for (int n=1; n<=nmax; n++) {
			
			double coef = 1.0 / mathDoubleFactorial(2*n-1);
			
			if (f) {
				f[n] *= coef;
			}

			if (dfdr) {
				dfdr[n] *= coef;
			}
		}
	}
}

// F  = I(r,kappa)
// dF = dI/dr
void reexpRecurRegular(int nmax, double r, double kappa,
	double f[], double dfdr[]) 
{
	assert(nmax >= 0);
	assert(r >= 0);
	assert(kappa > 0);

	// 
	const double x = kappa * r;
	const double sinhx = sinh(x);
	const double coshx = cosh(x);

	const double xsmall = 1.0e-8;

	// recurrent variables
	double fm = 0;
	double fn = 0;
	double fp = 0;

	{
		// initial value for n=0,1
		fm = 0;
		if (x > xsmall) {
			fn = sinhx / x;
			fp = (x*coshx - sinhx) / (x*x);
		} else { // lim x->0
			fn = 1;
			fp = 0;
		}

		if (f) {
			f[0] = fn;
		}

		if (dfdr) {
			// derivative n=0 is I1
			dfdr[0] = fp;
		}
	}
	
	for (int n=1; n<=nmax; n++) {

		fm = fn;
		fn = fp;

		// next function
		if (x > xsmall) {
			fp = fm - (double) (2*n+1) / x * fn;
		} else {
			fp = 0;
		}

		if (f != NULL) {
			f[n] = fn;
		}

		if (dfdr != NULL) {
			double dfn = (double) n * fm + (double) (n+1) * fp;
			dfn /= (2*n+1);
			dfdr[n] = dfn;
		}
	}

	if (1) {
		// 
		for (int n=0; n<=nmax; n++) {
			
			double coef = mathDoubleFactorial(2*n+1) / pow(kappa,n);
			
			if (f) {
				f[n] *= coef;
			}

			if (dfdr) {
				dfdr[n] *= coef * kappa;
			}
		}
	}
}


void reexpSetNcutoff(ReExpansion &reexp, const int nmax) {
	// multipole number for given n_max
	const int npole = sh_numpole(nmax);

	reexp.nmax = nmax;
	reexp.npole = npole;

	reexp.R.resize(npole,npole);
	reexp.S.resize(npole,npole);
}


void reexpCoordTrans(const Vector3d &pab,
	ReExpansion &rd)
{
	// original cartesian coord
	Vector3d ex(1.0, 0.0, 0.0);
	Vector3d ey(0.0, 1.0, 0.0);
	Vector3d ez(0.0, 0.0, 1.0);
	
	Matrix3d E; E << ex,ey,ez;

	// rotated coaxial coord
	// use center connection as new z-axis
	Vector3d ezhat = pab;
	ezhat.normalize();

	// new x-axis
	Vector3d exhat = ezhat.cross(ez);
	if (exhat.norm() < 1.0e-10) {
		exhat = ex;
	}
	exhat.normalize();

	// new y-axis
	Vector3d eyhat = ezhat.cross(exhat);
	eyhat.normalize();

	Matrix3d Ehat; Ehat << exhat,eyhat,ezhat;

	// transform matrix for xhat = Q.x
	Matrix3d Qmat = Ehat.transpose() * E;

	// original z-axis in new coord
	Vector3d ezprime = Qmat * ez;
	double rprime, thetaprime, phiprime;
	sh_cart2sph(ezprime(0),ezprime(1),ezprime(2), rprime,thetaprime,phiprime);

	// new z-axis in original coord
	double anggamma, angchi;
	sh_cart2sph(ezhat(0),ezhat(1),ezhat(2), rprime, anggamma, angchi);

	//
	rd.dist = pab.norm();
	rd.theta = thetaprime;
	rd.phi = phiprime;
	rd.chi = angchi;
	rd.Q = Qmat;

}

static inline double reexp_anm(double n, double m) {
	return sqrt((n+m+1) * (n-m+1) / (n*2+1) / (n*2+3));
}
static inline double reexp_bnm(double n, double m) {
	double s = (m>=0 ? 1 : -1);
	return s * sqrt((n-m-1) * (n-m) / (n*2-1) / (n*2+1));
}

//static inline complex_t& REntry(MatrixXcd &R, int n, int m, int s) {
//	return R(sh_sub2ind(n,m),sh_sub2ind(n,s));
//}

void reexpRotate(ReExpansion &reexp) 
{
	// cutoff
	const int nmax = reexp.nmax;
	const int npole = reexp.npole;
	const int nmax1 = nmax + 1;
	
	// rotation
	const double theta = reexp.theta;
	const double phi = reexp.phi;
	const double chi = reexp.chi;

	complex_t epphi = std::exp(ImUnit * phi);
	complex_t emphi = std::exp(-ImUnit * phi);
	double cost = cos(theta);
	double sint = sin(theta);
	complex_t echi = std::exp(ImUnit * chi);

	// buffer for recurrent relation
	MatrixXcd Rmat(npole*4,npole*4);
	Rmat.setZero();

#define R(n,m,s) Rmat(sh_sub2ind((n),(m)),sh_sub2ind((n),(s)))

	// step 1, initial value with m=0
	for (int n=0; n<=2*nmax1-1; n++) {
		for (int s=-n; s<=n; s++) {
			const int m = 0;
			R(n,m,s) = mathSphHarmonic(n,-s,theta,phi);
		}
	}

	// step 2
	for (int m=0; m<=nmax; m++) {
		for (int n=m+2; n<=2*nmax1-m-1; n++) {
			double bnm = reexp_bnm(n,m);

			for (int s=-n+1; s<=n-1; s++) {
				double b1 = reexp_bnm(n,s-1);
				complex_t R1 = R(n,m,s-1);
				R1 *= 0.5 * emphi * (1+cost) * b1;

				double b2 = reexp_bnm(n,-s-1);
				complex_t R2 = R(n,m,s+1);
				R2 *= 0.5 * epphi * (1-cost) * b2;

				double a3 = reexp_anm(n-1,s);
				complex_t R3 = R(n,m,s);
				R3 *= sint * a3;

				complex_t Rnew = -echi / bnm * (R1-R2+R3);
				R(n-1,m+1,s) = Rnew;
			}
		}
	}

	// step 3, fill rest
	for (int n=0; n<=nmax; n++) {
		for (int m=-n; m<=-1; m++) {
			for (int s=-n; s<=n; s++) {
				R(n,m,s) = std::conj(R(n,-m,-s));
			}
		}
	}

#undef R

	reexp.R = Rmat.block(0,0,npole,npole);

}

static inline double reexp_alpha(double n, double m) {
	return sqrt((n+m+1) * (n-m+1));
}
static inline double reexp_beta(double n, double m, double kappa) {
	return kappa*kappa * reexp_alpha(n,m) / (n*2+1) / (n*2+3);
}
static inline double reexp_eta(double n, double m) {
	double s = (m>=0 ? 1 : -1);
	return s * sqrt((n-m-1)*(n-m));
}
static inline double reexp_mu(double n, double m, double kappa) {
	return kappa*kappa * reexp_eta(n,m) / (n*2-1) / (n*2+1);
}

void reexpScale(ReExpansion &reexp)
{
	//
	const double kappa = reexp.kappa;
	const double r = reexp.dist;
	const double kr = kappa * r;
	
	// cutoff
	const int nmax = reexp.nmax;
	const int npole = reexp.npole;
	const int nmax1 = nmax + 1;

	// buffer
	MatrixXd Smat(npole*4,npole*4);
	Smat.setZero();

#define S(n,l,m) Smat(sh_sub2ind(n,m),sh_sub2ind(l,m))

	// step 1a, n=m=0
	for (int l=0; l<2*nmax1; l++) {
		const int n = 0;
		const int m = 0;
		S(n,l,m) = 1.0 / pow(r,l) * mathAdaptModSphBesselK(l,kr) * exp(-kr)/r;
	}

	// step 1b, l=m=0
	for (int n=0; n<2*nmax1; n++) {
		const int l = 0;
		const int m = 0;
		S(n,l,m) = pow(-1.0,n) * S(l,n,m);
	}

	// step 2, m=0
	for (int n=0; n<=nmax-1; n++) {
		for (int l=0; l<=2*nmax-n; l++) {
			const int m = 0;

			double beta1 = reexp_beta(l-1,m,kappa);
			double S1 = S(n,l-1,m);

			double beta2 = 0;
			double S2 = 0;
			if (n > 0) {
				beta2 = reexp_beta(n-1,m,kappa);
				S2 = S(n-1,l,0);
			}

			double alpha3 = reexp_alpha(l,m);
			double S3 = S(n,l+1,m);

			double alphanew = reexp_alpha(n,m);
			double Snew = -1.0/alphanew * (beta1*S1 + beta2*S2 + alpha3*S3);

			S(n+1,l,m) = Snew;
		}
	}

	// step 3
	for (int m=1; m<=nmax; m++) {
		
		// n=m
		for (int l=m; l<=2*nmax-m; l++) {
			double mu1 = reexp_mu(l,-m, kappa);
			double S1 = S(m-1,l-1,m-1);

			double eta2 = reexp_eta(l+1,m-1);
			double S2 = S(m-1,l+1,m-1);

			double etanew = reexp_eta(m,-m);
			double Snew = -1.0/etanew * (mu1*S1 + eta2*S2);

			S(m,l,m) = Snew;
		}

		//
		for (int n=m; n<=nmax-1; n++) {
			for (int l=n+1; l<=2*nmax-n-1; l++) {
				double beta1 = reexp_beta(l-1,m, kappa);
				double S1 = S(n,l-1,m);

				double beta2 = 0;
				double S2 = 0;
				if (n > m) {
					beta2 = reexp_beta(n-1,m, kappa);
					S2 = S(n-1,l,m);
				}

				double alpha3 = reexp_alpha(l,m);
				double S3 = S(n,l+1,m);

				double alphanew = reexp_alpha(n,m);
				double Snew = -1.0/alphanew * (beta1*S1 + beta2*S2 + alpha3*S3);

				S(n+1,l,m) = Snew;
			}
		}
	}

	// step 4, fill rest part
	for (int n=1; n<=nmax; n++) {
		for (int m=-n; m<=n; m++) {
			for (int l=abs(m); l<=nmax; l++) {
				if (n > l) {
					S(n,l,m) = pow(-1.0,n+l) * S(l,n,m);
				}
			}
		}
	}

	for (int n=1; n<=nmax; n++) {
		for (int m=-n; m<=-1; m++) {
			for (int l=abs(m); l<=nmax; l++) {
				S(n,l,m) = S(n,l,-m);
			}
		}
	}

#undef S

	reexp.S = Smat.block(0,0,npole,npole);
}


void reexpMatrix(MatrixXcd &T,
	const Vector3d &psrc, const Vector3d &pdst, double kappa, int nmax)
{
	// src -> dst
	const Vector3d pab = pdst - psrc;

	// descriptor for expansion
	ReExpansion rd;
	rd.kappa = kappa;

	//
	reexpSetNcutoff(rd, nmax);

	// 
	reexpCoordTrans(pab, rd);

	//std::cout << rd.chi << std::endl;
	//std::cout << rd.dist << std::endl;
	//std::cout << rd.theta << std::endl;
	//std::cout << rd.phi << std::endl;
	//std::cout << rd.Q << std::endl;


	//
	reexpRotate(rd);
	//std::cout << rd.R << std::endl;

	// 
	reexpScale(rd);
	//std::cout << rd.S << std::endl;

	T.resize(rd.npole,rd.npole);
	T = rd.R * rd.S * rd.R.adjoint();
}


void solveMultiBody(MultiSphereProblem &prob) 
{
	// parameters
	const double kappa = prob.invDebyeLength;

	const int nmax = prob.ncutoff;
	const int npole = sh_numpole(nmax);
	std::cout << "nmax=" << nmax << std::endl;
	std::cout << "npole=" << npole << std::endl;

	// spheres 
	const int nsph = prob.numSphere;
	std::cout << "nsph=" << nsph << std::endl;
	const std::vector<Vector3d> &sphCenter = prob.position;
	const std::vector<double> &sphRadius = prob.radius;
	const std::vector<double> &sphPotent = prob.zetaPotent;

	// build linear system
	const int ndof = npole * nsph;
	MatrixXcd A(ndof,ndof);
	VectorXcd rhs(ndof);
	A.setZero();
	rhs.setZero();

	for (int isph=0; isph<nsph; isph++) {
		// the current body
		const Vector3d ipos = sphCenter[isph];
		const double irad = sphRadius[isph];

		// offset in matrix
		const int ioff = isph * npole;

		// radial function
		VectorXd khat(npole);
		VectorXd ihat(npole);
		for (int n=0; n<=nmax; n++) {
			double kn = reexpSingularFunc(n,irad,kappa);
			double in = reexpRegularFunc(n,irad,kappa);
			for (int m=-n; m<=n; m++) {
				int inm = sh_sub2ind(n,m);
				khat(inm) = kn;
				ihat(inm) = in;
			}
		}

		// build matrix
		for (int jsph=0; jsph<nsph; jsph++) {
			if (isph == jsph) { // self
				// 
				A.block(ioff,ioff,npole,npole) = khat.asDiagonal();
			} else { // (i<-j) reexpansion pair
				// 
				const Vector3d jpos = sphCenter[jsph];
				// offset in matrix
				const int joff = jsph * npole;

				// reexpansion matrix
				MatrixXcd Tij(npole,npole);
				reexpMatrix(Tij, jpos, ipos, kappa, nmax);

				//
				A.block(ioff,joff,npole,npole) = ihat.asDiagonal() * Tij.conjugate();
			}
		}

		// build rhs according to BC
		rhs(ioff) = sphPotent[isph];
	}

	// solve 
	VectorXcd sol = A.colPivHouseholderQr().solve(rhs);

	// output coefficients
	prob.coef.resizeLike(sol);
	//prob.coef.resize(ndof,1);
	prob.coef = sol;
	prob.coef.resize(npole,nsph);
}

void test_multisphere() {

	//const double kappa = 1.0;
	const double kappa = 0.5;

	const int nmax = 6;
	const int npole = sh_numpole(nmax);
	std::cout << "nmax=" << nmax << std::endl;
	std::cout << "npole=" << npole << std::endl;

	// setup problem
	MultiSphereProblem prob;

	prob.invDebyeLength = kappa;

	prob.position.push_back(Vector3d(0,0,0));
	prob.radius.push_back(1.0);
	prob.zetaPotent.push_back(1.0);

	prob.position.push_back(Vector3d(3,1,0));
	prob.radius.push_back(1.0);
	prob.zetaPotent.push_back(1.0);

	prob.position.push_back(Vector3d(1,3,0));
	prob.radius.push_back(1.0);
	prob.zetaPotent.push_back(1.0);

	const int nsph = prob.position.size();
	prob.numSphere = nsph;
	std::cout << "nsph=" << nsph << std::endl;

	prob.ncutoff = nmax;

	//
	solveMultiBody(prob);

	////
	//const int ndof = npole * nsph;
	//MatrixXcd A(ndof,ndof);
	//VectorXcd rhs(ndof);
	//A.setZero();
	//rhs.setZero();

	//for (int isph=0; isph<nsph; isph++) {
	//	// 
	//	const Vector3d ipos = sphCenter[isph];
	//	const double irad = sphRadius[isph];

	//	//
	//	const int ioff = isph * npole;

	//	//
	//	VectorXd khat(npole);
	//	VectorXd ihat(npole);
	//	for (int n=0; n<=nmax; n++) {
	//		double kn = reexpSingularFunc(n,irad,kappa);
	//		double in = reexpRegularFunc(n,irad,kappa);
	//		for (int m=-n; m<=n; m++) {
	//			int inm = sh_sub2ind(n,m);
	//			khat(inm) = kn;
	//			ihat(inm) = in;
	//		}
	//	}

	//	// build matrix
	//	for (int jsph=0; jsph<nsph; jsph++) {
	//		if (isph == jsph) { // self
	//			// 
	//			A.block(ioff,ioff,npole,npole) = khat.asDiagonal();
	//		} else { // (i<-j) reexpansion pair
	//			// 
	//			const Vector3d jpos = sphCenter[jsph];
	//			// offset in matrix
	//			const int joff = jsph * npole;

	//			// reexpansion matrix
	//			MatrixXcd Tij(npole,npole);
	//			reexpMatrix(Tij, jpos, ipos, kappa, nmax);

	//			//
	//			A.block(ioff,joff,npole,npole) = ihat.asDiagonal() * Tij.conjugate();
	//		}
	//	}

	//	// build rhs according to BC
	//	rhs(ioff) = 1.0;
	//}

	//// solve 
	//VectorXcd sol = A.colPivHouseholderQr().solve(rhs);

	//// coefficients
	//MatrixXcd coef = sol;
	//coef.resize(npole,nsph);



	// plot
	if (1) {
		const MatrixXcd &coef = prob.coef;
		const std::vector<Vector3d> &sphCenter = prob.position;
		const std::vector<double> &sphRadius = prob.radius;


		double xlo = -2;
		double xhi = 6;
		double ylo = -2;
		double yhi = 4;
		//int nx = 81;
		//int ny = 61;
		const int nx = 128 + 1;
		const int ny = 96 + 1;
		double dx = (xhi-xlo)/(nx-1);
		double dy = (yhi-ylo)/(ny-1);

		FILE *fp = fopen("tmp0.csv","w");
		if (!fp) {
			fprintf(stderr, "failed to save\n");
			exit(1);
		}
		fprintf(fp, "x,y,z,ok,real,imag\n");

		for (int j=0; j<ny; j++) {
		for (int i=0; i<nx; i++) {
			Vector3d ptest(xlo+dx*i,ylo+dy*j,0.0);

			int ok = 1;
			std::vector<Vector3d> pcoord;
			for (int isph=0; isph<nsph; isph++) {
				Vector3d psph = sh_cart2sph(ptest-sphCenter[isph]);
				pcoord.push_back(psph);

				if (psph(0) < sphRadius[isph]) {
					ok = 0;
				}
			}

			complex_t val = 0;

			if (ok) {
				for (int isph=0; isph<nsph; isph++) {

					double r = pcoord[isph](0);
					double theta = pcoord[isph](1);
					double phi = pcoord[isph](2);

					complex_t ival = 0;
					for (int n=0; n<=nmax; n++) {
						for (int m=-n; m<=n; m++) {
							int inm = sh_sub2ind(n,m);
							ival += coef(inm,isph) * reexpSingularFunc(n,r,kappa) * mathSphHarmonic(n,m,theta,phi);
						}
					}

					val += ival;
				}
			}

			fprintf(fp,"%lf,%lf,%lf,%d,%lf,%lf\n",ptest(0),ptest(1),ptest(2), ok,val.real(),val.imag());
		}
		}

		fclose(fp);

	}

}



void hoge() {
	if (0) {
		const int nmax = 8;
		VectorXd ks(nmax);
		VectorXd is(nmax);
		for (int n=0; n<nmax; n++) {
			ks(n) = mathAdaptModSphBesselK(n,1.0);
			is(n) = mathAdaptModSphBesselI(n,1.0);
		}
		std::cout << ks << std::endl;
		std::cout << is << std::endl;
	}

	if (0) {
		const int nmax = 3;
		const int npole = sh_numpole(nmax);
		VectorXcd ys(npole);
		for (int n=0; n<=nmax; n++) {
			for (int m=-n; m<=n; m++) {
				int nm = sh_sub2ind(n,m);
				ys(nm) = mathSphHarmonic(n,m,M_PI/1.5,M_PI/2.5);
			}
		}

		std::cout << ys << std::endl;
	}

	if (0) {
		const double kappa = 1.0;
		const int nmax = 3;
		const int npole = sh_numpole(nmax);
		std::cout << "nmax=" << nmax << std::endl;
		std::cout << "npole=" << npole << std::endl;

		Vector3d pa(0,0,0);
		Vector3d pb(-1,-1,4);

		//reexpMatrix(pb,pa, kappa, nmax);
	}

	if (0) {
		const double kappa = 0.1;
		const int nmax = 6;
		const int npole = sh_numpole(nmax);
		std::cout << "nmax=" << nmax << std::endl;
		std::cout << "npole=" << npole << std::endl;

		const Vector3d pa(0,0,0);
		const Vector3d pb(3,1,0);
		const double ra = 1.0;
		const double rb = 1.0;

		MatrixXcd Tab(npole,npole);
		reexpMatrix(Tab, pb, pa, kappa, nmax);

		MatrixXcd Tba(npole,npole);
		reexpMatrix(Tba, pa, pb, kappa, nmax);

		VectorXd khata(npole);
		VectorXd ihata(npole);
		VectorXd khatb(npole);
		VectorXd ihatb(npole);
		for (int n=0; n<=nmax; n++) {
			double kna = reexpSingularFunc(n,ra,kappa);
			double ina = reexpRegularFunc(n,ra,kappa);
			double knb = reexpSingularFunc(n,rb,kappa);
			double inb = reexpRegularFunc(n,rb,kappa);

			for (int m=-n; m<=n; m++) {
				int inm = sh_sub2ind(n,m);
				khata(inm) = kna;
				ihata(inm) = ina;
				khatb(inm) = knb;
				ihatb(inm) = inb;
			}
		}

		//
		MatrixXcd La = khata.asDiagonal();
		MatrixXcd Ma = ihata.asDiagonal() * Tab.conjugate();
		MatrixXcd Lb = khatb.asDiagonal();
		MatrixXcd Mb = ihatb.asDiagonal() * Tba.conjugate();

		VectorXcd rhsa(npole);
		rhsa.setZero();
		rhsa(0) = 1.0;
		VectorXcd rhsb(npole);
		rhsb.setZero();
		rhsb(0) = 1.0;
		
		MatrixXcd A(npole*2,npole*2);
		A << La,Ma,Mb,Lb;
		
		VectorXcd rhs(npole*2);
		rhs << rhsa,rhsb;

		VectorXcd sol = A.colPivHouseholderQr().solve(rhs);

		VectorXcd coefa = sol.block(0,0,npole,1);
		VectorXcd coefb = sol.block(npole,0,npole,1);

		//std::cout << coefa << std::endl;

		{
			double xlo = -2;
			double xhi = 6;
			double ylo = -2;
			double yhi = 4;
			//int nx = 81;
			//int ny = 61;
			const int nx = 128 + 1;
			const int ny = 96 + 1;
			double dx = (xhi-xlo)/(nx-1);
			double dy = (yhi-ylo)/(ny-1);

			FILE *fp = fopen("tmp0.csv","w");
			if (!fp) {
				fprintf(stderr, "failed to save\n");
				exit(1);
			}
			fprintf(fp, "x,y,z,ok,real,imag\n");

			for (int j=0; j<ny; j++) {
			for (int i=0; i<nx; i++) {
				Vector3d ptest(xlo+dx*i,ylo+dy*j,0.0);

				double rada,thetaa,phia;
				sh_cart2sph(ptest(0)-pa(0),ptest(1)-pa(1),ptest(2)-pa(2),rada,thetaa,phia);

				double radb,thetab,phib;
				sh_cart2sph(ptest(0)-pb(0),ptest(1)-pb(1),ptest(2)-pb(2),radb,thetab,phib);

				int ok = rada>=ra && radb>=rb;
				complex_t val = 0;

				if (ok) {
					complex_t vala = 0;
					complex_t valb = 0;
					for (int n=0; n<=nmax; n++) {
						for (int m=-n; m<=n; m++) {
							int inm = sh_sub2ind(n,m);
							vala += coefa(inm) * reexpSingularFunc(n,rada,kappa) * mathSphHarmonic(n,m,thetaa,phia);
							valb += coefb(inm) * reexpSingularFunc(n,radb,kappa) * mathSphHarmonic(n,m,thetab,phib);
						}
					}
					val = vala + valb;
				}

				fprintf(fp,"%lf,%lf,%lf,%d,%lf,%lf\n",ptest(0),ptest(1),ptest(2), ok,val.real(),val.imag());
			}
			}

			fclose(fp);

		}


	}

	if (0) {
		test_multisphere();
	}

	if (1) {
		const int nmax = 4;
		VectorXd ks(nmax+1), kprime(nmax+1);
		VectorXd is(nmax+1), iprime(nmax+1);

		const double r = 1.5;
		const double kappa = 2.0;

		reexpRecurSingular(nmax, r, kappa, ks.data(), kprime.data());

		reexpRecurRegular(nmax, r, kappa, is.data(), iprime.data());

		std::cout << ks << std::endl;
		std::cout << kprime << std::endl;
		std::cout << is << std::endl;
		std::cout << iprime << std::endl;
	}
}




