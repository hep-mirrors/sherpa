#include "METOOLS/HadronCurrents/FormFactors/RChL_Functions.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include <cmath>

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

double RChL::H(const double & x,const double & y,
	      const double & lambda_p,const double & lambda_pp) {
  double lambda0 = (lambda_p+lambda_pp)/4.;
  return -lambda0*y + lambda_p*x + lambda_pp;
}

double RChL::Lambda_p(const double & F,const double & FA,const double & GV) {
  return sqr(F)/(2.*sqrt(2.)*FA*GV);
}

double RChL::Lambda_pp(const double & F,const double & GV,
		      const double & lambda_p) {
  // Eq.(8): lambda'' = -(1 - 2 GV^2/F^2) lambda'.  F, not FV - see the header.
  return -(1.-2.*sqr(GV)/sqr(F))*lambda_p;
}

Complex RChL::alpha2(const double & q2,const double & s1,const double & s2,
		    const double & s3,const double & mpi2,
		    const double & Mrho,const double & GammaRho_s1,
		    const double & GV,const double & FV) {
  Complex denom(sqr(Mrho)-s1, -Mrho*GammaRho_s1);
  return (3.*GV/FV) * (s1/q2) * (mpi2/(q2-mpi2)) * (s3-s2)/denom;
}

double RChL::AR(const double & q2,const double & x,const double & y,
	       const double & m1_2,const double & m2_2,const double & m3_2,
	       const double & GV,const double & FV) {
  return ( 3.*x + m1_2 - m3_2 +
	   (1.-2.*GV/FV) * (2.*q2 - 2.*x - y + m3_2 - m2_2) );
}

double RChL::BR(const double & x,const double & y,
	       const double & m1_2,const double & m2_2,
	       const double & GV,const double & FV) {
  return ( 2.*(m2_2-m1_2) + (1.-2.*GV/FV) * (y - x + m1_2 - m2_2) );
}

double RChL::ARR(const double & q2,const double & x,const double & y,
		const double & m1_2,const double & m2_2,const double & m3_2,
		const double & lambda_p,const double & lambda_pp) {
  return ( (lambda_p+lambda_pp)*(-3.*x+m3_2-m1_2) +
	   (2.*q2+x-y+m1_2-m2_2) * H(x/q2,m2_2/q2,lambda_p,lambda_pp) );
}

double RChL::BRR(const double & q2,const double & x,const double & y,
		const double & z,const double & m1_2,const double & m2_2,
		const double & m3_2,
		const double & lambda_p,const double & lambda_pp) {
  return ( 2.*(lambda_p+lambda_pp)*(m1_2-m2_2) +
	   (y-x+m2_2-m1_2) * H(z/q2,m3_2/q2,lambda_p,lambda_pp) );
}

double RChL::CR(const double & q2,const double & x,
	       const double & m1_2,const double & m2_2,const double & m3_2,
	       const double & c125,const double & c1256,const double & c1238,
	       const double & c4) {
  return ( c125*q2 - c1256*x + c1238*m3_2 + 8.*c4*(m1_2-m2_2) );
}

double RChL::CRR(const double & q2,const double & x,const double & m2,
		const double & d3,const double & d123) {
  return ( d3*(q2+x) + d123*m2 );
}

double RChL::DR(const double & q2,const double & x,const double & y,
	       const double & mK2,const double & mpi2,
	       const double & g13,const double & g2,
	       const double & g4,const double & g5) {
  double g123 = g13 + 2.*g2; // g1+2g2-g3
  return ( g123*(x+y) - 2.*g2*(q2+mK2) -
	   g13*(3.*mK2+mpi2) + 2.*g4*(mK2+mpi2) + 2.*g5*mK2 );
}

double RChL::ER(const double & x,const double & y,
	       const double & g13,const double & g2) {
  double g123 = g13 + 2.*g2;
  return g123*(x-y);
}

Complex RChL::BWsigma(const double & x,const double & Msigma2,
		     const double & Msigma,const double & Gamma_x) {
  return Msigma2/Complex(Msigma2-x,-Msigma*Gamma_x);
}

double RChL::SigmaPhaseSpace(const double & q2,const double & mP2) {
  double arg = 1.-4.*mP2/q2;
  return (arg>0. ? sqrt(arg) : 0.);
}

double RChL::GammaSigma(const double & x,const double & Msigma2,
		       const double & GammaSigma0,const double & mpi2) {
  double sigPoleP = SigmaPhaseSpace(Msigma2,mpi2);
  if (sigPoleP<=0.) return 0.;
  return GammaSigma0 * SigmaPhaseSpace(x,mpi2)/sigPoleP;
}

double RChL::Fsigma(const double & q2,const double & x,const double & mpi2,
		   const double & Rsigma2) {
  double lam = Tools::Lambda(q2,x,mpi2);
  return exp(-lam*Rsigma2/(8.*q2));
}

/////////////////////////////////////////////////////////////////////////
// Two-meson loop function A_PQ(s), Eq.(45)-(49) of 1203.3955.
// See the KNOWN LIMITATION note in the header: Jtilde_PQ(s) is
// approximated as Jbar_PQ(s) (dropping the -s*Jbar'(0) term).
/////////////////////////////////////////////////////////////////////////
Complex RChL::A_PQ(const double & s,const double & mP2,const double & mQ2,
		  const double & F,const double & mu2) {
  double Sigma = mP2+mQ2, Delta = mP2-mQ2;
  double muP = mP2/(32.*sqr(M_PI)*sqr(F)) * log(mP2/mu2);
  double muQ = mQ2/(32.*sqr(M_PI)*sqr(F)) * log(mQ2/mu2);
  double kPQ = (fabs(Delta)>1.e-12 ? sqr(F)/Delta * (muP-muQ) : 0.);
  // nu^2 = lambda(s,mP2,mQ2); below threshold this is negative, so we
  // continue nu into the complex plane rather than discarding its
  // imaginary part (cf. the paper's explicit warning on this point).
  double nu2 = Tools::Lambda(s,mP2,mQ2);
  Complex nu = (nu2>=0. ? Complex(sqrt(nu2),0.) : Complex(0.,sqrt(-nu2)));
  Complex logterm(0.,0.);
  if (abs(nu)>1.e-12) {
    Complex numer = (Complex(s,0.)+nu)*(Complex(s,0.)+nu) - Complex(sqr(Delta),0.);
    Complex denomr= (Complex(s,0.)-nu)*(Complex(s,0.)-nu) - Complex(sqr(Delta),0.);
    logterm = log(numer/denomr);
  }
  Complex Jbar = 1./(32.*sqr(M_PI)) *
    ( 2. + (Delta/s - Sigma/Delta) * log(mQ2/mP2) - (nu/s) * logterm );
  Complex Jtilde = Jbar; // approximation, see header note
  Complex M = (1./(12.*s))*(s-2.*Sigma)*Jbar +
              (sqr(Delta)/(3.*sqr(s)))*Jtilde -
              (1./6.)*kPQ + 1./(288.*sqr(M_PI));
  Complex L = sqr(Delta)/(4.*s) * Jbar;
  return -192.*sqr(M_PI) * (s*M - L) / s;
}

double RChL::Gamma_a1_PionFit(const double & q2,const double & Mrho,
			     const double & mpi2) {
  static const double a=1.54712, b=3.83256, c=4.52798, d=0.30997,
                       e=1.56106, f=3.73605, g=2.00856, h=0.38688,
                       p=-0.00108;
  double thr = 9.*mpi2;
  if (q2<=thr) return 0.;
  double MpiSum2 = sqr(Mrho+sqrt(mpi2));
  if (q2<MpiSum2) {
    double x = q2-thr;
    return pow(x,3.)*(a - b*x + c*sqr(x));
  }
  else if (q2<3.*MpiSum2) {
    return q2*(d - e/q2 + f/sqr(q2) - g/pow(q2,3.));
  }
  return h + 2.*p*(q2-3.*MpiSum2)/MpiSum2;
}
