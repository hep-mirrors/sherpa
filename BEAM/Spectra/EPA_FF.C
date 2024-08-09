#include "BEAM/Spectra/EPA_FF.H"

#include "ATOOLS/Math/Special_Functions.H"
#include "ATOOLS/Math/Bessel_Integrator.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/Settings.H"

#include <fstream>
#include <string>

using namespace BEAM;
using namespace ATOOLS;
using namespace std;

ATOOLS::Special_Functions ATOOLS::SF;

///////////////////////////////////////////////////////////////////////////////////
//
// Base class for different "form factors" that will enter the EPA spectra.
// By far and large they incorporate any information about internal structure.
//
// When dealing with the impact parameters b not
// - that we obtain the particle radii in fm and that we also use fm for any
//   potential external setting of impact parameters;
// - that we will give back impact parameters in fm as well; but
// - that the impact parameters b and radii R usually show up in conjunction
//   with photon energies, transverse momenta or particle masses, all in GeV.
//   The products enter functions such as exponentials or Bessel functions, 
//   hence we have to divide by hbar*c = 0.197 GeV fm to convert distances into
//   units of 1/GeV.
///////////////////////////////////////////////////////////////////////////////////

// Initialisation of relevant beam parameters: 
// note that the particle radius is in fm and transformed into 1/GeV
EPA_FF_Base::EPA_FF_Base(const ATOOLS::Flavour & beam) :
  m_beam(beam), m_mass(m_beam.Mass(true)), m_mass2(ATOOLS::sqr(m_mass)),
  m_R(m_beam.Radius()/rpa->hBar_c()),
  m_q2min(0.), m_q2max(1.e99), m_pt2max(-1.),
  p_Nred_x(nullptr), p_N_xb(nullptr),
  m_approx(0), m_analytic(0) {}

const double EPA_FF_Base::SelectB(double & x) {
  if (x<m_xmin) x = m_xmin;
  if (x>m_xmin) x = m_xmax;
  double maxval = (*p_max_xb)(x), b = m_bmin;
  size_t trials = 10000;
  while (trials>0) {
    b = m_bmin*pow(m_bmax/m_bmin,ran->Get());
    if ((*p_N_xb)(x,b)/maxval>ran->Get()) break;
  }
  // Return b in fm instead of the internally used 1/GeV.
  return b*rpa->hBar_c();
}

// We assume that the units in the b axis are in 1/GeV
void EPA_FF_Base::Fill_Nxb_Table(axis & xaxis,axis & baxis) {
  msg_Out()<<METHOD<<" in "<<xaxis.m_nbins<<" * "<<baxis.m_nbins<<" bins:\n"
  	   <<"   x in ["<<xaxis.m_xmin<<", "<<xaxis.m_xmax<<"], "
  	   <<"b in ["<<baxis.m_xmin<<", "<<baxis.m_xmax<<"], "
  	   <<"from R = "<<m_R<<" 1/GeV = "<<(m_R*rpa->hBar_c())<<" fm.\n";
  m_xmin   = xaxis.m_xmin; m_xmax   = xaxis.m_xmax;
  m_bmin   = baxis.m_xmin; m_bmax   = baxis.m_xmax;
  p_N_xb   = new TwoDim_Table(xaxis,baxis);
  p_max_xb = new OneDim_Table(xaxis);
  N_xb_int * kernel = new N_xb_int(this);
  Bessel_Integrator bessel(kernel,1);
  for (size_t i=0;i<xaxis.m_nbins;i++) {
    double maxval = 0.;
    for (size_t j=0;j<baxis.m_nbins;j++) {
      kernel->SetXB(xaxis.x(i),baxis.x(j));
      double value = sqr(bessel())/M_PI;
      if (!(i%10) && !(j%10)) {
	SetSwitch("analytic",1);
	double ref = N(xaxis.x(i),baxis.x(j));
	msg_Out()<<METHOD<<"(x = "<<xaxis.x(i)<<", "
		 <<"b = "<<baxis.x(j)<<" 1/GeV = "
		 <<(baxis.x(j)*rpa->hBar_c())<<" fm): num = "<<value<<" vs. "
		 <<"ana = "<<ref<<" [1], ratio = "<<value/ref<<"\n";
	SetSwitch("analytic",0);
      }
      p_N_xb->Fill(i,j,value);
      if (dabs(value)>maxval) maxval = dabs(value); 
    }
    p_max_xb->Fill(i,maxval);
  }
}

double N_xb_int::operator()(double y) {
  // Integration argument y here is bT*qT as mandated by the Bessel function:
  // - argument of form factor Q^2 = qT^2+x^2m^2 with qT^2 = y^2/bT^2
  // - we assume that m_b is in units of 1/GeV, qT is in GeV, overall results
  //   are in GeV.
  double qT   = y/m_b, qT2 = sqr(qT), xm2 = sqr(m_x*p_ff->Mass()); 
  double res = ( 1./m_b * qT2/(qT2+xm2) * (*p_ff)(m_x,qT2+xm2) ); 
  return res;
}

void EPA_FF_Base::Fill_Nredx_Table(axis & xaxis) {
  // The b axis is defined by the particle radius taken in 1/GeV 
  p_Nred_x            = new OneDim_Table(xaxis);
  Nred_x_int * kernel = new Nred_x_int(this);
  Gauss_Integrator gauss(kernel);
  double bmin = m_R, bmax = p_N_xb->GetAxis(1).m_xmax;
  for (size_t i=0;i<xaxis.m_nbins;i++) {
    kernel->SetX(xaxis.x(i));
    double value = gauss.Integrate(bmin,bmax,1.e-5);
    //SetSwitch("analytic",1);
    //double Nred  = ReducedN(xaxis.x(i)); 
    //SetSwitch("analytic",0);
    p_Nred_x->Fill(i,value);
  }
}

double Nred_x_int::operator()(double b) {
  // Result is in 1/GeV, will yield 1/GeV^2 upon integration
  return 2.*M_PI*b*(*p_ff->GetN_xb())(m_x,b);
}


///////////////////////////////////////////////////////////////////////////////////
//
// Point-like form factors and related functions in different approximations, all
// based on Budnev et al., Phys. Rep. C15 (1974) 181.
//
///////////////////////////////////////////////////////////////////////////////////

EPA_Point::EPA_Point(const ATOOLS::Flavour & beam) :
  EPA_FF_Base(beam) { m_analytic = 1; m_approx = 2; }

const double EPA_Point::operator()(const double & x,const double & Q2) {
  double wt = 1.;
  if (m_approx>=1) wt *= (1.+sqr(1.-x))/2.; 
  if (m_approx>=2) wt -= (1.-x)*Q2min(x)/Q2; 
  return wt;
}

const double EPA_Point::N(const double & x) {
  // Budnev et al., Phys. Rep. C15 (1974) 181, Eq. (6.17b)
  // This is in units of [1]
  double q2min = Q2min(x), q2max = Q2max(x);
  if (x<=0. || x>=1. || q2max<=q2min) return 0.;
  double wt            = log(q2max/q2min);
  if (m_approx>=1) wt *= (1.+sqr(1.-x))/2.; 
  if (m_approx>=2) wt -= sqr(m_mass*x) * (1./q2min - 1./q2max);
  return wt;
}

const double EPA_Point::SelectB(const double & x) {
  double bmin = 1.e-12, bmax = 1.e12;
  if (m_R>=1.e-12) bmin = 1.e-6*m_R, bmax = 1.e12*m_R;
  double b = bmin, maxval = N(x,b);
  size_t trials = 10000;
  while (trials>0) {
    b = bmin*pow(bmax/bmin,ran->Get());
    if (N(x,b)/maxval>ran->Get()) break;
  }
  return b*rpa->hBar_c();
}

const double EPA_Point::N(const double & x,const double & b) {
  // Result is in GeV^2, inherited from the mass in the argument.
  // Integration over the (2D) impact parameter plane (in units of 1/GeV^2)
  // will yield a result in units of [1].
  if (m_analytic>0) {
    double arg = m_mass*x*SF.Kn(1,m_mass*x*b);
    return sqr(arg)/M_PI;
  }
  return (*p_N_xb)(x,b);
}

const double EPA_Point::ReducedN(const double & x) {
  // Result is in units of [1].
  if (m_analytic>0) {
    double mxR = m_mass*x*m_R, mxR2 = sqr(mxR);
    double K0  = SF.Kn(0,mxR), K1 = SF.Kn(1,mxR);
    return mxR2*(sqr(K0)-sqr(K1)) + 2.*mxR*K0*K1;
  }
  return (*p_Nred_x)(x);
}

///////////////////////////////////////////////////////////////////////////////////
//
// Dipole form factors and related functions in different approximations, all
// based on Budnev et al., Phys. Rep. C15 (1974) 181.
// In its present form, this form factor makes sense for nucleons only.
//
///////////////////////////////////////////////////////////////////////////////////

EPA_Dipole::EPA_Dipole(const ATOOLS::Flavour & beam) :
  EPA_FF_Base(beam), m_Lambda2(0.71) {
  m_approx = 1; 
  if (!m_beam.IsNucleon()) THROW(fatal_error,"Wrong form factor for "+m_beam.IDName());
  m_mu2 = m_beam.Charge()!=0. ? 2.79*2.79 : 0.;
  SetABC();
}

void EPA_Dipole::FillTables(const size_t & nxbins,const size_t & nbbins) {
  size_t nx = nxbins==0 ? 120 : nxbins;
  size_t nb = nbbins==0 ? 150 : nbbins;
  axis xaxis(nx,     1.e-6,        1., axis_mode::log); 
  axis baxis(nb, 1.e-3*m_R, 1.e12*m_R, axis_mode::log);
  Fill_Nxb_Table(xaxis,baxis);
  Fill_Nredx_Table(xaxis);
}
  
const double EPA_Dipole::operator()(const double & x,const double & Q2) {
  // c.f. V.M. Budnev et al., Phys. Rep. C15(1974)181, Eq. (D.7)
  double q2min = Q2min(x), q2max = Q2max(x), prefC = m_mu2, prefD = 1.;
  // taking into account Q^2-dependence of form factors by over-riding their
  // value at Q^2 = 0.
  if (m_approx>0) {
    prefC = m_mu2/sqr(1.+Q2/m_Lambda2);
    prefD = (4.*m_mass2+Q2*m_mu2)/(4.*m_mass2+Q2)/sqr(1.+Q2/m_Lambda2);
  }
  return ( sqr(x)/2.*prefC + (1.-x)*(1.+q2min/Q2)*prefD );
}

void EPA_Dipole::SetABC() {
  // a, b, c coeffients from Budnev et al., Eq. (D.7)
  m_aDip = (1.+m_mu2)/4. + 4.*m_mass2/m_Lambda2;  // should be  7.16
  m_bDip = 1.-4.*m_mass2/m_Lambda2;               // should be -3.96
  m_cDip = (m_mu2-1.)/pow(m_bDip,4.);             // should be  0.028
  msg_Out()<<METHOD<<": a = "<<m_aDip<<", b = "<<m_bDip<<", c = "<<m_cDip<<"\n";
}

const double EPA_Dipole::phi(const double & y,const double & arg) {
  // phi_i(x) from Budnev et al., Eq. (D.7)
  double one_arg = 1.+arg, sum1=0., sum2=0.;
  sum1 += 1./one_arg + 1./(2.*one_arg*one_arg) + 1./(3.*one_arg*one_arg*one_arg);
  sum2 += ( m_bDip/one_arg + m_bDip*m_bDip/(2.*one_arg*one_arg) +
	    m_bDip*m_bDip*m_bDip/(3.*one_arg*one_arg*one_arg) );
  return ((1.+m_aDip*y)*(-log(1.+1./arg) + sum1) -
	  (1.-m_bDip)*y/(4.*arg*pow(1.+arg,3)) +
	  m_cDip*(1.+y/4.)*(log((1.+arg-m_bDip)/(1.+arg))+sum2) );
}

const double EPA_Dipole::N(const double & x) {
  // c.f. V.M. Budnev et al., Phys. Rep. C15(1974)181, Eq. (D.7)
  double q2min = Q2min(x), q2max = Q2max(x);
  if (q2max<=q2min) return 0.;
  // taking into account the Q^2-dependence of form factors ...
  if (m_approx>0) {
    double y = sqr(x)/(1.-x);
    return (1.-x)*(phi(y,q2max/m_Lambda2)-phi(y,q2min/m_Lambda2)); 
  }
  // ... or ignoring it.
  return ( (1.-x+m_mu2*sqr(x)/2.)*log(q2max/q2min) -
	   (1.-x)*(1.-q2min/q2max) );
}

const double EPA_Dipole::N(const double & x,const double & b) {
  if (m_analytic>0) return -1.;
  return (*p_N_xb)(x,b);
}

const double EPA_Dipole::ReducedN(const double & x) {
  if (m_analytic>0) return -1.;
  return (*p_Nred_x)(x);
}



///////////////////////////////////////////////////////////////////////////////////
//
// Gaussian form factors replacing the dipole form factors of Budnev et al.,
// Phys. Rep. C15 (1974) 181.
//
///////////////////////////////////////////////////////////////////////////////////

EPA_Gauss::EPA_Gauss(const ATOOLS::Flavour & beam) : EPA_FF_Base(beam), m_Q02(1.) {
  m_approx = 1;
  if (m_beam==ATOOLS::Flavour(kf_p_plus)) { m_mu2 = 2.79*2.79; m_Q02 = 0.71; }
  else if (m_beam.IsIon()) { m_mu2 = 0.; m_Q02 = sqr(2./m_R); m_pt2max = sqr(1./m_R); }
}

void EPA_Gauss::FillTables(const size_t & nxbins,const size_t & nbbins) {
  size_t nx   = nxbins==0 ? 120 : nxbins;
  size_t nb   = nbbins==0 ? 150 : nbbins;
  double xmin = 1.e-6/m_beam.GetAtomicNumber(), xmax = 1./m_beam.GetAtomicNumber();
  double bmin = 1.e-3*m_R, bmax = 1.*m_R;
  axis xaxis(nx, xmin, xmax, axis_mode::log); 
  axis baxis(nb, bmin, bmax, axis_mode::log);
  Fill_Nxb_Table(xaxis,baxis);
  Fill_Nredx_Table(xaxis);
}
  
const double EPA_Gauss::operator()(const double & x,const double & Q2) {
  // c.f. V.M. Budnev et al., Phys. Rep. C15(1974)181, Eq. (D.7)
  // but modifying the dipole form to a Gaussian form
  double q2min = Q2min(x), q2max = Q2max(x), prefC = m_mu2, prefD = 1.;
  // taking into account Q^2-dependence of form factors by over-riding their
  // value at Q^2 = 0 ...
  if (m_approx>0) {
    prefC = m_mu2*exp(-Q2/m_Q02);
    prefD = exp(-Q2/m_Q02);
  }
  // ... or ignoring it
  return ( sqr(x)/2.*prefC + (1.-x)*(1.+q2min/Q2)*prefD);
}

const double EPA_Gauss::N(const double & x) {
  // c.f. V.M. Budnev et al., Phys. Rep. C15(1974)181, Eq. (D.7)
  double q2min = Q2min(x), q2max = Q2max(x), term1, term2;
  // taking into account the Q^2-dependence of form factors and using that
  // Ei(-x) = -IncompleteGamma(0,x)
  if (m_approx>0) {
    term1 = SF.IncompleteGamma(0,q2min/m_Q02) - SF.IncompleteGamma(0,q2max/m_Q02);
    term2 = term1 -q2min * ( exp(-q2max/m_Q02)/q2max - exp(-q2min/m_Q02)/q2min +
			     1./m_Q02 * term1 );
  }
  // ... or ignoring it.
  else {
    term1 = log(q2max/q2min);
    term2 = term1+q2min*(1./q2max-1./q2min);
  }
  double res = m_mu2*sqr(x)*term1 + (1.-x)*term2; 
  msg_Out()<<"\n"<<METHOD<<"("<<m_approx<<", x = "<<x<<" "
	   <<"Q^2 in ["<<q2min<<", "<<q2max<<"] :"
  	   <<term1<<" & "<<term2<<" for Q0 = "<<sqrt(m_Q02)<<": "<<res<<".\n";
  return res;
}

const double EPA_Gauss::N(const double & x,const double & b) {
  if (m_analytic>0) return -1.;
  return (*p_N_xb)(x,b);
}

const double EPA_Gauss::ReducedN(const double & x) {
  if (m_analytic>0) return -1.;
  return (*p_Nred_x)(x);
}

///////////////////////////////////////////////////////////////////////////////////
//
// Form factor for ions in a Wood-Saxon potential.  A few comments:
//
// - Output for the density modifed to recover the atomic number A on integration,
//   and in units of fm^-3 or GeV^3.   Internally we normalise it to unity,
//   because we will later multiply the square of the form factor with the square
//   of the charge Z.
// - Radius in 1/GeV, therefore "nucelar skin" m_d also in 1/GeV
//
///////////////////////////////////////////////////////////////////////////////////

EPA_WoodSaxon::EPA_WoodSaxon(const ATOOLS::Flavour & beam) :
  EPA_FF_Base(beam), m_d(0.5/rpa->hBar_c()), p_FF_Q2(nullptr), p_N(nullptr) 
{
  if (!m_beam.IsIon()) THROW(fatal_error,"Wrong form factor for "+m_beam.IDName());
  m_pt2max = sqr(1./m_R);
}
  
void EPA_WoodSaxon::FillTables(const size_t & nxbins,const size_t & nbbins) {
  m_rho0 = CalculateDensity();
  InitFFTable(1.e-12,1.e4);
  InitNTable(1.e-10,1.);
  size_t nx   = nxbins==0 ? 120 : nxbins;
  size_t nb   = nbbins==0 ? 150 : nbbins;
  double xmin = 1.e-6/m_beam.GetAtomicNumber(), xmax = 1./m_beam.GetAtomicNumber();
  double bmin = 1.e-3*m_R, bmax = 1.*m_R;
  axis xaxis(nx, xmin, xmax, axis_mode::log); 
  axis baxis(nb, bmin, bmax, axis_mode::log);
  Fill_Nxb_Table(xaxis,baxis);
  Fill_Nredx_Table(xaxis);
}
  
void EPA_WoodSaxon::InitFFTable(const double & q2min,const double & q2max) {
  p_FF_Q2 = new OneDim_Table(axis(100000,q2min,q2max,axis_mode::log));
  WS_potential * ws = new WS_potential(m_R,m_d);
  Gauss_Integrator gauss(ws);
  for (size_t i=0;i<p_FF_Q2->GetAxis().m_nbins;i++) {
    ws->SetQ(sqrt(p_FF_Q2->GetAxis().x(i)));
    double rmin = 0., rmax = m_R;
    double res  = m_rho0*gauss.Integrate(rmin,rmax,1.e-6,0), inc = 0.;
    do {
      rmin = rmax; rmax *= 2.;
      res += inc  = m_rho0*gauss.Integrate(rmin,rmax,1.e-6,0);
    } while (rmax<4.*m_R || dabs(inc/res)>1.e-6);
    p_FF_Q2->Fill(i,res);
  }
}

void EPA_WoodSaxon::InitNTable(const double & xmin,const double & xmax) {
  p_N = new OneDim_Table(axis(10000,xmin,xmax,axis_mode::log));
  N_argument * nx = new N_argument(this);
  Gauss_Integrator gauss(nx);
  for (size_t i=0;i<p_N->GetAxis().m_nbins;i++) {
    double x     = p_N->GetAxis().x(i);
    double q2min = Q2min(x); 
    double q2max = Q2min(x)+m_pt2max; 
    nx->SetX(x);
    double res   = gauss.Integrate(q2min,q2max,1.e-6,0);
    if (!(i%1000)) msg_Out()<<METHOD<<"("<<q2min<<", "<<q2max<<"): N("<<x<<") = "<<res<<"\n";
    p_N->Fill(i,res);
  }
}

double EPA_WoodSaxon::CalculateDensity() {
  Rho_argument * rho = new Rho_argument(m_R,m_d);
  Gauss_Integrator gauss(rho);
  double rmin = 0., rmax = m_R;
  double res  = gauss.Integrate(rmin,rmax,1.e-3), inc = 0.;
  do {
    rmin = rmax; rmax *= 2.;
    res += inc  = gauss.Integrate(rmin,rmax,1.e-6);
  } while (inc>1.e-6*res);
  return 1./res;
}

const double EPA_WoodSaxon::operator()(const double & x,const double & Q2) {
  return (1.-x)*(1.-Q2min(x)/Q2)*(*p_FF_Q2)(Q2);
}

const double EPA_WoodSaxon::N(const double & x) {
  if (m_analytic>0) return -1.;
  return (*p_N)(x);
}

const double EPA_WoodSaxon::N(const double & x,const double & b) {
  if (m_analytic>0) return -1.;
  return (*p_N_xb)(x,b);
}

const double EPA_WoodSaxon::ReducedN(const double & x) {
  if (m_analytic>0) return -1.;
  return (*p_Nred_x)(x);
}
