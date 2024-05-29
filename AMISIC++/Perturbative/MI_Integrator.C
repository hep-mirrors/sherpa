#include "AMISIC++/Perturbative/MI_Integrator.H"
#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Tools/Lookup_Tables.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMISIC;
using namespace ATOOLS;

MI_Integrator::MI_Integrator(MI_Processes * procs) :
  p_procs(procs), m_pt2min(1.), m_bmax(0.), m_bvol(1.),
  m_xsmax(0.), m_xsmaxMO(0.), 
  m_MCpoints(10000) { }

void MI_Integrator::Initialize(PDF::ISR_Handler * isr) {
  ////////////////////////////////////////////////////////////////////////////
  // Minimal pt2, maximal b, and minimal and maximal x from the PDFs
  ////////////////////////////////////////////////////////////////////////////
  m_pt2min    = p_procs->PT2Min();
  m_Emin      = p_procs->EMin();
  m_bmax      = p_procs->Bmax();
  for (size_t i=0;i<2;i++) {
    m_xmin[i] = Max(   1.e-6,isr->PDF(i)->XMin());
    m_xmax[i] = Min(1.-1.e-6,isr->PDF(i)->XMax());
  }
  ////////////////////////////////////////////////////////////////////////////
  // 10000 points will yield errors of about 1%.
  ////////////////////////////////////////////////////////////////////////////
  m_MCpoints = ATOOLS::Min(size_t(10000),size_t((*mipars)("nMC_points")));
}

double MI_Integrator::TrialEvent(const double & s) {
  ////////////////////////////////////////////////////////////////////////////
  // Producing a trial event to be used in Impact_Parameter to define
  // kinematics and impact parameter in the first scatter for impact-parameter
  // dependent matter overlaps.  The impact parameter dependence is captured
  // by the matter overlap and in the explicitly impact-parameter dependent
  // Sudakov argument.
  //////////////////////////////////////////////////////////////////////////// 
  if (4.*m_pt2min>=(1.-1.e-6)*s) return 0.;
  double invpt2min  = 1./m_pt2min, invpt2max = 4./s;
  double invpt2diff = invpt2min-invpt2max;
  double wt = 0.;
  do {
    m_pt2 = 1./(invpt2min-ran->Get()*invpt2diff);
    MakeKinematics(m_pt2,s);
    wt    = (*p_procs)(false) * m_yvol * sqr(m_pt2) * invpt2diff /m_xsmax;
  } while (wt<ran->Get());
  return m_pt2;
}

double MI_Integrator::XSecInPT2Bin(const double & pt2,const double & s) {
  ////////////////////////////////////////////////////////////////////////////
  // This method is used to fill the pt^2 look-up table of the Sudakov form
  // factor argument which we use in the selection of the first pt^2 and b
  // in min bias events and in the determination of b in MPI events.
  // 
  // MC-integrated differential cross setion dsigma/dpt^2 in 1/GeV^4 for a
  // given transverse momentum pt^2:
  // 1/(16 pi) int_{-ymax}^{+ymax} dy_1 dy_2
  //      [  x_1 f(x_1, muF_1) x_2 f(x_2, muF_2)
  //        |M(shat,that,uhat)|^2 / shat^2 ]
  // Integrating over the impact parameter (usually set externally through the
  // MI_Processes) recovers the original result (the impact-parameter integral
  // over the matter overlap should integrate to unity!).
  // 1/(16 pi) int_{-ymax}^{+ymax} dy_1 dy_2
  //      [  x_1 f(x_1, muF_1) x_2 f(x_2, muF_2)
  //         |M(shat,that,uhat)|^2 / shat^2 O(x1, x2, muF1, muF2) ].
  ////////////////////////////////////////////////////////////////////////////
  if (pt2<=m_pt2min || 4.*pt2>=(1.-1.e-6)*s) return 0.;
  double sum = 0., sum2 = 0., xsec = 0., uncert = 0., xs;
  unsigned long int sumtrials = 0;
  do {
    xs = MakeKinematics(pt2,s) ? (*p_procs)(true) * m_yvol : 0.;
    sum  += xs;
    sum2 += sqr(xs);
    if (++sumtrials%m_MCpoints==0) {
      xsec   = sum/double(sumtrials);
      uncert = sqrt(sum2 - sqr(xsec))/double(sumtrials);
    }
  } while (xsec==0 || (uncert/xsec>5.e-3));
  m_xsec   = sum/double(sumtrials);
  m_uncert = sqrt(sum2 - sqr(m_xsec))/double(sumtrials);
  return m_xsec;
}

double MI_Integrator::
operator()(const double & s,const bool & withMO,const bool & intB) {
  ////////////////////////////////////////////////////////////////////////////
  // MC-integrated total cross setion sigma_hard in 1/GeV^2
  // 1/(16 pi) int dpt^2 int_{-ymax}^{+ymax} dy_1 dy_2  dsigma/dpt^2
  // with dsigma/dpt^2 given as below as
  // dsigma/dpt^2 = x_1 f(x_1, muF_1) x_2 f(x_2, muF_2)
  //                |M(shat,that,uhat)|^2 / shat^2
  // The flag withMO includes the dynamic matter overlap O(x1,x2,muF1,muF2) 
  // as an additional factor, to be integrated over the impact parameter,
  // if intB is set to true.
  // 1/(16 pi) int d^2b int dpt^2 int_{-ymax}^{+ymax} dy_1 dy_2
  //                     [ dsigma/dpt^2 O(x1,x2,muF1,muF2) ]  
  // which should not change the overall result.
  ////////////////////////////////////////////////////////////////////////////
  if (4.*m_pt2min>=(1.-1.e-6)*s) return 0.;
  double invpt2min  = 1./m_pt2min, invpt2max = 4./s;
  double invpt2diff = invpt2min-invpt2max;
  double pt2vol, R02, b2, bvol = 1.;
  ////////////////////////////////////////////////////////////////////////////
  // Maximal radius estimate if integrating over impact parameter.
  ////////////////////////////////////////////////////////////////////////////
  if (withMO && intB) {
    p_procs->GetOverlap()->FixDynamicRadius();
    R02 = p_procs->GetOverlap()->DynamicRadius2();
  }
  double sum = 0., sum2 = 0., xsec = 0., uncert = 0., xs;
  unsigned int sumtrials = 0.;
  do {
    ////////////////////////////////////////////////////////////////////////
    // Select pt^2 according to 1/pt^4 and calculate the Jacobean
    ////////////////////////////////////////////////////////////////////////
    double rand = ran->Get();
    m_pt2  = 1./(invpt2min-ran->Get()*invpt2diff);
    pt2vol = invpt2diff * sqr(m_pt2);
    ////////////////////////////////////////////////////////////////////////
    // Select b^2 and calculate the Jacobean for integrating over
    // impact parameter.
    ////////////////////////////////////////////////////////////////////////
    if (withMO && intB) {
      b2   = -R02*log(Max(1.e-12,ran->Get()));
      bvol = M_PI*R02*exp(b2/R02);
      p_procs->SetB(sqrt(b2));
    }
    ////////////////////////////////////////////////////////////////////////
    // Fix the rest of the kinematics: rapidities, x's, Mandelstams. 
    ////////////////////////////////////////////////////////////////////////
    xs    = (MakeKinematics(m_pt2,s) ?
	     (*p_procs)(withMO) * m_yvol * pt2vol * bvol : 0. );
    sum  += xs;
    sum2 += sqr(xs);
    if (++sumtrials%m_MCpoints==0) {
      xsec   = sum/double(sumtrials);
      uncert = sqrt(sum2 - sqr(xsec))/double(sumtrials);
    }
    if ( withMO && xs>m_xsmaxMO) m_xsmaxMO = xs;
    if (!withMO && xs>m_xsmax)   m_xsmax   = xs;
  } while (xsec==0 || (uncert/xsec>5.e-3));
  m_xsec   = sum/double(sumtrials);
  m_uncert = sqrt(sum2 - sqr(m_xsec))/double(sumtrials);
  return m_xsec;
}
  
bool MI_Integrator::MakeKinematics(const double & pt2,const double & s)
{
  ////////////////////////////////////////////////////////////////////////////
  // We select the two rapidities of the two outgoing massless particles
  // flat in the full interval and hit-or-miss by making sure the x values
  // are inside the allowed range:
  //  x_{1,2} = xT/2 * [exp(+/- y1) + exp(+/- y2)] with xT = (4pt^2/S)^0.5.
  ////////////////////////////////////////////////////////////////////////////
  m_pt2     = pt2;
  double xt = sqrt(4.*m_pt2/s);
  if (xt>1.) return false;
  //////////////////////////////////////////////////////////////////////////
  // Generate two trial rapidities and keep the integration volume.
  ////////////////////////////////////////////////////////////////////////////
  double ymax = log(1./xt*(1.+sqrt(1.-xt*xt)));
  for (size_t i=0;i<2;i++) m_y[i] = ymax*(-1.+2.*ran->Get());
  m_yvol = sqr(2.*ymax);
  ////////////////////////////////////////////////////////////////////////////
  // Obtain Bjorken-x from the trial rapidities and make sure they are
  // still smaller than the residual x in the remnants.  Will enter the
  // exact PDF calculation in the exact differential cross section.
  ////////////////////////////////////////////////////////////////////////////
  m_x[0] = xt * (exp( m_y[0]) + exp( m_y[1]))/2.;
  m_x[1] = xt * (exp(-m_y[0]) + exp(-m_y[1]))/2.;
  ////////////////////////////////////////////////////////////////////////////
  // TODO: Misses term for incoming masses, we will hope for the best by
  //       treating them as massless - that could be improved in the future.
  ////////////////////////////////////////////////////////////////////////////
  if (m_x[0]<=m_xmin[0] || m_x[0]>=m_xmax[0] ||
      m_x[1]<=m_xmin[1] || m_x[1]>=m_xmax[1] ||
      m_x[0]*m_x[1]<=xt*xt)      return false; 
  if (m_x[0]*sqrt(s)/2.<=m_Emin ||
      m_x[1]*sqrt(s)/2.<=m_Emin) return false;
  ////////////////////////////////////////////////////////////////////////////
  // Mandelstams for the exact matrix element.
  ////////////////////////////////////////////////////////////////////////////
  m_cost = sqrt(1.-Min(1.,(xt*xt)/(m_x[0]*m_x[1])));
  m_shat = m_x[0] * m_x[1] * s;
  m_that = -0.5 * m_shat * (1.-m_cost);
  m_uhat = -0.5 * m_shat * (1.+m_cost);
  return true;
}
