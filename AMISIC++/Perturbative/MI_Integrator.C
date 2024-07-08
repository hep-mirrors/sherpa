#include "AMISIC++/Perturbative/MI_Integrator.H"
#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Tools/Matter_Overlap.H"
#include "AMISIC++/Tools/Lookup_Tables.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMISIC;
using namespace ATOOLS;

MI_Integrator::MI_Integrator(MI_Processes * procs) :
  p_procs(procs), m_pt2min(1.), m_xsmax(-1.),
  m_MCpoints(10000) { }

void MI_Integrator::Initialize(PDF::ISR_Handler * isr) {
  ////////////////////////////////////////////////////////////////////////////
  // Minimal pt2, maximal b, and minimal and maximal x from the PDFs
  ////////////////////////////////////////////////////////////////////////////
  m_pt2min    = p_procs->PT2Min();
  m_Emin      = p_procs->EMin();
  for (size_t i=0;i<2;i++) {
    m_xmin[i] = Max(   1.e-6,isr->PDF(i)->XMin());
    m_xmax[i] = Min(1.-1.e-6,isr->PDF(i)->XMax());
  }
  ////////////////////////////////////////////////////////////////////////////
  // 10000 points will yield errors of about 1%.
  ////////////////////////////////////////////////////////////////////////////
  m_MCpoints = ATOOLS::Min(size_t(10000),size_t((*mipars)("nMC_points")));
}

bool MI_Integrator::TrialEvent(const double & s,Matter_Overlap * mo) {
  ////////////////////////////////////////////////////////////////////////////
  // Producing a trial event to define kinematics and impact parameter in
  // the first scatter.  The impact parameter dependence is captured
  // by the matter overlap and in the explicitly impact-parameter dependent
  // Sudakov rejection algorithm.
  ////////////////////////////////////////////////////////////////////////////
  if (4.*m_pt2min>=(1.-1.e-6)*s) return false;
  double invpt2min = 1./m_pt2min, invpt2max = 4./s;
  double pt2vol    = invpt2min-invpt2max;
  double bwt, wt;
  size_t trials    = 0;
  do {
    m_pt2 = 1./(invpt2min-ran->Get()*pt2vol);
    if (MakeKinematics(m_pt2,s)) 
      wt  = (*p_procs)()*m_yvol*sqr(m_pt2)*pt2vol/m_xsmax;
    else wt = 0.;
    if (trials>100) return false;
  } while (wt<ran->Get());
  return true;
}

double MI_Integrator::
operator()(const double & s,Matter_Overlap * mo,const double & b) {
  ////////////////////////////////////////////////////////////////////////////
  // MC-integrated total cross setion sigma_hard in 1/GeV^2
  // 1/(16 pi) int dpt^2 int_{-ymax}^{+ymax} dy_1 dy_2  dsigma/dpt^2
  // with dsigma/dpt^2 given as below as
  // dsigma/dpt^2 = x_1 f(x_1, muF_1) x_2 f(x_2, muF_2)
  //                |M(shat,that,uhat)|^2 / shat^2
  ////////////////////////////////////////////////////////////////////////////
  if (4.*m_pt2min>=(1.-1.e-6)*s) return 0.;
  double invpt2min = 1./m_pt2min, invpt2max = 4./s;
  double pt2vol    = invpt2min-invpt2max;
  double sum  = 0., sum2 = 0., xsec = 0., uncert = 1.e12, xs;
  double mowt = mo ? (mo->IsDynamic() ? 0. : (*mo)(b)) : 1.;
  unsigned int sumtrials = 0.;
  do {
    ////////////////////////////////////////////////////////////////////////
    // Select pt^2 according to 1/pt^4 and calculate the Jacobean
    ////////////////////////////////////////////////////////////////////////
    m_pt2  = 1./(invpt2min-ran->Get()*pt2vol);
    ////////////////////////////////////////////////////////////////////////
    // Fix the rest of the kinematics: rapidities, x's, Mandelstams. 
    ////////////////////////////////////////////////////////////////////////
    if (MakeKinematics(m_pt2,s)) {
      if (mo && mo->IsDynamic()) {
	mo->FixDynamicRadius(m_x[0],m_x[1]);
	mowt = (*mo)(b);
      }
      xs = (*p_procs)() * m_yvol * pt2vol * sqr(m_pt2);
      if (xs>m_xsmax) m_xsmax = xs;
      sum  += xs * mowt;
      sum2 += sqr(xs * mowt);
    }
    if (++sumtrials%(100*m_MCpoints)==0) {
      xsec   = sum/double(sumtrials);
      uncert = sqrt(sum2 - sqr(xsec))/double(sumtrials);
    }
  } while (xsec==0 || (uncert/xsec>5.e-3));
  m_xsec   = sum/double(sumtrials);
  m_uncert = sqrt(sum2 - sqr(m_xsec))/double(sumtrials);
  return m_xsec;
}


bool MI_Integrator::
MakeKinematics(const double & pt2,const double & s) {
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
      m_x[0]*m_x[1]<=xt*xt) return false;
  if (m_x[0]*sqrt(s)/2.<=m_Emin || m_x[1]*sqrt(s)/2.<=m_Emin) return false;
  ////////////////////////////////////////////////////////////////////////////
  // Mandelstams for the exact matrix element.
  ////////////////////////////////////////////////////////////////////////////
  m_cost = sqrt(1.-Min(1.,(xt*xt)/(m_x[0]*m_x[1])));
  m_shat = m_x[0] * m_x[1] * s;
  m_that = -0.5 * m_shat * (1.-m_cost);
  m_uhat = -0.5 * m_shat * (1.+m_cost);
  return true;
}
