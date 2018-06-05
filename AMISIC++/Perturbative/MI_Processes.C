#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "EXTRA_XS/Main/Single_Process.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/Message.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;


MI_Processes::MI_Processes() : ME_Generator_Base("Amisic::Processes"),
			       m_sigmaND(1.), m_integral(0.) {}

MI_Processes::~MI_Processes() {
  while (!m_groups.empty()) {
    delete m_groups.back();
    m_groups.pop_back();
  }
  m_groups.clear();
}

bool MI_Processes::Initialize(const std::string &path,const std::string &file,
			      MODEL::Model_Base *const model,
			      BEAM::Beam_Spectra_Handler *const beam,
			      PDF::ISR_Handler *const isr) {
  // Get PDFs and couplings
  p_model    = model;
  p_isr      = isr;
  p_pdf[0]   = p_isr->PDF(0);
  p_pdf[1]   = p_isr->PDF(1);
  p_alphaS   = dynamic_cast<MODEL::Running_AlphaS *>
    (model->GetScalarFunction("alpha_S"));
  p_alpha    = dynamic_cast<MODEL::Running_AlphaQED *>
    (model->GetScalarFunction("alpha_QED"));
  // Initialize model parameters:
  // - pt_0, the IR regulator in the propagator and in the strong coupling
  // - pt_min, the IR cut-off for the 2->2 scatters
  // - Ecms, the cms energy of the hadron collision
  m_pt0      = (*mipars)("pt_0");   m_pt02   = m_pt0*m_pt0;
  m_ptmin    = (*mipars)("pt_min"); m_ptmin2 = m_ptmin*m_ptmin;
  m_ecms     = rpa->gen.Ecms();     m_S      = m_ecms*m_ecms;
  // will have to make this part of an external input scheme
  m_scale_scheme   = "MPI"; 
  m_kfactor_scheme = "MPI";
  // These parameters are for the hard part of the Sudakov form factor.
  // It is given by int_{pt^2}^{s/4} dpt^2 dsigma/dpt^2, and we tabulate the
  // integral in nbins between s/4 and pt_min^2, resulting in a stepsize of
  // pt^2_step.  In each of the bins we MC integrate over the rapidities of the
  // two outgoing particles with MC_points points.  This trick is borrowed
  // from Sjostrand's implementation, as explicit integration will be too slow.
  // 10000 points will yield errors of about 1%.
  m_nbins    = int((*mipars)("nPT_bins"));
  m_MCpoints = int((*mipars)("nMC_points"));
  m_intbins.resize(m_nbins);
  m_pt2step  = log(m_S/(4.*m_ptmin2))/double(m_nbins);
  // Mass scheme for the subsequent parton shower.
  m_massmode       = 1;
  SetPSMasses(p_defaultreader);
  // Now initialize the 2->2 scatters and prepare the integral for the
  // "Sudakov form factor", Eq. (37) of Sjostrand-van der Zijl
  return (InitializeAllProcesses() && PrepareSudakovFactor());
}

bool MI_Processes::InitializeAllProcesses() {
  // The following pure QCD processes should be included by default
  // - gg-initial states:  gg->gg, gg->qqb
  // - qq', qbqb', qqb'-initial states: qq'->qq' and friends
  // - qqb initial states: qqb -> qqb, qqb->gg
  // - qq, qbqb initial states: qq->qq, qbqb->qbqb
  // - gq, gqb initial states: gq->gq, gqb->gqb
  // - qg, qbg initial states: qg->qg, qbg->qbg
  // They are clustered according to their initial states to optimize
  // the calls to PDFs and to make sure we have to evaluate the |ME|^2
  // that depend on the partons only once
  m_groups.push_back(new MI_GG_Processes());
  m_groups.push_back(new MI_Q1Q2_Processes());
  m_groups.push_back(new MI_QQB_Processes());
  m_groups.push_back(new MI_QQ_Processes());
  m_groups.push_back(new MI_QG_Processes());
  // The following processes should depend on switches.  At the moment we
  // just add them without further ado.
  // - qg-initiated photon production
  // - qqbar-initiated photon production
  // We are missing di-photon production here: 
  // - gg->gamma gamma and qqbar->gamma gamma
  m_groups.push_back(new MI_QG_QGamma_Processes());
  m_groups.push_back(new MI_QQ_GGamma_Processes());
  // We are missing the production of (heavy quarkonia) mesons MQQ in  
  // - singlet production qqbar -> MQQ, gg -> MQQ, gq -> MQQ+q etc.
  // - octet production gg -> MQQ^(8)+g, qqbar->MQQ^(8)+g etc.
  // We could also add production of gauge bosons:
  // - qqbar->Z, qqbar'->W
  // - gq->Zq, gq->Wq', etc.
  SetPDFs();
  SetAlphaS();
  return true;
}

void MI_Processes::SetPDFs() {
  // The PDFs for the process groups.
  for (list<MI_Process_Group *>::iterator mig = m_groups.begin();
       mig!=m_groups.end();mig++)
    (*mig)->SetPDFs(p_pdf[0],p_pdf[1]);
}

void MI_Processes::SetAlphaS() {
  // The couplings for the process groups.
  for (list<MI_Process_Group *>::iterator mig = m_groups.begin();
       mig!=m_groups.end();mig++) {
    (*mig)->SetAlphaS(p_alphaS);
    (*mig)->SetAlpha(p_alpha);
  }
}

void MI_Processes::CalcPDFs(const double & x1,const double & x2,
			    const double & scale) {
  // Calculate both sets of PDFs at the relevant x and Q^2
  p_pdf[0]->Calculate(x1,Max(scale,p_pdf[0]->Q2Min()));
  p_pdf[1]->Calculate(x2,Max(scale,p_pdf[1]->Q2Min()));
}

const double MI_Processes::XSec(const Vec4D_Vector & momenta) {
  // Return the total parton-level scattering cross section, summed over all
  // contributing processes.  In contrast to the other XSec method the input
  // here are the momenta and not the Mandelstam variables, and the PDFs are
  // calculated here.
  double shat = (momenta[0]+momenta[1]).Abs2();
  double that = (momenta[0]-momenta[2]).Abs2();
  double uhat = (momenta[0]-momenta[3]).Abs2();
  double x1   = 2.*momenta[0][0]/m_ecms, x2 = 2.*momenta[1][0]/m_ecms;  
  double pt2  = that*uhat/shat, scale = pt2;
  CalcPDFs(x1,x2,scale);
  return (*this)(shat,that,uhat);
}

const double MI_Processes::XSec(const double & shat,const double & that,
				const double & uhat) {
  // Return the total parton-level scattering cross section, summed over all
  // contributing processes.  This implicitly assumes that the PDFs have already
  // been set.  It seems this could be a doubling of methods - will have to check.
  return (*this)(shat,that,uhat);
}

const double MI_Processes::operator()(const double & shat,const double & that,
				      const double & uhat) {
  // Return the total parton-level scattering cross section, summed over all
  // contributing processes.  This implicitly assumes that the PDFs have already
  // been set.  
  m_lastxs = 0.;
  double pt2 = that*uhat/shat;
  for (list<MI_Process_Group *>::iterator mig = m_groups.begin();
       mig!=m_groups.end();mig++) {
    (*mig)->SetScale(pt2);
    m_lastxs += (**mig)(shat,that,uhat);
  }
  return m_lastxs;  
}

MI_Process * MI_Processes::SelectProcess() {
  // Sum over all cross sections of all groups and select one of the groups.
  // Then select one of the processes within the group.
  double diff = m_lastxs * ran->Get();
  list<MI_Process_Group *>::iterator mig = m_groups.begin();
  while (mig!=m_groups.end()) {
    diff -= (*mig)->LastXS();
    if (diff<0.) break;
    mig++;
  }
  if (mig==m_groups.end()) mig--;
  return (*mig)->SelectProcess();
}

bool MI_Processes::PrepareSudakovFactor() {
  // In this method we bin pt^2 in nbins, distributed logarithmically
  // and accumulate the integral in steps, the result is being stored in intbins.
  // intbins[i] = Sum_{pt_i^2}^{pt_nmax^2} dSigma(pt_i^2)/dpt^2 * pt_{i+1}^2-pt_i^2
  // where pt_nmax^2 = s/4, the maximal pt^2.
  // N.B.: I use left steps, thereby somewhat overestimating the integral, this
  // could be improved by going trapezoid or similar.
  double pt2last = m_ptmin2*exp(m_pt2step*m_nbins), pt2, dpt2;
  for (int bin=m_nbins-1;bin>=0;bin--) {
    pt2            = m_ptmin2*exp(m_pt2step*bin);
    dpt2           = pt2last-pt2;
    m_intbins[bin] = m_integral += dSigma(pt2) * dpt2/m_sigmaND;
    pt2last        = pt2;
  }
  m_integral *= m_sigmaND;
  msg_Info()<<METHOD<<" calculates integral for Sudakov form factor starting at pt = "
	    <<sqrt(pt2last)<<" in "<<m_nbins<<" steps,\n   sigma = "<<m_integral
	    <<" 1/Gev^2 = "<<(m_integral*rpa->Picobarn()/1.e9)<<" mb.\n";
  return true;
}

double MI_Processes::dSigma(const double & pt2) {
  // Estimated cross setion for a given transverse momentum:
  // It is given by 
  // 1/(16 pi) int_{-ymax}^{+ymax} dy_1 dy_2  [  x_1 f(x_1, pt^2) x_2 f(x_2, pt^2)
  //                                             |M(shat,that,uhat)|^2 / shat^2    ]
  // Here we extracted the factor g^4 = alpha_S^2(pt^2)/16 pi^2 put of the matrix element
  // such that the prefactor is pi instead of 1/(16 pi).
  // We select the two rapidities of the two outgoing massless particles flat in the
  // full interval and hit-or-miss by making sure the x values are inside the allowed
  // range.  x_{1,2} = xT/2 * [exp(+/- y1) + exp(+/- y2)] with xT = (4pt^2/S)^0.5.
  double xt       = sqrt(4.*pt2/m_S);
  double ymax     = log(1./xt*(1.+sqrt(1.-xt*xt)));
  double PSfac    = sqr(2.*ymax);
  double res      = 0., res2 = 0.;
  for (size_t i=0;i<m_MCpoints;i++) {
    double y1     = ymax*(-1.+2.*ran->Get());
    double y2     = ymax*(-1.+2.*ran->Get());
    double x1     = xt * (exp(y1)  + exp(y2))/2.;
    double x2     = xt * (exp(-y1) + exp(-y2))/2.;
    if (x1<1.e-6 || x1>1. || x2<1.e-6 || x2>1.) continue;
    double cost   = sqrt(1.-(xt*xt)/(x1*x2));
    double shat   = x1 * x2 * m_S;
    double that   = -0.5 * shat * (1.-cost);
    double uhat   = -0.5 * shat * (1.+cost);
    CalcPDFs(x1,x2,pt2);
    double dsigma = (*this)(shat,that,uhat) * PSfac;
    res  += dsigma;
    res2 += dsigma*dsigma;
  }
  double result = res/double(m_MCpoints);
  //double uncert = sqrt((res2/double(m_MCpoints)-sqr(result))/double(m_MCpoints));
  return result;
}

const double MI_Processes::SudakovArgument(const double & pt2) const {
  // Linear interpolation between the pre-calculated points for the Sudakov form factor
  int bin     = int(1./m_pt2step*log(pt2/m_ptmin2));
  double pt21 = m_ptmin2*exp(m_pt2step*(bin)), pt22 = m_ptmin2*exp(m_pt2step*(bin+1));
  double val1 = m_intbins[bin],                val2 = m_intbins[bin+1];
  double val  = (val1*(pt22-pt2)+val2*(pt2-pt21))/(pt22-pt21);
  return val;
}

/*
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Function_Base.H"

class SudakovTestArgument: public ATOOLS::Function_Base {
  MI_Processes * p_procs;
public:
  SudakovTestArgument(MI_Processes * procs) : p_procs(procs) {}
  double operator()(double pt2) { return p_procs->dSigma(pt2); }
};

void MI_Processes::TestSudakovFactor() {
  SudakovTestArgument sudarg(this);
  Gauss_Integrator sudint(&sudarg);
  double pt2(10000.);
  double integral = sudint.Integrate(pt2,m_S/4.,0.003,1);
  msg_Out()<<METHOD<<"("<<pt2<<", "<<(m_S/4.)<<") = "<<SudakovArgument(pt2)<<".\n";
}
*/

