#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "EXTRA_XS/Main/Single_Process.H"
#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;


MI_Processes::MI_Processes(const bool & variable_s) :
  ME_Generator_Base("Amisic::Processes"), m_variable_s(variable_s) {}

MI_Processes::~MI_Processes() {
  while (!m_groups.empty()) {
    delete m_groups.back();
    m_groups.pop_back();
  }
  m_groups.clear();
}

bool MI_Processes::Initialize(MODEL::Model_Base *const model,
			      BEAM::Beam_Spectra_Handler *const beam,
			      PDF::ISR_Handler *const isr) {
  // Get PDFs and couplings
  p_isr       = isr;
  m_muFfac    = sqr((*mipars)("FacScale_Factor"));
  for (size_t i=0;i<2;i++) {
    p_pdf[i]  = p_isr->PDF(i);
    m_xmin[i] = Max(1.e-6,p_pdf[i]->XMin());
  }
  p_alphaS    = dynamic_cast<MODEL::Running_AlphaS *>
    (model->GetScalarFunction("alpha_S"))->GetAs(PDF::isr::hard_subprocess);
  p_alpha     = dynamic_cast<MODEL::Running_AlphaQED *>
    (model->GetScalarFunction("alpha_QED"));
  // Initialize model parameters:
  // - pt_0, the IR regulator in the propagator and in the strong coupling
  // - pt_min, the IR cut-off for the 2->2 scatters
  // - Ecms, the cms energy of the hadron collision
  m_pt02        = sqr((*mipars)("pt_0"));
  m_ptmin2      = sqr((*mipars)("pt_min"));
  m_ecms        = rpa->gen.Ecms();
  m_S = m_S_lab = m_ecms*m_ecms;
  m_ptmax2      = m_S/4.;
  // TODO: will have to make this part of an external input scheme
  m_scale_scheme   = "MPI";
  m_kfactor_scheme = "MPI";
  // These parameters are for the hard part of the Sudakov form factor.
  // It is given by int_{pt^2}^{s/4} dpt^2 dsigma/dpt^2, and we tabulate the
  // integral in nbins between s/4 and pt_min^2, resulting in a stepsize of
  // pt^2_step.  In each of the bins we MC integrate over the rapidities of the
  // two outgoing particles with MC_points points.
  // 10000 points will yield errors of about 1%.
  m_pt2bins  = size_t((*mipars)("nPT_bins"));
  m_sbins    = size_t((*mipars)("nS_bins"));
  m_MCpoints = size_t((*mipars)("nMC_points"));
  // Mass scheme for the subsequent parton shower.
  m_massmode = 1;
  SetPSMasses();
  // Now initialize the 2->2 scatters and prepare the integral for the
  // "Sudakov form factor", Eq. (37) of Sjostrand-van Zijl
  InitializeAllProcesses();
  return PrepareSudakovFactor();
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
  p_pdf[0]->Calculate(x1,Min(Max(m_muFfac*scale,p_pdf[0]->Q2Min()),p_pdf[0]->Q2Max()));
  p_pdf[1]->Calculate(x2,Min(Max(m_muFfac*scale,p_pdf[1]->Q2Min()),p_pdf[1]->Q2Max()));
}

double MI_Processes::
operator()(const double & shat,const double & that,const double & uhat,
	   const double & x1,const double & x2) {
  // Return the total parton-level scattering cross section, summed over all
  // contributing processes.  This implicitly assumes that the PDFs have already
  // been set.
  double pt2 = that*uhat/shat;
  CalcPDFs(x1,x2,pt2);
  m_lastxs   = 0.;
  for (auto mig : m_groups) {
    mig->SetScale(pt2);
    m_lastxs += (*mig)(shat,that,uhat);
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
  // In this method we bin pt^2 in nbins, distributed logarithmically, and
  // accumulate the normalised integral in steps, the result is being stored in m_intbins
  // (the differential un-normalised integral in the pt^2 bin is stored in m_diffbins and
  // is used in the test routines only)
  // intbins[i] = Sum_{pt_i^2}^{pt_nmax^2} dSigma(pt_i^2)/dpt^2 * [pt_{i+1}^2-pt_i^2]
  // where pt_nmax^2 = s/4, the maximal pt^2.
  // N.B.: Note that we count the bins "down", i.e. bin 0 is at pt_nmax^2.
  // N.B.: Using left steps, thereby somewhat overestimating the integral.
  axis sbins    = (m_variable_s ?
		   axis(m_sbins, 4.*m_ptmin2, m_S, axis_mode::linear) :
		   axis(1, m_S , m_S, axis_mode::linear) );
  axis pt2bins  = axis(m_pt2bins,m_ptmin2,m_ptmax2,axis_mode::log);
  p_sudakov     = new Sudakov_Argument(this,sbins,pt2bins);
  return true;
}

double MI_Processes::dSigma(const double & pt2) {
  // Estimated differnetial cross setion dsigma/dpt^2 in 1/GeV^4 for a given transverse
  // momentum:
  // 1/(16 pi) int_{-ymax}^{+ymax} dy_1 dy_2  [  x_1 f(x_1, pt^2) x_2 f(x_2, pt^2)
  //                                             |M(shat,that,uhat)|^2 / shat^2    ]
  // Here we pulled the factor g^4 = alpha_S^2(pt^2)/16 pi^2 out of the matrix element
  // such that the prefactor is pi instead of 1/(16 pi).
  // We select the two rapidities of the two outgoing massless particles flat in the
  // full interval and hit-or-miss by making sure the x values are inside the allowed
  // range:  x_{1,2} = xT/2 * [exp(+/- y1) + exp(+/- y2)] with xT = (4pt^2/S)^0.5.
  if (pt2<m_ptmin2 || 4.*pt2>m_S) return 0.;
  double xt       = sqrt(4.*pt2/m_S);
  double ymax     = log(1./xt*(1.+sqrt(1.-xt*xt)));
  double yvolume  = sqr(2.*ymax);
  double res      = 0.;
  for (size_t i=0;i<m_MCpoints;i++) {
    double y1     = ymax*(-1.+2.*ran->Get());
    double y2     = ymax*(-1.+2.*ran->Get());
    double x1     = xt * (exp(y1)  + exp(y2))/2.;
    double x2     = xt * (exp(-y1) + exp(-y2))/2.;
    if (x1>m_xmin[0] && x1<1. && x2>m_xmin[1] && x2<1. &&
	xt*xt<x1*x2) {
      double cost = sqrt(1.-Min(1.,(xt*xt)/(x1*x2)));
      double shat = x1 * x2 * m_S;
      double that = -0.5 * shat * (1.-cost);
      double uhat = -0.5 * shat * (1.+cost);
      // this is just the sum over the matrix elements, grouped by parton content:
      // [x_1 f_i(x_1, mu^2) x_2 f(x_2, mu^2)]  [pi / shat^2] *
      // {alpha_S, alpha alpha_S, alpha}^2 |M_ij(shat,that,uhat)|^2
      // where the couplings alpha reflect the string/eletromagnetic coupling and the
      // |M_{ij}|^2 are given by the operators of the underlying XS_Base
      // realised e.g. in QCD_Processes.C, and including colour factors.
      res += (*this)(shat,that,uhat,x1,x2) * yvolume;
    }
  }
  return res/double(m_MCpoints);
}

void MI_Processes::UpdateS(const double & s) {
  // Update c.m. energy for variable centre-of-mass energies:
  // relevant for processes involving EPA photons etc..
  // Recalculate the non-diffractive and other cross sections
  m_S      = s;
  m_ecms   = sqrt(m_S);
  m_pt02   = mipars->CalculatePT02(m_S);
  (*p_xsecs)(m_S);
  // need to upate pt02 and ptmin2 for new s as well.
  for (list<MI_Process_Group *>::iterator mig = m_groups.begin();
       mig!=m_groups.end();mig++)  (*mig)->SetPT02(m_pt02);
}
