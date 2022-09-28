#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "EXTRA_XS/Main/Single_Process.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;


MI_Processes::MI_Processes(bool variable_s) : ME_Generator_Base("Amisic::Processes"),
			       m_ptmax2(1.e12), m_sigmaND(1.), m_integral(0.),
			       m_test(false), m_variable_s(variable_s) {}

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
    (model->GetScalarFunction("alpha_S"));
  p_alpha     = dynamic_cast<MODEL::Running_AlphaQED *>
    (model->GetScalarFunction("alpha_QED"));
  // Initialize model parameters:
  // - pt_0, the IR regulator in the propagator and in the strong coupling
  // - pt_min, the IR cut-off for the 2->2 scatters
  // - Ecms, the cms energy of the hadron collision
  m_pt0         = (*mipars)("pt_0");
  m_pt02        = m_pt0*m_pt0;
  m_ptmin       = (*mipars)("pt_min");
  m_ptmin2      = m_ptmin*m_ptmin;
  m_ecms        = rpa->gen.Ecms();
  m_S = m_S_lab = m_ecms*m_ecms;
  m_ptmax2      = sqr(m_ecms/2.);
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
  m_sbins    = int((*mipars)("nS_bins"));
  m_intbins.resize(m_nbins);
  m_diffbins.resize(m_nbins);
  m_pt2step  = log(m_S/(4.*m_ptmin2))/double(m_nbins);
  // Try the same integration etc. in xT = 2pT/E
  m_xTmin    = 2.*m_ptmin/m_ecms;   m_xTmax    = 1.;
  m_xTstep   = log(m_xTmax/m_xTmin)/double(m_nbins);
  // Mass scheme for the subsequent parton shower.
  m_massmode       = 1;
  SetPSMasses();
  // Now initialize the 2->2 scatters and prepare the integral for the
  // "Sudakov form factor", Eq. (37) of Sjostrand-van Zijl
  InitializeAllProcesses();
  return m_variable_s ? FillCaches() : PrepareSudakovFactor();
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
  p_pdf[0]->Calculate(x1,Max(m_muFfac*scale,p_pdf[0]->Q2Min()));
  p_pdf[1]->Calculate(x2,Max(m_muFfac*scale,p_pdf[1]->Q2Min()));
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
  //double xTlast = m_xTmin*exp(m_xTstep*m_nbins), sigmalast = 0.;
  //double dxT, xT, sigma;
  double pt2last = m_ptmin2*exp(m_pt2step*m_nbins);
  double sigma, pt2, dpt2, sigmalast;
  for (int bin=m_nbins-1;bin>=0;bin--) {
    pt2             = m_ptmin2*exp(m_pt2step*bin);
    dpt2            = pt2last-pt2;
    sigma           = dSigma(pt2);
    m_diffbins[bin] = sigma;
    m_intbins[bin]  = m_integral += sigma * dpt2/m_sigmaND;
    msg_Debugging()<<"   Sudakov(pt = "<<sqrt(pt2)<<") = "<<m_intbins[bin]<<" from "
              <<((sigmalast + sigma)/2./m_sigmaND)<<" * "<<dpt2<<".\n";
    pt2last        = pt2;
    sigmalast = sigma;
  }
  m_integral *= m_sigmaND;
  if (m_test) Test();
  return true;
}

void MI_Processes::Test() {
  double pt2last = m_ptmin2*exp(m_pt2step*m_nbins);
  msg_Out()<<METHOD<<" calculated integral for Sudakov form factor starting at pt = "
	   <<sqrt(pt2last)<<" in "<<m_nbins<<" steps,\n"
	   <<"   sigma = "<<m_integral<<" 1/Gev^2 = "<<(m_integral*rpa->Picobarn()/1.e9)<<" mb, "
	   <<" sigma/sigmaND = "<<m_integral/m_sigmaND<<".\n";
  double pt = 5.;
  while (pt<1000.) {
    msg_Out()<<" Log[Sud(pt = "<<pt<<")] = "<<SudakovArgument(sqr(pt))
	     <<"  --> Int_pt2^s dqt2 dsigma/dqt2 (pt = "<<pt<<") = "
	     <<(SudakovArgument(sqr(pt))*m_sigmaND*rpa->Picobarn())
	     <<" pb for sigmaND = "<<m_sigmaND<<" 1/GeV^2\n"
	     <<"   Test interpolation: "<<dSigma(sqr(pt))<<" vs "
	     <<SudakovDiffArgument(sqr(pt))<<" = "
	     <<(2.*(dSigma(sqr(pt))-SudakovDiffArgument(sqr(pt)))/
		(dSigma(sqr(pt))+SudakovDiffArgument(sqr(pt))) * 100)<<"%.\n";
    pt*=10.;
  }
  //exit(1);
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
  double res      = 0.;
  for (size_t i=0;i<m_MCpoints;i++) {
    double y1     = ymax*(-1.+2.*ran->Get());
    double y2     = ymax*(-1.+2.*ran->Get());
    double x1     = xt * (exp(y1)  + exp(y2))/2.;
    double x2     = xt * (exp(-y1) + exp(-y2))/2.;
    if (x1<1.e-6 || x1>1. || x2<1.e-6 || x2>1. || xt*xt>x1*x2) continue;
    double cost   = sqrt(1.-Min(1.,(xt*xt)/(x1*x2)));
    double shat   = x1 * x2 * m_S;
    double that   = -0.5 * shat * (1.-cost);
    double uhat   = -0.5 * shat * (1.+cost);
    double dsigma = 0.;
    if (x1>m_xmin[0] && x2>m_xmin[1]) {
      CalcPDFs(x1,x2,pt2);
      dsigma = (*this)(shat,that,uhat) * PSfac;
    }
    res  += dsigma;
  }
  double result = res/double(m_MCpoints);
  msg_Debugging()<<"dSigma(pt = "<<sqrt(pt2)<<")/dpt^2 = "
            <<(result*rpa->Picobarn())<<" pb GeV^-2\n";
  return result;
}

const double MI_Processes::SudakovArgument(const double & pt2) const {
  if (m_variable_s)
    return SudakovArgumentForVariableS(pt2);
  else
    return SudakovArgumentForConstantS(pt2);
}

double MI_Processes::SudakovArgumentForVariableS(const double &pt2) const {
  // Linear interpolation between the pre-calculated points for the Sudakov form factor
  if (pt2>m_ptmax2 || pt2<m_ptmin2 || m_S < 4*m_ptmin2) return 0.;
  int sbin      = log(m_S/m_S_lab) / log(m_sstep);
  int ptbin     = int(1./m_pt2step*log(pt2/m_ptmin2));

  double s1   = m_S_lab * std::pow(m_sstep, sbin), s2 = m_S_lab * std::pow(m_sstep, sbin+1);
  double pt21 = m_ptmin2*exp(m_pt2step*(ptbin)), pt22 = m_ptmin2*exp(m_pt2step*(ptbin+1));
  double val11 = m_cache_intbins[sbin][ptbin],   val12 = m_cache_intbins[sbin][ptbin+1];
  double val21 = m_cache_intbins[sbin+1][ptbin], val22 = m_cache_intbins[sbin+1][ptbin+1];
  double val = val11 * (s2 - m_S) * (pt22 - pt2) +
               val12 * (s2 - m_S) * (pt2 - pt21) +
               val21 * (m_S - s1) * (pt22 - pt2) +
               val22 * (m_S - s1) * (pt2 - pt21);
  val *= 1. / (s2 - s1) / (pt22 - pt21);
  return val;
}

double MI_Processes::SudakovArgumentForConstantS(const double &pt2) const {
  // Linear interpolation between the pre-calculated points for the Sudakov form factor
  if (pt2>m_ptmax2 || pt2<m_ptmin2) return 0.;
  int bin     = int(1./m_pt2step*log(pt2/m_ptmin2));
  double pt21 = m_ptmin2*exp(m_pt2step*(bin)), pt22 = m_ptmin2*exp(m_pt2step*(bin+1));
  double val1 = m_intbins[bin],                val2 = m_intbins[bin+1];
  double val  = (val1*(pt22-pt2)+val2*(pt2-pt21))/(pt22-pt21);
  return val;
}

const double MI_Processes::SudakovDiffArgument(const double & pt2) const {
  // Linear interpolation between the pre-calculated points for the Sudakov form factor
  if (pt2>m_ptmax2 || pt2<m_ptmin2) return 0.;
  int bin     = int(1./m_pt2step*log(pt2/m_ptmin2));
  double pt21 = m_ptmin2*exp(m_pt2step*(bin)), pt22 = m_ptmin2*exp(m_pt2step*(bin+1));
  double val1 = m_diffbins[bin],               val2 = m_diffbins[bin+1];
  double val  = (val1*(pt22-pt2)+val2*(pt2-pt21))/(pt22-pt21);
  return val;
}

void MI_Processes::Update(double s) {
  // taken from the initialization
  m_S       = s;
  m_ecms    = sqrt(s);
  m_ptmax2  = sqr(m_ecms/2.);
  m_pt2step = log(m_S/(4.*m_ptmin2))/double(m_nbins);
  m_xTmin   = 2.*m_ptmin/m_ecms;
  CalculateIntegralFromCache();
}

bool MI_Processes::FillCaches() {
  msg_Out() << METHOD << ": Filling cache for multi-parton interactions, for " << m_sbins << " bins: \n";
  m_test = false;
  m_sstep = std::pow(4*m_ptmin2/m_S_lab, 1./m_sbins);
  m_cache_diffbins.resize(m_sbins,std::vector<double>(m_nbins));
  m_cache_intbins.resize(m_sbins,std::vector<double>(m_nbins));
  m_cache_integral.resize(m_sbins);
  for (int sbin = 0; sbin < m_sbins; ++sbin) {
    msg_Info() << "  Integrating bin " << sbin+1 << " of " << m_sbins << ". \n";
    if (!(rpa->gen.BatchMode()&2) && sbin != m_sbins-1) msg_Info() << mm(1,mm::up);
    m_S = m_S_lab * std::pow(m_sstep, sbin);
    (*p_xsecs)(m_S);
    m_sigmaND = p_xsecs->XSnd();
    double pt2last = m_ptmin2*exp(m_pt2step*m_nbins);
    double sigma, pt2, dpt2;
    for (int ptbin=m_nbins-1;ptbin>=0;ptbin--) {
      pt2             = m_ptmin2*exp(m_pt2step*ptbin);
      dpt2            = pt2last-pt2;
      sigma           = dSigma(pt2);
      m_cache_diffbins[sbin][ptbin] = sigma;
      m_cache_intbins[sbin][ptbin] = sigma * dpt2/m_sigmaND;
      m_cache_integral[sbin] += sigma * dpt2;
      pt2last        = pt2;
    }
  }
  msg_Info() << "  Caching successfully completed. \n";
  m_integral = m_cache_integral[m_sbins-1];
  return true;
}

void MI_Processes::CalculateIntegralFromCache() {
  int sbin = log(m_S/m_S_lab) / log(m_sstep);
  double s1 = m_S_lab * std::pow(m_sstep, sbin);
  double s2 = m_S_lab * std::pow(m_sstep, sbin+1);
  double val1 = m_cache_integral[sbin], val2 = m_cache_integral[sbin+1];
  m_integral = (val1 * (s2 - m_S) + val2 * (m_S - s1)) / (s2 - s1);
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

