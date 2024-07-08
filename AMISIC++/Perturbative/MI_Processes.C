#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "EXTRA_XS/Main/Single_Process.H"
#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;


MI_Processes::MI_Processes() :
  ME_Generator_Base("Amisic"),
  m_integrator(MI_Integrator(this)),
  p_xsecs(NULL), m_fm2GeV2(sqr(rpa->hBar()*rpa->c()*1.e12)) {}

MI_Processes::~MI_Processes() {
  while (!m_groups.empty()) {
    delete m_groups.back();
    m_groups.pop_back();
  }
  m_groups.clear();
}

bool MI_Processes::Initialize(MODEL::Model_Base* const          model,
                              BEAM::Beam_Spectra_Handler* const beam,
                              PDF::ISR_Handler* const           isr,
                              YFS::YFS_Handler* const           yfs)
{
  ///////////////////////////////////////////////////////////////////////////
  // Get PDFs and couplings
  m_muF_scheme  = mipars->GetScaleFScheme();
  m_muF_fac     = sqr((*mipars)("FacScale_Factor"));
  for (size_t i=0;i<2;i++) {
    p_pdf[i]      = isr->PDF(i);
    m_xmin[i]     = p_pdf[i]->XMin();
    m_xmax[i]     = p_pdf[i]->XMax();
    p_remnants[i] = isr->GetRemnant(i);
    m_resx[i]     = 1.;
  }
  m_muR_scheme  = mipars->GetScaleRScheme();
  m_muR_fac     = sqr((*mipars)("RenScale_Factor"));
  p_alphaS      = dynamic_cast<MODEL::Running_AlphaS *>
    (model->GetScalarFunction("alpha_S"))->GetAs(PDF::isr::hard_subprocess);
  p_alpha       = dynamic_cast<MODEL::Running_AlphaQED *>
    (model->GetScalarFunction("alpha_QED"));
  ///////////////////////////////////////////////////////////////////////////
  // Initialize model parameters:
  // - pt_0, the IR regulator in the propagator and in the strong coupling
  // - pt_min, the IR cut-off for the 2->2 scatters
  // - Ecms, the cms energy of the hadron collision
  ///////////////////////////////////////////////////////////////////////////
  m_pt02        = sqr((*mipars)("pt_0"));
  m_ptmin2      = sqr((*mipars)("pt_min"));
  m_ecms        = rpa->gen.Ecms();
  m_Emin        = (*mipars)("E_min");
  m_S = m_S_lab = m_ecms*m_ecms;
  m_ptmax2      = m_S/4.;
  ///////////////////////////////////////////////////////////////////////////
  // Now initialize the 2->2 scatters and the integrator ...
  ///////////////////////////////////////////////////////////////////////////
  InitializeAllProcesses();
  m_integrator.Initialize(isr);
  ///////////////////////////////////////////////////////////////////////////
  // ... and integrate away.
  ///////////////////////////////////////////////////////////////////////////
  (*p_xsecs)(m_S);
  msg_Info()<<"   "<<std::string(77,'-')<<"\n\n";
  ///////////////////////////////////////////////////////////////////////////
  // Mass scheme for the subsequent parton shower.
  ///////////////////////////////////////////////////////////////////////////
  m_massmode = 1;
  SetPSMasses();
  return true;
}

bool MI_Processes::InitializeAllProcesses() {
  ///////////////////////////////////////////////////////////////////////////
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
  ///////////////////////////////////////////////////////////////////////////
  m_groups.push_back(new MI_GG_Processes());
  m_groups.push_back(new MI_Q1Q2_Processes());
  m_groups.push_back(new MI_QQB_Processes());
  m_groups.push_back(new MI_QQ_Processes());
  m_groups.push_back(new MI_QG_Processes());
  ///////////////////////////////////////////////////////////////////////////
  // The following processes should depend on switches.  At the moment we
  // just add them without further ado.
  // - qg-initiated photon production
  // - qqbar-initiated photon production
  // We are missing di-photon production here:
  // - gg->gamma gamma and qqbar->gamma gamma
  ///////////////////////////////////////////////////////////////////////////
  m_groups.push_back(new MI_QG_QGamma_Processes());
  m_groups.push_back(new MI_QQ_GGamma_Processes());
  ///////////////////////////////////////////////////////////////////////////
  // We are missing the production of (heavy quarkonia) mesons MQQ in
  // - singlet production qqbar -> MQQ, gg -> MQQ, gq -> MQQ+q etc.
  // - octet production gg -> MQQ^(8)+g, qqbar->MQQ^(8)+g etc.
  // We could also add production of gauge bosons:
  // - qqbar->Z, qqbar'->W
  // - gq->Zq, gq->Wq', etc.
  ///////////////////////////////////////////////////////////////////////////
  SetPDFs();
  SetAlphaS();
  return true;
}

double MI_Processes::operator()() {
  ///////////////////////////////////////////////////////////////////////////
  // Assuming the mandelstam parameters etc. have already been determined
  // in the integrator, we can just call extract them.
  ///////////////////////////////////////////////////////////////////////////
  double shat = m_integrator.SHat();
  double that = m_integrator.THat();
  double uhat = m_integrator.UHat();
  double x1   = m_integrator.X(0);
  double x2   = m_integrator.X(1);
  return (*this)(shat,that,uhat,x1,x2);
}

double MI_Processes::
operator()(const double & shat,const double & that,const double & uhat,
	   const double & x1,const double & x2) {
  ///////////////////////////////////////////////////////////////////////////
  // Return the total parton-level scattering cross section, summed over all
  // contributing processes.  This implicitly assumes that the PDFs have
  // already been set.
  ///////////////////////////////////////////////////////////////////////////
  m_lastxs   = 0.;
  if (x1<=m_xmin[0]/m_resx[0] || x1>=m_xmax[0] ||
      x2<=m_xmin[1]/m_resx[1] || x2>=m_xmax[1]) return 0.;
  CalcScales(shat,that,uhat);
  CalcPDFs(x1,x2);
  /////////////////////////////////////////////////////////////////////////
  // This is the sum over the matrix elements, grouped by parton content:
  // [x_1 f_i(x_1, mu^2) x_2 f(x_2, mu^2)]  [pi / shat^2] *
  //         {alpha_S, alpha alpha_S, alpha}^2 |M_ij(shat,that,uhat)|^2
  // where the couplings alpha reflect the string/eletromagnetic coupling
  // and the |M_{ij}|^2 are given by the operators of the underlying
  // XS_Base realised e.g. in QCD_Processes.C, and include colour factors.
  /////////////////////////////////////////////////////////////////////////
  for (auto mig : m_groups) {
    mig->SetScale(m_muR2);
    m_lastxs += (*mig)(shat,that,uhat);
  }
  return m_lastxs;
}

int MI_Processes::FillHardScatterBlob(Blob *&  blob,const double & pt2veto) {
  ///////////////////////////////////////////////////////////////////////////
  // Selecting a process for the hard collision and, if allowed, returning
  // a new blob for the perturbative scatter, filled with scales and incoming
  // and outgoing particles.
  ///////////////////////////////////////////////////////////////////////////
  MI_Process * proc = SelectProcess();
  if (proc==nullptr) msg_Error()<<"Warning in "<<METHOD<<": no proc selected.\n";
  if (proc==nullptr ||
      !proc->MakeKinematics(&m_integrator,p_remnants) ||
      !proc->SetColours()) return 0;
  if (pt2veto>0. && m_integrator.PT2()>pt2veto) return -1;
  array<int,2> inflavs;
  for (size_t i=0;i<2;i++) {
    Particle * part = proc->GetParticle(i);
    blob->AddToInParticles(part);
    inflavs[i] = (part->Flav().IsAnti() ? -1 : 1) * part->Flav().Kfcode();
  }
  for (size_t i=2;i<4;i++) blob->AddToOutParticles(proc->GetParticle(i));
  blob->AddData("WeightsMap",new Blob_Data<Weights_Map>({}));
  blob->AddData("Renormalization_Scale",new Blob_Data<double>(m_muR2));
  blob->AddData("Factorization_Scale",new Blob_Data<double>(m_muF2));
  blob->AddData("Resummation_Scale",new Blob_Data<double>(Max(m_muR2,m_muF2)));
  PDF_Info info(inflavs[0],inflavs[1],
		m_integrator.X(0),m_integrator.X(1),
		m_muF2,m_muF2);
  blob->AddData("PDFInfo",new Blob_Data<PDF_Info>(info));
  return 1;
}

double MI_Processes::TotalCrossSection(const double & s,const bool & output) {
  ///////////////////////////////////////////////////////////////////////////
  // Calculate the hard cross section first, by iterating over the pt2 bins
  ///////////////////////////////////////////////////////////////////////////
  m_xshard      = m_integrator(s,nullptr,0.);
  if (output) {
    msg_Info()<<"   "<<std::string(77,'-')<<"\n"
	      <<"   | "<<METHOD<<": xs_pert = "
	      <<std::setprecision(4)<<std::setw(10)
	      <<(m_xshard*rpa->Picobarn()/1.e9)<<" mb "
	      <<"+- "<<std::setprecision(0)<<std::setw(3)
	      <<(100.*m_integrator.Uncertainty()/m_xshard)
	      <<"%."<<std::string(9,' ')<<"|\n";
  }
  return m_xshard;
}

void MI_Processes::
CalcScales(const double & shat,const double & that,const double & uhat) {
  double pt2 = that*uhat/shat;
  switch (m_muR_scheme) {
  case scale_scheme::PT_with_Raps:
    THROW(fatal_error,"Scale scheme PT_with_Raps not implemented yet!")
  case scale_scheme::PT:
  default:
    m_muR2 = m_muR_fac*(pt2 + m_pt02);
  }

  switch (m_muF_scheme) {
  case scale_scheme::PT_with_Raps:
    THROW(fatal_error,"Scale scheme PT_with_Raps not implemented yet!")
  case scale_scheme::PT:
  default:
    m_muF2 = m_muF_fac*(pt2 + m_pt02);
  }
}

void MI_Processes::CalcPDFs(const double & x1,const double & x2) {
  ///////////////////////////////////////////////////////////////////////////
  // Calculate both sets of PDFs at the relevant x and Q^2
  ///////////////////////////////////////////////////////////////////////////
  p_pdf[0]->Calculate(x1,m_muF2);
  p_pdf[1]->Calculate(x2,m_muF2);
}

void MI_Processes::SetAlphaS() {
  ///////////////////////////////////////////////////////////////////////////
  // The couplings for the process groups.
  ///////////////////////////////////////////////////////////////////////////
  for (list<MI_Process_Group *>::iterator mig = m_groups.begin();
       mig!=m_groups.end();mig++) {
    (*mig)->SetAlphaS(p_alphaS);
    (*mig)->SetAlpha(p_alpha);
  }
}

void MI_Processes::SetPDFs() {
  ///////////////////////////////////////////////////////////////////////////
  // The PDFs for the process groups.
  ///////////////////////////////////////////////////////////////////////////
  for (list<MI_Process_Group *>::iterator mig = m_groups.begin();
       mig!=m_groups.end();mig++)
    (*mig)->SetPDFs(p_pdf[0],p_pdf[1]);
}

MI_Process * MI_Processes::SelectProcess() {
  ///////////////////////////////////////////////////////////////////////////
  // Sum over all cross sections of all groups and select one of the groups.
  // Then select one of the processes within the group.
  ///////////////////////////////////////////////////////////////////////////
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

void MI_Processes::UpdateS(const double & s) {
  ///////////////////////////////////////////////////////////////////////////
  // Update c.m. energy for variable centre-of-mass energies:
  // relevant for processes involving EPA photons etc..
  // Recalculate the non-diffractive and other cross sections
  ///////////////////////////////////////////////////////////////////////////
  m_S      = s;
  m_ecms   = sqrt(m_S);
  m_pt02   = mipars->CalculatePT02(m_S);
  (*p_xsecs)(m_S);
  ///////////////////////////////////////////////////////////////////////////
  // need to upate pt02 and ptmin2 for new s as well.
  ///////////////////////////////////////////////////////////////////////////
  for (list<MI_Process_Group *>::iterator mig = m_groups.begin();
       mig!=m_groups.end();mig++)  (*mig)->SetPT02(m_pt02);
}
