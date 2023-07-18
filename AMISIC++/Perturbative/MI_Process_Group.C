#include "AMISIC++/Perturbative/MI_Process_Group.H"
#include "AMISIC++/Perturbative/QCD_Processes.H"
#include "AMISIC++/Perturbative/QED_Processes.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Base class
// A container for sufficiently similar processes to be clustered for joint
// integration and generation of scatters
///////////////////////////////////////////////////////////////////////////////

MI_Process_Group::MI_Process_Group(const std::string & name) :
  m_name(name), m_lastxs(0.), m_pref(M_PI)
{
  m_muR_scheme = mipars->GetScaleScheme();
  m_muR_fac    = sqr((*mipars)("RenScale_Factor"));
  m_muF_fac    = sqr((*mipars)("FacScale_Factor"));
  m_pt02       = sqr((*mipars)("pt_0"));
}

MI_Process_Group::~MI_Process_Group() {
  while (!m_processes.empty()) {
    delete m_processes.back();
    m_processes.pop_back();
  }
  m_me2s.unique();
  while (!m_me2s.empty()) {
    delete m_me2s.back();
    m_me2s.pop_back();
  }
  m_processes.clear();
  m_me2s.clear();
}

double MI_Process_Group::
operator()(const double & shat,const double & that,const double & uhat) {
  // This calculates dsigma/dpt^2 in units of 1/GeV^4
  // It is given by
  // pi/shat^2 * f_1(x_1) f_2(x_2) * alpha_{S,QED}^2 |ME|^2 * [pt^2/(pt^2+pt0^2)]^2
  // where the last term regularises the t/u-channel singluarities
  PreCalculate(shat,that,uhat);
  double tot  = 0.;
  for (list<MI_Process * >::iterator mit=m_processes.begin();
       mit!=m_processes.end();mit++) {
    tot += ( p_pdf[0]->GetXPDF((*mit)->Flav(0)) *
		  p_pdf[1]->GetXPDF((*mit)->Flav(1)) ) * (**mit)();
    if (std::isnan(tot))
      tot = 0.;
  }
  m_lastxs = m_pref/sqr(shat) *
             Coupling(Scale(m_scale)) *
             SoftCorrection(m_scale) * tot;
  return m_lastxs;
}

void MI_Process_Group::
PreCalculate(const double & shat,const double & that,const double & uhat) {
  // Calculating the squared MEs depending on the Mandelstam variables.
  // They do not include couplings, flux factors, pdf's, and are in units of 1.
  for (list<XS_Base * >::iterator xsit=m_me2s.begin();
       xsit!=m_me2s.end();xsit++) 
    (*xsit)->Calc(shat,that,uhat);
}

double MI_Process_Group::SoftCorrection(const double & pt2) const {
  // Getting rid of the t/u-channel singularities
  return sqr(pt2/(pt2+m_pt02));
}

double MI_Process_Group::Scale(const double & pt2) const {
  // Default scale, including an IR regularisation - maybe we should get more choices.
  switch (m_muR_scheme) {
  case scale_scheme::code::PT_with_Raps:
    THROW(fatal_error,"Scale scheme PT_with_Rapidities not implemented yet.");
  case scale_scheme::code::PT:
  default:
    return (pt2+m_pt02);
  }
}

void MI_Process_Group::Output() const {
  msg_Out()<<"Group: "<<Name()<<":\n";
  for (list<MI_Process * >::const_iterator mit=m_processes.begin();
       mit!=m_processes.end();mit++) 
    msg_Out()<<"   *** "<<(*mit)->Name()<<"\n";
}

MI_Process * MI_Process_Group::SelectProcess() {
  // Selects a process according to the relative differential cross sections
  // at the given kinematic configuration.
  double tot = 0.;
  for (list<MI_Process *>::iterator mipit=m_processes.begin();
       mipit!=m_processes.end();mipit++) {
    MI_Process * mip = (*mipit);
    tot += (p_pdf[0]->GetXPDF(mip->Flav(0)) *
	    p_pdf[1]->GetXPDF(mip->Flav(1)) * (*mip)());
  }
  tot *= ran->Get();
  for (auto mipit : m_processes) {
    tot -= (p_pdf[0]->GetXPDF(mipit->Flav(0)) *
	    p_pdf[1]->GetXPDF(mipit->Flav(1)) * (*mipit)());
    if (tot<=0.) return mipit;
  }
  return m_processes.back();
}

///////////////////////////////////////////////////////////////////////////////
// GG initiated processes
///////////////////////////////////////////////////////////////////////////////

MI_GG_Processes::MI_GG_Processes() :
  MI_Process_Group("MPI_gg_processes"), m_Nqq(0) {
  XS_Base * gg2gg(new gg_gg()), * gg2qqbar(new gg_qqbar());
  m_me2s.push_back(gg2gg);
  m_me2s.push_back(gg2qqbar);

  vector<Flavour> flavs;
  Flavour gluon(kf_gluon);
  for (size_t i=0;i<4;i++) flavs.push_back(gluon);
  m_processes.push_back(new MI_Process(flavs));
  m_processes.back()->SetME2(gg2gg);
  // add five identical copies for gg -> qq.
  for (size_t i=1;i<6;i++) {
    if (Flavour(i).Mass()>0.) continue;
    flavs[2] = Flavour(i); flavs[3] = flavs[2].Bar();
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(gg2qqbar);
    m_Nqq++;
  }
}

double MI_GG_Processes::Coupling(const double & scale) const {
  return sqr((*p_alphaS)(Max(m_pt02,m_muR_fac*Scale(m_scale))));
}

MI_Process * MI_GG_Processes::SelectProcess() {
  double tot  = 0.;
  list<MI_Process * >::iterator mit;
  for (mit=m_processes.begin();mit!=m_processes.end();mit++) tot += (**mit)();
  tot *= ran->Get();
  mit = m_processes.begin();
  while (tot>0.) {
    tot -= (**mit)();
    mit++;
    if (mit==m_processes.end()) {
      mit--;
      break;
    }
  }
  return (*mit);
}



///////////////////////////////////////////////////////////////////////////////
// qqbar and qbarq initiated processes
///////////////////////////////////////////////////////////////////////////////

MI_QQB_Processes::MI_QQB_Processes():
  MI_Process_Group("MPI_qqb_processes") {
  XS_Base * qqbar2qqbar(new qqbar_qqbar()), * qqbar2gg(new qqbar_gg()),
    *q1q1bar2q2q2bar(new q1q1bar_q2q2bar());
  m_me2s.push_back(qqbar2qqbar);
  m_me2s.push_back(qqbar2gg);
  m_me2s.push_back(q1q1bar2q2q2bar);

  vector<Flavour> flavs;
  flavs.resize(4);
  for (size_t i=2;i<3;i++) {
    if (Flavour(i).Mass()>0.) continue;
    // q qbar -> q qbar
    flavs[0] = flavs[2] = Flavour(i);
    flavs[1] = flavs[3] = flavs[0].Bar();
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qqbar2qqbar);
    // qbar q -> qbar q
    for (size_t j=0;j<4;j++) flavs[j] = flavs[j].Bar();
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qqbar2qqbar);
    // q qbar -> gg
    flavs[2] = flavs[3] = Flavour(kf_gluon);
    flavs[0] = Flavour(i);
    flavs[1] = flavs[0].Bar();
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qqbar2gg);
    // qbar q -> g g
    for (size_t j=0;j<2;j++) flavs[j] = flavs[j].Bar();
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qqbar2gg);
    // now loop over all flavours, where i and j are different
    for (size_t j=1;j<2;j++) {
      if (i==j || Flavour(j).Mass()>0.) continue;
      // q1 q1bar -> q2 q2bar
      flavs[0] = Flavour(i);
      flavs[1] = flavs[0].Bar();
      flavs[2] = Flavour(j);
      flavs[3] = flavs[2].Bar();
      m_processes.push_back(new MI_Process(flavs));
      m_processes.back()->SetME2(q1q1bar2q2q2bar);
      // q1bar q1 -> q2bar q2
      for (size_t k=0;k<4;k++) flavs[k] = flavs[k].Bar();
      m_processes.push_back(new MI_Process(flavs));
      m_processes.back()->SetME2(q1q1bar2q2q2bar);
    }
  }
}

double MI_QQB_Processes::Coupling(const double & scale) const {
  return sqr((*p_alphaS)(Max(m_pt02,m_muR_fac*Scale(m_scale))));
}

///////////////////////////////////////////////////////////////////////////////
// qq and qbar qbar initiated processes
///////////////////////////////////////////////////////////////////////////////

MI_QQ_Processes::MI_QQ_Processes():
  MI_Process_Group("MPI_qq_processes") {
  XS_Base * qq2qq(new qq_qq());
  m_me2s.push_back(qq2qq);

  vector<Flavour> flavs;
  flavs.resize(4);
  for (size_t i=1;i<6;i++) {
    if (Flavour(i).Mass()>0.) continue;
    flavs[0] = flavs[1] = flavs[2] = flavs[3] = Flavour(i);
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qq2qq);
    flavs[0] = flavs[1] = flavs[2] = flavs[3] = Flavour(i).Bar();
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qq2qq);
  }
}

double MI_QQ_Processes::Coupling(const double & scale) const {
  return sqr((*p_alphaS)(Max(m_pt02,m_muR_fac*Scale(m_scale))));
}

MI_Process * MI_QQ_Processes::SelectProcess() {
  double pdf = 0.;
  for (list<MI_Process *>::iterator mipit=m_processes.begin();
       mipit!=m_processes.end();mipit++) {
    MI_Process * mip = (*mipit);
    pdf += p_pdf[0]->GetXPDF(mip->Flav(0)) * p_pdf[1]->GetXPDF(mip->Flav(1));
  }
  pdf *= ran->Get();
  for (auto mipit : m_processes) {
    pdf -= p_pdf[0]->GetXPDF(mipit->Flav(0)) *
           p_pdf[1]->GetXPDF(mipit->Flav(1));
    if (pdf<=0.) return mipit;
  }
  return m_processes.back();
}

///////////////////////////////////////////////////////////////////////////////
// qg and gq initiated processes
///////////////////////////////////////////////////////////////////////////////


MI_QG_Processes::MI_QG_Processes():
  MI_Process_Group("MPI_qg_processes") {
  XS_Base * qg2qg(new qg_qg());
  m_me2s.push_back(qg2qg);

  vector<Flavour> flavs;
  flavs.resize(4);
  for (size_t i=1;i<6;i++) {
    if (Flavour(i).Mass()>0.) continue;
    flavs[0] = flavs[2] = Flavour(kf_gluon);
    flavs[1] = flavs[3] = Flavour(i);
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qg2qg);
    flavs[1] = flavs[3] = Flavour(i).Bar();
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qg2qg);

    flavs[1] = flavs[3] = Flavour(kf_gluon);
    flavs[0] = flavs[2] = Flavour(i);
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qg2qg);
    flavs[0] = flavs[2] = Flavour(i).Bar();
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qg2qg);
  }
}

double MI_QG_Processes::Coupling(const double & scale) const {
  return sqr((*p_alphaS)(Max(m_pt02,m_muR_fac*Scale(m_scale))));
}

MI_Process * MI_QG_Processes::SelectProcess() {
  double pdf = 0.;
  for (list<MI_Process *>::iterator mipit=m_processes.begin();
       mipit!=m_processes.end();mipit++) {
    MI_Process * mip = (*mipit);
    pdf += p_pdf[0]->GetXPDF(mip->Flav(0))*p_pdf[1]->GetXPDF(mip->Flav(1));
  }
  pdf *= ran->Get();
  for (auto mipit : m_processes) {
    pdf -= p_pdf[0]->GetXPDF(mipit->Flav(0)) *
           p_pdf[1]->GetXPDF(mipit->Flav(1));
    if (pdf<=0.) return mipit;
  }
  return m_processes.back();
}

///////////////////////////////////////////////////////////////////////////////
// qq', qqbar', qbar q', and qbar qbar' initiated processes
///////////////////////////////////////////////////////////////////////////////

MI_Q1Q2_Processes::MI_Q1Q2_Processes():
  MI_Process_Group("MPI_q1q2_processes") {
  XS_Base * q1q22q1q2(new q1q2_q1q2());
  m_me2s.push_back(q1q22q1q2);

  vector<Flavour> flavs;
  flavs.resize(4);
  for (size_t i=1;i<6;i++) {
    if (Flavour(i).Mass()>0.) continue;
    for (size_t j=1;j<6;j++) {
      if (i==j || Flavour(j).Mass()>0.) continue;
      flavs[0] = flavs[2] = Flavour(i);
      flavs[1] = flavs[3] = Flavour(j);
      m_processes.push_back(new MI_Process(flavs));
      m_processes.back()->SetME2(q1q22q1q2);
      flavs[0] = flavs[2] = Flavour(i).Bar();
      flavs[1] = flavs[3] = Flavour(j);
      m_processes.push_back(new MI_Process(flavs));
      m_processes.back()->SetME2(q1q22q1q2);
      flavs[0] = flavs[2] = Flavour(i);
      flavs[1] = flavs[3] = Flavour(j).Bar();
      m_processes.push_back(new MI_Process(flavs));
      m_processes.back()->SetME2(q1q22q1q2);
      flavs[0] = flavs[2] = Flavour(i).Bar();
      flavs[1] = flavs[3] = Flavour(j).Bar();
      m_processes.push_back(new MI_Process(flavs));
      m_processes.back()->SetME2(q1q22q1q2);
    }
  }
}

double MI_Q1Q2_Processes::Coupling(const double & scale) const {
  return sqr((*p_alphaS)(Max(m_pt02,m_muR_fac*Scale(m_scale))));
}

MI_Process * MI_Q1Q2_Processes::SelectProcess() {
  double pdf = 0.;
  for (list<MI_Process *>::iterator mipit=m_processes.begin();
       mipit!=m_processes.end();mipit++) {
    MI_Process * mip = (*mipit);
    pdf += p_pdf[0]->GetXPDF(mip->Flav(0))*p_pdf[1]->GetXPDF(mip->Flav(1));
  }
  pdf *= ran->Get();
  for (auto mipit : m_processes) {
    pdf -= p_pdf[0]->GetXPDF(mipit->Flav(0)) *
           p_pdf[1]->GetXPDF(mipit->Flav(1));
    if (pdf<=0.) return mipit;
  }
  return m_processes.back();
}

///////////////////////////////////////////////////////////////////////////////
// qg initiated photon production
///////////////////////////////////////////////////////////////////////////////

MI_QG_QGamma_Processes::MI_QG_QGamma_Processes() :
  MI_Process_Group("MPI_qg_qgamma_processes") {
  XS_Base * qg_qgam(new qg_qgamma());
  m_me2s.push_back(qg_qgam);

  vector<Flavour> flavs;
  flavs.resize(4);
  Flavour gluon(kf_gluon), photon(kf_photon);
  for (size_t i=1;i<6;i++) {
    if (Flavour(i).Mass()>0.) continue;
    flavs[0] = flavs[2] = Flavour(i);
    flavs[1] = gluon;
    flavs[3] = photon;
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qg_qgam);
    flavs[0] = flavs[2] = Flavour(i).Bar();
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qg_qgam);
    flavs[1] = flavs[3] = Flavour(i);
    flavs[0] = gluon;
    flavs[2] = photon;
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qg_qgam);
    flavs[1] = flavs[3] = Flavour(i).Bar();
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qg_qgam);
  }
}
  
double MI_QG_QGamma_Processes::Coupling(const double & scale) const {
  return ( (*p_alphaS)(m_muR_fac*Scale(Max(m_pt02,m_scale))) *
	   (*p_alpha)(Max(m_pt02,Scale(m_scale)) ));
}

MI_Process * MI_QG_QGamma_Processes::SelectProcess() {
  double pdf = 0.;
  for (list<MI_Process *>::iterator mipit=m_processes.begin();
       mipit!=m_processes.end();mipit++) {
    MI_Process * mip = (*mipit);
    pdf += p_pdf[0]->GetXPDF(mip->Flav(0))*p_pdf[1]->GetXPDF(mip->Flav(1)) *
      sqr(mip->Flav(0).IsQuark()?mip->Flav(0).Charge():mip->Flav(1).Charge());
  }
  pdf *= ran->Get();
  for (auto mipit : m_processes) {
    pdf -= p_pdf[0]->GetXPDF(mipit->Flav(0)) *
           p_pdf[1]->GetXPDF(mipit->Flav(1)) *
           sqr(mipit->Flav(0).IsQuark()?mipit->Flav(0).Charge()
                                           :mipit->Flav(1).Charge());
    if (pdf<=0.) return mipit;
  }
  return m_processes.back();
}

///////////////////////////////////////////////////////////////////////////////
// qqbar initated photon proudction
///////////////////////////////////////////////////////////////////////////////

MI_QQ_GGamma_Processes::MI_QQ_GGamma_Processes() :
  MI_Process_Group("MPI_qqbar_ggamma_processes") {
  XS_Base * qq_ggam(new qqbar_ggamma());
  m_me2s.push_back(qq_ggam);

  vector<Flavour> flavs;
  flavs.resize(4);
  flavs[2] = Flavour(kf_gluon);
  flavs[3] = Flavour(kf_photon);
  for (size_t i=1;i<6;i++) {
    if (Flavour(i).Mass()>0.) continue;
    flavs[0] = Flavour(i);
    flavs[1] = flavs[0].Bar();
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qq_ggam);
    flavs[0] = Flavour(i).Bar();
    flavs[1] = flavs[0].Bar();
    m_processes.push_back(new MI_Process(flavs));
    m_processes.back()->SetME2(qq_ggam);
  }
}

double MI_QQ_GGamma_Processes::Coupling(const double & scale) const {
  return sqr((*p_alpha)(Scale(m_scale)));
}

MI_Process * MI_QQ_GGamma_Processes::SelectProcess() {
  double pdf = 0.;
  for (list<MI_Process *>::iterator mipit=m_processes.begin();
       mipit!=m_processes.end();mipit++) {
    MI_Process * mip = (*mipit);
    pdf += p_pdf[0]->GetXPDF(mip->Flav(0))*p_pdf[1]->GetXPDF(mip->Flav(1)) *
      sqr(mip->Flav(0).Charge());
  }
  pdf *= ran->Get();
  for (auto mipit : m_processes) {
    pdf -= p_pdf[0]->GetXPDF(mipit->Flav(0)) *
           p_pdf[1]->GetXPDF(mipit->Flav(1)) *
           sqr(mipit->Flav(0).Charge());
    if (pdf<=0.) return mipit;
  }
  return m_processes.back();
}
