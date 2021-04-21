#include "SHRiMPS/Event_Generation/Soft_Diffractive_Event_Generator.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace SHRIMPS;

Soft_Diffractive_Event_Generator::
Soft_Diffractive_Event_Generator(Sigma_SD * sigma,const int & test) :
  Event_Generator_Base(sigma),
  p_sigma(sigma), 
  m_sigma(0.)
{
  for (size_t i=0;i<3;i++) m_sigma   += m_rate[i] = p_sigma->GetXSec(i);  
  for (size_t i=0;i<3;i++) m_rate[i] /= m_sigma;
  for (size_t i=0;i<2;i++) {
    m_p[i] = ATOOLS::rpa->gen.PBeam(i);
    m_E[i] = m_p[i][0];
  }
  m_sign1 = -1+2*int(m_p[0][3]>0);
  InitialiseHadronMaps();
  m_histomap[std::string("Q_SD")] = new Histogram(0,0.0,1.0,1000);
  m_xsec = m_sigma;
}

void Soft_Diffractive_Event_Generator::InitialiseHadronMaps() {
  // Assume pp/ppbar collisions only
  // this can easily be extended to mesons/neutrons if necessary
  for (size_t i=0;i<2;i++) {
    m_beam[i] = (i==0)?ATOOLS::rpa->gen.Beam1():ATOOLS::rpa->gen.Beam2();
    if (m_beam[i]==ATOOLS::Flavour(kf_p_plus)) {
      m_hadronmaps[i][ATOOLS::Flavour(kf_N_1440_plus)] = 1.;
    }
    if (m_beam[i]==ATOOLS::Flavour(kf_p_plus).Bar()) {
      m_hadronmaps[i][ATOOLS::Flavour(kf_N_1440_plus).Bar()] = 1.;
    }
  }
}

Soft_Diffractive_Event_Generator::~Soft_Diffractive_Event_Generator() {
  if (!m_histomap.empty()) {
    Histogram * histo;
    std::string name;
    for (std::map<std::string,Histogram *>::iterator 
	   hit=m_histomap.begin();hit!=m_histomap.end();hit++) {
      histo = hit->second;
      name  = std::string("QE_Analysis/")+hit->first+std::string(".dat");
      histo->Finalize();
      histo->Output(name);
      delete histo;
    }
    m_histomap.clear();
  }
}

int Soft_Diffractive_Event_Generator::
GenerateEvent(ATOOLS::Blob_List * blobs,const bool & flag) {
  ATOOLS::Blob * blob(blobs->FindFirst(ATOOLS::btp::Soft_Collision));
  if (!blob || blob->Status()!=ATOOLS::blob_status::needs_minBias) return 0;
  if (blob->NInP()>0)  {
    msg_Error()<<"Error in "<<METHOD<<": blob has particles."<<std::endl
	       <<(*blob)<<std::endl;
    blob->DeleteInParticles();
  }
  if (blob->NOutP()>0) {
    msg_Error()<<"Error in "<<METHOD<<": blob has particles."<<std::endl
	       <<(*blob)<<std::endl;
    blob->DeleteOutParticles();
  }
  SelectModeAndFSHadrons();
  FixKinematics();
  FillBlob(blob);
  return 1;
}


void Soft_Diffractive_Event_Generator::SelectModeAndFSHadrons() {
  double disc = ran->Get(), help = disc;
  for (m_mode=0;m_mode<3;++m_mode) {
    disc -= m_rate[m_mode];
    if (disc<0.) break;
  }
  switch (m_mode) {
  case 0:
  case 1:
    m_out[m_mode]   = m_beam[m_mode];
    m_out[1-m_mode] = SelectFromMap(1-m_mode);
    break;
  case 2:
    for (size_t i=0;i<2;i++) m_out[i] = SelectFromMap(i);
    break;
  }
}

ATOOLS::Flavour Soft_Diffractive_Event_Generator::SelectFromMap(const size_t & hadron) {
  if (m_hadronmaps[hadron].empty())   return m_beam[hadron];
  if (m_hadronmaps[hadron].size()==1) return m_hadronmaps[hadron].begin()->first;
  double disc = ran->Get();
  for (std::map<ATOOLS::Flavour, double>::iterator hmit=m_hadronmaps[hadron].begin();
       hmit!=m_hadronmaps[hadron].end();hmit++) {
    disc -= hmit->second;
    if (disc<0.) return hmit->first;
  }
  return m_hadronmaps[hadron].begin()->first;
}

void Soft_Diffractive_Event_Generator::FixKinematics() {
  double Etot = m_p[0][0]+m_p[1][0], E[2];
  for (size_t i=0;i<2;i++) E[i] = (sqr(Etot)+sqr(m_out[i].Mass())+-sqr(m_out[1-i].Mass()))/(2.*Etot);
  double p = sqrt(sqr(E[0])-sqr(m_out[0].Mass())), p2 = p*p;
  m_abs_t  = p_sigma->SelectT(m_mode);
  double costheta = 1.-m_abs_t/(2.*p2), sintheta = sqrt(1.-sqr(costheta));
  double pt = p*sintheta, pt2 = sqr(pt);
  double phi(2.*M_PI*ran->Get()), ptx(pt*cos(phi)), pty(pt*sin(phi));
  double pl1(m_sign1*sqrt(p2-pt2)), pl2(-m_sign1*sqrt(p2-pt2));
  m_pout[0] = Vec4D(E[0], ptx, pty,pl1);
  m_pout[1] = Vec4D(E[1],-ptx,-pty,pl2);
}

void Soft_Diffractive_Event_Generator::FillBlob(ATOOLS::Blob * blob) {
  Particle * partin[2], * partout[2];
  for (size_t i=0;i<2;i++) {
    partin[i]  = new Particle(-1,m_beam[i],m_p[i]);
    partin[i]->SetNumber();
    partin[i]->SetBeam(i);
    partout[i] = new Particle(-1,m_out[i],m_pout[i]);
    partout[i]->SetNumber();
    blob->AddToInParticles(partin[i]);
    blob->AddToOutParticles(partout[i]);
  }
  blob->UnsetStatus(ATOOLS::blob_status::needs_minBias);
  blob->AddStatus(ATOOLS::blob_status::needs_hadrondecays);
  blob->AddStatus(ATOOLS::blob_status::needs_beams);
  blob->SetType(ATOOLS::btp::Soft_Diffractive_Collision);
  //msg_Out()<<(*blob)<<"\n";
}
