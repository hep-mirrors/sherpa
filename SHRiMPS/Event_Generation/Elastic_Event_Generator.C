#include "SHRiMPS/Event_Generation/Elastic_Event_Generator.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace SHRIMPS;

Elastic_Event_Generator::
Elastic_Event_Generator(Sigma_Elastic * sigma,const int & test) :
  Event_Generator_Base(sigma), p_sigma(sigma), 
  m_beam1(ATOOLS::rpa->gen.Beam1()), m_beam2(ATOOLS::rpa->gen.Beam2()),
  m_p1(ATOOLS::rpa->gen.PBeam(0)),   m_p2(ATOOLS::rpa->gen.PBeam(1)),
  m_E12(ATOOLS::sqr(m_p1[0])), m_pl12(Vec3D(m_p1).Sqr()), m_pl1(sqrt(m_pl12)),
  m_E22(ATOOLS::sqr(m_p2[0])), m_pl22(Vec3D(m_p2).Sqr()), m_pl2(sqrt(m_pl22)),         
  m_sign1(-1+2*int(m_p1[3]>0)),
  m_ana(true)
{
  if (m_ana) m_histomap[std::string("Q_elastic")] = new Histogram(0,0.0,1.0,1000);
  m_xsec = p_sigma->XSec();
}

Elastic_Event_Generator::~Elastic_Event_Generator() {
  if (m_ana && !m_histomap.empty()) {
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

int Elastic_Event_Generator::InitEvent(ATOOLS::Blob_List * blobs) {
  p_blob = blobs->FindFirst(ATOOLS::btp::Soft_Collision);
  if (!p_blob || p_blob->Status()!=ATOOLS::blob_status::needs_minBias) return -1;
  if (p_blob->NInP()>0)  {
    msg_Error()<<"Error in "<<METHOD<<": blob has particles.\n"<<(*p_blob)<<"\n";
    p_blob->DeleteInParticles();
  }
  if (p_blob->NOutP()>0) {
    msg_Error()<<"Error in "<<METHOD<<": blob has particles.\n"<<(*p_blob)<<"\n";
    p_blob->DeleteOutParticles();
  }
  FixKinematics();
  FillBlob();    
  return 1;
}

void Elastic_Event_Generator::FixKinematics() {
  m_abs_t = p_sigma->SelectT();
  double costheta = 1.-m_abs_t/(2.*m_pl12), sintheta = sqrt(1.-sqr(costheta));
  double pt = m_pl1*sintheta, pt2 = sqr(pt);
  double phi(2.*M_PI*ran->Get()), ptx(pt*cos(phi)), pty(pt*sin(phi));
  double pl1(m_sign1*sqrt(m_pl12-pt2)), pl2(-m_sign1*sqrt(m_pl22-pt2));
  m_p1out = Vec4D(m_p1[0],ptx,pty,pl1);
  m_p2out = Vec4D(m_p2[0],-ptx,-pty,pl2);
  if (m_ana) m_histomap[std::string("Q_elastic")]->Insert(m_abs_t);
}

void Elastic_Event_Generator::FillBlob() {
  Particle * part1in(new Particle(-1,m_beam1,m_p1));
  part1in->SetNumber();
  part1in->SetBeam(0);
  Particle * part2in(new Particle(-1,m_beam2,m_p2));
  part2in->SetNumber();
  part2in->SetBeam(1);
  Particle * part1out(new Particle(-1,m_beam1,m_p1out));
  part1out->SetNumber();
  Particle * part2out(new Particle(-1,m_beam2,m_p2out));
  part2out->SetNumber();
  
  p_blob->AddToInParticles(part1in);
  p_blob->AddToInParticles(part2in);
  p_blob->AddToOutParticles(part1out);
  p_blob->AddToOutParticles(part2out);
  p_blob->UnsetStatus(ATOOLS::blob_status::needs_minBias);
  p_blob->SetStatus(ATOOLS::blob_status::needs_beams);
  p_blob->SetType(ATOOLS::btp::Elastic_Collision);
  p_blob->AddData("Weight",new Blob_Data<double>(m_xsec));
}
