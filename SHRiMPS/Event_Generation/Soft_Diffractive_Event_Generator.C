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
  msg_Out() << "Hadron maps initialized <<<<<<-------- 1" << std::endl;
  m_histomap[std::string("Q_SD")] = new Histogram(0,0.0,1.0,1000);
  m_xsec = m_sigma;
}

void Soft_Diffractive_Event_Generator::InitialiseHadronMaps() {
  // Assume pp/ppbar collisions only
  // this can easily be extended to mesons/neutrons if necessary
  for (int i = 0; i<2; i++) {
    m_hadronmaps[i].clear();
  }
  msg_Out() << "Got into the initialise hadron maps" << std::endl;
  for (size_t i=0;i<2;i++) {
    double R = ran -> Get(); 
    m_beam[i] = (i==0)?ATOOLS::rpa->gen.Beam1():ATOOLS::rpa->gen.Beam2();
    msg_Out() << "R = " << R << std::endl;
    if (R < m_Prob1440) {
      if (m_beam[i]==ATOOLS::Flavour(kf_p_plus)) {
        m_hadronmaps[i][ATOOLS::Flavour(kf_N_1440_plus)] = 1.;
        msg_Out() << "N(1440)" << std::endl;
      }
      if (m_beam[i]==ATOOLS::Flavour(kf_p_plus).Bar()) {
        m_hadronmaps[i][ATOOLS::Flavour(kf_N_1440_plus).Bar()] = 1.;
        msg_Out() << "anti N(1440)" << std::endl;
      }
    }
    else if(R >= m_Prob1440 && R < (m_Prob1440 + m_Prob1710)) {
      if (m_beam[i]==ATOOLS::Flavour(kf_p_plus)) {
        m_hadronmaps[i][ATOOLS::Flavour(kf_N_1710_plus)] = 1.;
        msg_Out() << "!!!!!!!!!!!!!!!!!!!!!!N(1710)" << std::endl;
      }
      if (m_beam[i]==ATOOLS::Flavour(kf_p_plus).Bar()) {
        m_hadronmaps[i][ATOOLS::Flavour(kf_N_1710_plus).Bar()] = 1.;
        msg_Out() << "!!!!!!!!!!!!!!!!!!!!!!anti N(1710)" << std::endl;
      }
    }
    else {
      if (m_beam[i]==ATOOLS::Flavour(kf_p_plus)) {
        m_hadronmaps[i][ATOOLS::Flavour(kf_ud_0)] = .66;
        m_hadronmaps[i][ATOOLS::Flavour(kf_uu_1)] = .34;
        msg_Out() << "Continuous mass" << std::endl;
      }
      if (m_beam[i]==ATOOLS::Flavour(kf_p_plus).Bar()) {
        m_hadronmaps[i][ATOOLS::Flavour(kf_ud_0).Bar()] = .66;
        m_hadronmaps[i][ATOOLS::Flavour(kf_uu_1).Bar()] = .34;
        msg_Out() << "anti Continuous mass" << std::endl;
      }
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
  ATOOLS::Flavour no_particle;
  for (int i = 0; i < 4; i++) m_out[i] = no_particle;
  ATOOLS::Vec4D no_4vec;
  for (int i = 0; i < 4; i++) m_pout[i] = no_4vec;
  msg_Out() << "Hadron maps to be initialized again <<<<<<-------- 2" << std::endl;
  InitialiseHadronMaps();
  SelectModeAndFSHadrons();
  FixKinematics();
  FillBlob(blob);
  return 1;
}


void Soft_Diffractive_Event_Generator::SelectModeAndFSHadrons() {
  m_contMassRange0 = 0;
  m_contMassRange1 = 0;
  double disc = ran->Get(), help = disc;
  for (m_mode=0;m_mode<3;++m_mode) {
    disc -= m_rate[m_mode];
    if (disc<0.) break;
  }
  msg_Out() << "mode:   " << m_mode << std::endl;
  switch (m_mode) {
  case 0:
    m_out[m_mode]   = m_beam[m_mode];
    m_out[1-m_mode] = SelectFromMap(1-m_mode);
    if (m_out[1-m_mode] == ATOOLS::Flavour(kf_ud_0)){
      m_out[1-m_mode+2] = ATOOLS::Flavour(kf_u);
      m_contMassRange1 = 1;
    }
    else if(m_out[1-m_mode] == ATOOLS::Flavour(kf_ud_0).Bar()){
      m_out[1-m_mode+2] = ATOOLS::Flavour(kf_u).Bar();
      m_contMassRange1 = 1;
    }
    else if (m_out[1-m_mode] == ATOOLS::Flavour(kf_uu_1)){
      m_out[1-m_mode+2] = ATOOLS::Flavour(kf_d);
      m_contMassRange1 = 1;
    }
    else if(m_out[1-m_mode] == ATOOLS::Flavour(kf_uu_1).Bar()){
      m_out[1-m_mode+2] = ATOOLS::Flavour(kf_d).Bar();
      m_contMassRange1 = 1;
    }
    //msg_Out() << " Single Diff: "<< (1-m_mode) << std::endl;
    //msg_Out() << "here! " << m_out[0] << "\t" << m_out[1] << "\t" << m_out[2] << "\t" << m_out[3] << std::endl;
    break;
  case 1:
    m_out[m_mode]   = m_beam[m_mode];
    m_out[1-m_mode] = SelectFromMap(1-m_mode);
    msg_Out() << m_out[1-m_mode] << std::endl;
    if (m_out[1-m_mode] == ATOOLS::Flavour(kf_ud_0)){
      m_out[1-m_mode+2] = ATOOLS::Flavour(kf_u);
      m_contMassRange0 = 1;
    }
    else if(m_out[1-m_mode] == ATOOLS::Flavour(kf_ud_0).Bar()){
      m_out[1-m_mode+2] = ATOOLS::Flavour(kf_u).Bar();
      m_contMassRange0 = 1;
    }
    else if (m_out[1-m_mode] == ATOOLS::Flavour(kf_uu_1)){
      m_out[1-m_mode+2] = ATOOLS::Flavour(kf_d);
      m_contMassRange0 = 1;
    }
    else if(m_out[1-m_mode] == ATOOLS::Flavour(kf_uu_1).Bar()){
      m_out[1-m_mode+2] = ATOOLS::Flavour(kf_d).Bar();
      m_contMassRange0 = 1;
    }
    //msg_Out() << " Single Diff: "<< (1-m_mode) << std::endl;
    //msg_Out() << "here! " << m_out[0] << "\t" << m_out[1] << "\t" << m_out[2] << "\t" << m_out[3] << std::endl;
    break;
  case 2:
    for (size_t i=0;i<2;i++) {
      m_out[i] = SelectFromMap(i);
      msg_Out() << m_out[i] << std::endl;
      if (m_out[i] == ATOOLS::Flavour(kf_ud_0)){
        m_out[i+2] = ATOOLS::Flavour(kf_u);
        if (i == 0) m_contMassRange0 = 1;
        else m_contMassRange1 = 1;
      }
      else if(m_out[i] == ATOOLS::Flavour(kf_ud_0).Bar()){
        m_out[i+2] = ATOOLS::Flavour(kf_u).Bar();
        if (i == 0) m_contMassRange0 = 1;
        else m_contMassRange1 = 1;
      }
      else if (m_out[i] == ATOOLS::Flavour(kf_uu_1)){
        m_out[i+2] = ATOOLS::Flavour(kf_d);
        if (i == 0) m_contMassRange0 = 1;
        else m_contMassRange1 = 1;
      }
      else if(m_out[i] == ATOOLS::Flavour(kf_uu_1).Bar()){
        m_out[i+2] = ATOOLS::Flavour(kf_d).Bar();
        if (i == 0) m_contMassRange0 = 1;
        else m_contMassRange1 = 1;
      }
    }
    //msg_Out() << " Double diffractive " << std::endl;
    //msg_Out() << m_out[0] << "\t" << m_out[1] << "\t" << m_out[2] << "\t" << m_out[3] << std::endl;
    break;
  }
}

ATOOLS::Flavour Soft_Diffractive_Event_Generator::SelectFromMap(const size_t & hadron) {
  if (m_hadronmaps[hadron].empty()){
    //std::cout << "i = " << hadron << "\t" << m_beam[hadron] <<std::endl;
    return m_beam[hadron];
  }
  if (m_hadronmaps[hadron].size()==1) {
    std::cout << "i = " << hadron << "\t" << m_hadronmaps[hadron].begin()->first <<std::endl;
    return m_hadronmaps[hadron].begin()->first;
  }
  double disc = ran->Get();
  for (std::map<ATOOLS::Flavour, double>::iterator hmit=m_hadronmaps[hadron].begin();
       hmit!=m_hadronmaps[hadron].end();hmit++) {
    disc -= hmit->second;
    if (disc<0.) {
      //std::cout << "i = " << hadron << "\t" << hmit->first <<std::endl;
      return hmit->first;
    }
  }
  //std::cout << "i = " << hadron << "\t" << m_hadronmaps[hadron].begin()->first <<std::endl;
  return m_hadronmaps[hadron].begin()->first;
}

std::vector<ATOOLS::Vec4D> Soft_Diffractive_Event_Generator::SplitIntoQandQQ(ATOOLS::Vec4D pmu, double Mqq, double Mq) {
  ATOOLS::Poincare rot(pmu,Vec4D(0.,0.,0.,1.));
  rot.Rotate(pmu);
  double x =  2./3. + (ran -> GetGaussian())/10.;
  while (x > 1) {
    x =  2./3. + (ran -> GetGaussian())/10.;
  }
  double Ep = pmu[0], pz = pmu[3];
  double PprimeAbs = sqrt(sqr(x*Ep) - sqr(Mqq));
  double cosTheta = (sqr(Mq) - sqr((1-x)*Ep) + sqr(PprimeAbs) + sqr(pz))/(2*pz*PprimeAbs);
  if (cosTheta > 1) {
    msg_Out() << "cosine is larger than one" << std::endl;
  }
  double sinTheta = sqrt(1 - sqr(cosTheta));
  double phi = 2.*M_PI*ran->Get();
  double PprimeX = PprimeAbs*cos(phi)*sinTheta;
  double PprimeY = PprimeAbs*sin(phi)*sinTheta;
  double PprimeZ = PprimeAbs*cosTheta;
  std::vector<ATOOLS::Vec4D> pprime;
  ATOOLS::Vec4D Pp = Vec4D(x*Ep,PprimeX,PprimeY,PprimeZ);
  ATOOLS::Vec4D Kp = Vec4D((1-x)*Ep,-PprimeX,-PprimeY,pz-PprimeZ);
  rot.RotateBack(Pp);
  rot.RotateBack(Kp);
  pprime.push_back(Pp); //di quark
  pprime.push_back(Kp); //quark
  return pprime;
}

double Soft_Diffractive_Event_Generator::drawMass(){
  double Mcont_sqr_min = pow(2.,2); // GeV^2 
  double Mcont_sqr_max = pow(20.,2); // GeV^2
  double b = .25;
  double a = b*(-1 + m_Prob1710 + m_Prob1440)/(exp(-b*Mcont_sqr_max) - exp(-b*Mcont_sqr_min));
  double R = ran -> Get();
  double Msqr = (-1/b)*log((1-R)*exp(-b*Mcont_sqr_min)+R*exp(-b*Mcont_sqr_max));
  return Msqr;
}

std::vector<double> Soft_Diffractive_Event_Generator::ComputePxPyPz(double p, int sign1, int mode) {
  std::vector<double> momentum_vec;
  double p2 = p*p;
  m_abs_t = p_sigma->SelectT(mode);
  double costheta = 1.-m_abs_t/(2.*p2), sintheta = sqrt(1.-sqr(costheta));
  double pt = p*sintheta, pt2 = sqr(pt);
  double phi(2.*M_PI*ran->Get()), ptx(pt*cos(phi)), pty(pt*sin(phi));
  double pl1(sign1*sqrt(p2-pt2)), pl2(-sign1*sqrt(p2-pt2));
  momentum_vec.push_back(ptx);
  momentum_vec.push_back(pty);
  momentum_vec.push_back(pl1);
  return momentum_vec;
}

ATOOLS::Vec4D Soft_Diffractive_Event_Generator::Get4Vector(double M2[], double Etot){
  std::vector<ATOOLS::Vec4D> vec_4vec;
  double E[2];
  for (size_t i=0;i<2;i++) E[i] = (sqr(Etot)+M2[i]+-M2[1-i])/(2.*Etot);
  double p = sqrt(sqr(E[0])-M2[0]);
  std::vector<double> momentum_vec = ComputePxPyPz(p, m_sign1, m_mode);
  vec_4vec.push_back(Vec4D(E[0], momentum_vec[0], momentum_vec[1], momentum_vec[2]));
  vec_4vec.push_back(Vec4D(E[1],-momentum_vec[0],-momentum_vec[1],-momentum_vec[2]));
  return vec_4vec[0];
}

void Soft_Diffractive_Event_Generator::FixKinematics() {
  msg_Out() << "Fix Kinematics" << std::endl;
  msg_Out() << m_out[0] << "\t" << m_out[1] << "\t" << m_out[2] << "\t" << m_out[3] << std::endl;
  double Etot = m_p[0][0]+m_p[1][0];
  double M2[2];
  ATOOLS::Vec4D pmu[4];
  if (m_contMassRange0 == 0 && m_contMassRange1 == 0) {
    msg_Out() << "\nNone of them in the Continuous mass range !!! " << std::endl;
    msg_Out() << m_out[0] << "\t" << m_out[1] << std::endl;
    for (size_t i=0;i<2;i++) M2[i] = sqr(m_out[i].Mass());
    pmu[0] = Get4Vector(M2, Etot);
    pmu[1] = Vec4D(Etot-pmu[0][0],-pmu[0][1],-pmu[0][2],-pmu[0][3]);
    m_pout[0] = pmu[0];
    m_pout[1] = pmu[1];
  }
  else if(m_contMassRange0 ^ m_contMassRange1){
    msg_Out() << "\nEither Particle 0 OR Particle 1 is in the Continuous mass range !!! " << std::endl;
    if (m_contMassRange0 == 1) m_contMassIndex = 0;
    else m_contMassIndex = 1;
    M2[m_contMassIndex] = drawMass();
    M2[1-m_contMassIndex] = sqr(m_out[1-m_contMassIndex].Mass());
    double Mqq = m_out[m_contMassIndex].Mass(), Mq = m_out[m_contMassIndex+2].Mass();
    pmu[0] = Get4Vector(M2, Etot);
    pmu[1] = Vec4D(Etot-pmu[0][0],-pmu[0][1],-pmu[0][2],-pmu[0][3]);
    std::vector<ATOOLS::Vec4D> split = SplitIntoQandQQ(pmu[m_contMassIndex],Mqq,Mq);
    m_pout[m_contMassIndex] = split[0]; //di quark
    m_pout[m_contMassIndex+2] = split[1]; //quark
    m_pout[1-m_contMassIndex] = pmu[1-m_contMassIndex];
  }
  else {
    msg_Out() << "\nBoth of them in the Continuous mass !!! " << std::endl;
    for (size_t i = 0;i<2;i++) M2[i] = sqr(drawMass());
    pmu[0] = Get4Vector(M2, Etot);
    pmu[1] = Vec4D(Etot-pmu[0][0],-pmu[0][1],-pmu[0][2],-pmu[0][3]);
    double finalMasses[4];
    for (size_t i = 0;i<4;i++) finalMasses[i] = m_out[i].Mass();
    std::vector<ATOOLS::Vec4D> splits[2];
    for (int i = 0; i<2; i++){
      splits[i] = SplitIntoQandQQ(pmu[i],finalMasses[i],finalMasses[i+2]);
      m_pout[i] = splits[i][0];
      m_pout[i+2] = splits[i][1];
    }
  }
}

void Soft_Diffractive_Event_Generator::FillBlob(ATOOLS::Blob * blob) {
  msg_Out() << "FILLING BLOB!!!!" << std::endl;
  msg_Out() << "m_mode = " << m_mode << std::endl;
  Particle * partin[2], * partout[4]; 
  for (size_t i=0;i<2;i++) {
    partin[i]  = new Particle(-1,m_beam[i],m_p[i]);
    partin[i]->SetNumber();
    partin[i]->SetBeam(i);
    partin[i]->SetInfo('I');
    partout[i] = new Particle(-1,m_out[i],m_pout[i]);
    partout[i]->SetNumber();
    partout[i]->SetInfo('F');
    blob->AddToInParticles(partin[i]);
    blob->AddToOutParticles(partout[i]);
    msg_Out() << i <<  "\t" << m_out[i] << "\t" << m_pout[i] << std::endl;
  }
  if (m_contMassRange0 ^ m_contMassRange1){
    partout[m_contMassIndex+2] = new Particle(-1,m_out[m_contMassIndex+2],m_pout[m_contMassIndex+2]);
    partout[m_contMassIndex+2]->SetNumber();
    partout[m_contMassIndex+2]->SetFlow(1,-1);
    partout[m_contMassIndex+2]->SetInfo('F');
    blob->OutParticle(m_contMassIndex)->SetFlow(2,partout[m_contMassIndex+2]->GetFlow(1));
    msg_Out() << *(blob -> OutParticle(m_contMassIndex)) << std::endl;
    blob->AddToOutParticles(partout[m_contMassIndex+2]);
    blob->AddStatus(ATOOLS::blob_status::needs_hadronization);
    msg_Out() << m_contMassIndex+2 <<  "\t" << m_out[m_contMassIndex+2] << "\t" << m_pout[m_contMassIndex+2] << std::endl;
  }
  else if (m_contMassRange0 == 1 && m_contMassRange1 == 1) {
    for (int i = 2; i < 4; i++){
      partout[i] = new Particle(-1,m_out[i],m_pout[i]);
      partout[i]->SetNumber();
      partout[i]->SetFlow(1,-1);
      blob->OutParticle(i-2)->SetFlow(2,partout[i]->GetFlow(1));
      blob->AddToOutParticles(partout[i]);
      msg_Out() << i <<  "\t" << m_out[i] << "\t" << m_pout[i] << std::endl;
    }
    blob->AddStatus(ATOOLS::blob_status::needs_hadronization);
}
  msg_Out() << (m_pout[0] + m_pout[2]).Abs2() << std::endl;
  blob->UnsetStatus(ATOOLS::blob_status::needs_minBias);
  blob->AddStatus(ATOOLS::blob_status::needs_hadrondecays);
  blob->AddStatus(ATOOLS::blob_status::needs_beams);
  blob->SetType(ATOOLS::btp::Soft_Diffractive_Collision);
  msg_Out() << "something" << std::endl;
}
