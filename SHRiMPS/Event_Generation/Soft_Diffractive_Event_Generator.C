#include "SHRiMPS/Event_Generation/Soft_Diffractive_Event_Generator.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace SHRIMPS;
using namespace ATOOLS;

Soft_Diffractive_Event_Generator::
Soft_Diffractive_Event_Generator(Sigma_D * sigma,const int & test) :
  Event_Generator_Base(sigma),
  p_sigma(sigma), m_sigma(0.),
  m_massexp(0.5), m_Prob1440(.4),  m_Prob1710(0.2), m_runmode(MBpars.RunMode())
{
  msg_Out() << "RUN MODE: " << m_runmode << "\n";
  for (size_t i=0;i<3;i++) m_sigma   += m_rate[i] = p_sigma->GetXSec(i);
  for (size_t i=0;i<3;i++) {
    m_rate[i] /= m_sigma;
  }
  for (size_t i=0;i<2;i++) {
    m_p[i]        = rpa->gen.PBeam(i);
    m_E[i]        = m_p[i][0];
    m_beam[i]     = i==0? rpa->gen.Beam1() : rpa->gen.Beam2();
    m_minmass2[i] = sqr(m_beam[i].HadMass()+0.5);
    m_maxmass2[i] = sqr(20.);
  }
  m_sign1 = -1+2*int(m_p[0][3]>0);
  m_histomap[std::string("Q_SD")] = new Histogram(0,0.0,1.0,1000);
  m_xsec = m_sigma;
  InitialiseHadronMaps();
}

void Soft_Diffractive_Event_Generator::InitialiseHadronMaps() {
  // Assume pp/ppbar collisions only
  // this can easily be extended to mesons/neutrons if necessary
  // flavour wave function for proton (and its resonances):
  // |p> = sqrt(1/3) | d+uu_1> + sqrt(1/6) | u+ud_1> + sqrt(1/2) | u+ud_0>
  // need to square the weights for the probabilities - so ratio is 2:1:3.
  std::pair<Flavour, Flavour> flpair;
  for (size_t i=0;i<2;i++) {
    if (m_beam[i]==Flavour(kf_p_plus) ||
	m_beam[i]==Flavour(kf_p_plus).Bar()) {
      double rest = 1.;
      bool barred = (m_beam[i]==Flavour(kf_p_plus).Bar());
      if (m_Prob1440>1.e-12) {
	flpair.first  = barred?Flavour(kf_N_1440_plus).Bar():Flavour(kf_N_1440_plus);
	flpair.second = Flavour(kf_none);
	rest         -= m_hadronmaps[i][flpair] = m_Prob1440;
	m_minmass2[i] = Max(m_minmass2[i], sqr(flpair.first.HadMass()+2.*flpair.first.Width())); 
      }
      if (rest>1.e-6 && m_Prob1710>1.e-12) {
	flpair.first  = barred?Flavour(kf_N_1710_plus).Bar():Flavour(kf_N_1710_plus);
	flpair.second = Flavour(kf_none);
	rest         -= m_hadronmaps[i][flpair] = m_Prob1710;
	m_minmass2[i] = Max(m_minmass2[i], sqr(flpair.first.HadMass()+2.*flpair.first.Width())); 
      }
      if (rest>1.e-6) {
	double cont   = rest;
	flpair.first  = barred?Flavour(kf_d).Bar():Flavour(kf_d);
	flpair.second = barred?Flavour(kf_uu_1).Bar():Flavour(kf_uu_1);
	rest         -= m_hadronmaps[i][flpair] = cont*(1./3.);
	flpair.first  = barred?Flavour(kf_u).Bar():Flavour(kf_u);
	flpair.second = barred?Flavour(kf_ud_1).Bar():Flavour(kf_ud_1);
	rest         -= m_hadronmaps[i][flpair] = cont*(1./6.);
	flpair.first  = barred?Flavour(kf_u).Bar():Flavour(kf_u);
	flpair.second = barred?Flavour(kf_ud_0).Bar():Flavour(kf_ud_0);
	rest         -= m_hadronmaps[i][flpair] = rest;
	m_minmass2[i] = Max(m_minmass2[i], sqr(flpair.first.HadMass()+flpair.second.HadMass()+0.5));
      }
    }
    if (sqrt(m_maxmass2[i])<=sqrt(m_minmass2[i])+1.)
      m_maxmass2[i] = sqr(sqrt(m_minmass2[i])+1.);
    m_expargLow[i] = exp(-m_massexp*m_minmass2[i]);
    m_expargUp[i]  = exp(-m_massexp*m_maxmass2[i]);
  }

  for (size_t i=0; i<2; i++) {
    msg_Info()<<"Diffractive modes for "<<m_beam[i]<<":\n";
    for (std::map<std::pair<Flavour, Flavour>, double>::iterator hmit=m_hadronmaps[i].begin();
	 hmit!=m_hadronmaps[i].end();hmit++) {
      msg_Info()<<"  --> "<<hmit->first.first<<" + "<<hmit->first.second<<"  with P = "<<hmit->second<<"\n";
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

int Soft_Diffractive_Event_Generator::InitEvent(ATOOLS::Blob_List * blobs) {
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
  for (size_t i=0;i<4;i++) {
    m_out[i] = Flavour(kf_none); m_pout[i] = Vec4D(0.,0.,0.,0.);
  }
  for (size_t i=0;i<2;i++) m_contMassRange[i] = false;
  SelectMode();
  msg_Out() << "MODE PICKED = " << m_mode << "  " << int(m_runmode) << std::endl;
  SelectFS();
  FixKinematics();
  FillBlob();
  return 1;
}

void Soft_Diffractive_Event_Generator::SelectMode() {
  double disc = ran->Get();
  switch(m_runmode) {
    case run_mode::soft_diffractive_events:
      for (m_mode=0;m_mode<3;++m_mode) {
        disc -= m_rate[m_mode];
        if (disc<=0.) break;
      }
    case run_mode::single_diffractive_events:
      m_mode = 1;
      break;
    case run_mode::double_diffractive_events:
      m_mode = 2;
      break;
  }
}

void Soft_Diffractive_Event_Generator::SelectFS() {
  switch (m_mode) {
  case 0:
  case 1:
    // the two single-diffractive modes
    m_out[2*m_mode]   = m_beam[m_mode];
    SelectFlavours(1-m_mode);
    break;
  case 2:
    for (size_t beam=0;beam<2;beam++) SelectFlavours(beam);
    break;
  }
}

void Soft_Diffractive_Event_Generator::SelectFlavours(size_t beam) {
  double disc = 0.9999999999*ran->Get();
  for (std::map<std::pair<Flavour, Flavour>, double>::iterator hmit=m_hadronmaps[beam].begin();
       hmit!=m_hadronmaps[beam].end(); hmit++) {
    disc-=hmit->second;
    if (disc<=0.) {
      m_out[2*beam]   = hmit->first.first;
      m_out[2*beam+1] = hmit->first.second;
      if (m_out[2*beam+1]!=Flavour(kf_none)) m_contMassRange[beam] = true;
      return;
    }
  }
  msg_Error()<<METHOD<<" should never arrive here - this is a safety measure.\n";
  m_out[2*beam]   = m_hadronmaps[beam].begin()->first.first;
  m_out[2*beam+1] = m_hadronmaps[beam].begin()->first.second;
  if (m_out[2*beam+1]!=Flavour(kf_none)) m_contMassRange[beam] = true;
}

void Soft_Diffractive_Event_Generator::FixKinematics() {
  Vec4D outMom[2];
  FixBinarySystem(outMom);
  for (size_t beam=0;beam<2;beam++) {
    if (m_contMassRange[beam]) SplitQandQQ(beam,outMom[beam]);
    else {
      m_pout[2*beam]   = outMom[beam];
      m_pout[2*beam+1] = Vec4D(0.,0.,0.,0.);
    }
    //msg_Out()<<METHOD<<": mom = "<<outMom[beam]<<" ("<<sqrt(outMom[beam].Abs2())<<")-> "
    //	     <<m_out[2*beam]<<" + "<<m_out[2*beam+1]<<"\n"
    //	     <<"--> "<<m_pout[2*beam]<<" ("<<sqrt(m_pout[2*beam].Abs2())<<") + "
    //	     <<m_pout[2*beam+1]<<" ("<<sqrt(m_pout[2*beam+1].Abs2())<<")\n";
  }
}

void Soft_Diffractive_Event_Generator::FixBinarySystem(Vec4D (& moms)[2]) {
  double absT = p_sigma->SelectT(m_mode), Q2 = (m_p[0]+m_p[1]).Abs2(), Q = sqrt(Q2), outMass2[2];
  for (size_t beam=0;beam<2;beam++) {
    outMass2[beam] = m_contMassRange[beam] ? SelectMass2(beam) : sqr(m_out[2*beam].HadMass());
  }
  double E[2];
  for (size_t beam=0;beam<2;beam++) E[beam] = (Q2+outMass2[beam]-outMass2[1-beam])/(2.*Q);
  double p   = sqrt(sqr(E[0])-outMass2[0]), pt = sqrt(absT), pl = sqrt(sqr(p)-absT);
  double phi = 2.*M_PI*ran->Get();
  moms[0] = Vec4D(E[0],pt*cos(phi),pt*sin(phi),m_sign1*pl);
  moms[1] = Vec4D(Q,0.,0.,0.)-moms[0];
  //msg_Out()<<METHOD<<"(t = "<<absT<<"):\n";
  //for (size_t beam=0;beam<2;beam++) {
  //  msg_Out()<<"   "<<m_beam[beam]<<" --> "<<m_out[2*beam]<<" + "<<m_out[2*beam+1]<<", "
  //	     <<"continuous masses = "<<m_contMassRange[beam]
  //	     <<" --> "<<sqrt(outMass2[beam])
  //	     <<" --> "<<moms[beam]<<"  ("<<sqrt(moms[beam].Abs2())<<")\n";
  //}
}

void Soft_Diffractive_Event_Generator::SplitQandQQ(size_t beam,ATOOLS::Vec4D & mom) {
  double Q2  = mom.Abs2(), Q = sqrt(Q2);
  double m12 = sqr(m_out[2*beam].HadMass()), m22 = sqr(m_out[2*beam+1].HadMass());
  double E1  = (Q2+m12-m22)/(2.*Q), p = sqrt(sqr(E1)-m12);
  double costheta  = 1.-2.*ran->Get(), sintheta = sqrt(1.-sqr(costheta)), phi = 2.*M_PI*ran->Get();
  m_pout[2*beam]   = Vec4D(E1,p*sintheta*cos(phi),p*sintheta*sin(phi),p*costheta);
  Poincare rest    = Poincare(mom);
  rest.BoostBack(m_pout[2*beam]);
  m_pout[2*beam+1] = mom-m_pout[2*beam];
}

double Soft_Diffractive_Event_Generator::SelectMass2(size_t beam) {
  return -1./m_massexp* log(m_expargLow[beam] + ran->Get()*(m_expargUp[beam]-m_expargLow[beam]));
}

void Soft_Diffractive_Event_Generator::FillBlob() {
  Particle * part;
  for (size_t beam=0;beam<2;beam++) {
    part = new Particle(-1,m_beam[beam],m_p[beam]);
    part->SetNumber();
    part->SetBeam(beam);
    part->SetInfo('I');
    p_blob->AddToInParticles(part);
  }
  Vec4D Pfrombeams[2];
  for (size_t beam=0;beam<2;beam++) {
    if (m_contMassRange[beam]) {
      msg_Out()<<"  - "<<METHOD<<"(beam = "<<beam<<") selected continuous mass range.\n";
      p_blob->AddStatus(blob_status::needs_hadronization);
      for (size_t j=0;j<2;j++) {
	      part = new Particle(-1,m_out[2*beam+j],m_pout[2*beam+j]);
	      part->SetNumber();
	      part->SetBeam(beam);
	      part->SetInfo('F');
	    if ((m_out[2*beam+j].IsQuark() && !m_out[2*beam+j].IsAnti()) ||
	        (m_out[2*beam+j].IsDiQuark() && m_out[2*beam+j].IsAnti()))
	      part->SetFlow(1,500+beam);
	    else if ((m_out[2*beam+j].IsQuark() && m_out[2*beam+j].IsAnti()) ||
		          (m_out[2*beam+j].IsDiQuark() && !m_out[2*beam+j].IsAnti()))
	      part->SetFlow(2,500+beam);
	    p_blob->AddToOutParticles(part);
      Pfrombeams[beam]+=part->Momentum();
      //msg_Out() << part->Momentum() << "\t" << (part->Momentum() ).Abs2() << "\n";
      }
    }
    else {
      part = new Particle(-1,m_out[2*beam],m_pout[2*beam]);
      part->SetNumber();
      part->SetBeam(beam);
      part->SetInfo('F');
      Pfrombeams[beam]+=part->Momentum();
      //msg_Out() << part->Momentum() << "\t" << (part->Momentum() ).Abs2() << "\n";
      p_blob->AddToOutParticles(part);
    }
    //msg_Out() << "beam " << beam << "\t" << Pfrombeams[beam] << "\t" << Pfrombeams[beam].Abs2() << "\n";
  }
  p_blob->UnsetStatus(blob_status::needs_minBias);
  if (!p_blob->Has(blob_status::needs_hadronization))
    p_blob->AddStatus(blob_status::needs_hadrondecays);
  p_blob->AddStatus(blob_status::needs_beams);
  p_blob->SetType(btp::Soft_Diffractive_Collision);
}

/*
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


std::vector<double> Soft_Diffractive_Event_Generator::ComputePxPyPz(double p, int sign1, int mode) {
  std::vector<double> momentum_vec;
  double p2 = p*p;
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

*/
