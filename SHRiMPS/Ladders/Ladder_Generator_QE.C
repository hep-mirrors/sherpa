#include "SHRiMPS/Ladders/Ladder_Generator_QE.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <list>
#include "SHRiMPS/Cross_Sections/Sigma_Elastic.H"
#include "SHRiMPS/Cross_Sections/Sigma_Partonic.H"
#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"

using namespace SHRIMPS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Ladder_Generator_QE::Ladder_Generator_QE() : 
  Ladder_Generator_Base(),
  m_beam1(ATOOLS::rpa->gen.Beam1()), m_beam2(ATOOLS::rpa->gen.Beam2()),
  m_Pbeam1((ATOOLS::rpa->gen.PBeam(0))), m_Pbeam2(ATOOLS::rpa->gen.PBeam(1)),
  m_sign1(-1+2*int(m_Pbeam1[3]>0)), m_fraction(.1),
  m_partonic(Sigma_Partonic(xs_mode::perturbative))
  {}
  
Ladder_Generator_QE::~Ladder_Generator_QE() {}

Ladder * Ladder_Generator_QE::operator()(const Vec4D & pos,Sigma_Elastic * sigma_el,Sigma_D * sigma_d) {
  p_sigma_el = sigma_el;
  p_sigma_d = sigma_d;
  InitLadder(pos);
  for (int beam=0;beam<2;beam++) m_ladder_ends[beam] = ATOOLS::Flavour(kf_gluon);
  FillGluons();
  SelectPropagatorColours();
  if (p_ladder->GetEmissions()->size()==2) {
    m_interactionType = GetIntType();
    FixEmissionsKinematics_elastic();
    //ConstructSimpleLadder();
    ConstructISKinematics();
  }
  else {
    msg_Out() << "number of emissions is not 2!" << endl;
  }
  msg_Out() << *p_ladder << endl;
  return p_ladder;
}

void Ladder_Generator_QE::FixEmissionsKinematics_SD(int mode) {
  m_p1 = m_Pbeam1*m_fraction;
  m_p2 = m_Pbeam2*m_fraction;
  double etot(m_p1[0]+m_p2[0]);
}

void Ladder_Generator_QE::FixEmissionsKinematics_elastic() {
  m_p1 = m_Pbeam1*m_fraction;
  msg_Out() << "-----> p1 = " << m_p1 << endl;
  m_p2 = m_Pbeam2*m_fraction;
  msg_Out() << "-----> p2 = " << m_p2 << endl;
  m_pl12 = Vec3D(m_p1).Sqr();
  m_pl22 = Vec3D(m_p2).Sqr();
  m_pl1 = sqrt(m_pl12);
  m_pl2 = sqrt(m_pl22);
  switch (m_interactionType) {
    case 0: m_abs_t = p_sigma_d->SelectT(0);
    case 1: m_abs_t = p_sigma_d->SelectT(1);
    case 2: m_abs_t = p_sigma_d->SelectT(2);
    case 4: m_abs_t = p_sigma_el->SelectT();
  }
  //std::ofstream tvals;
  //tvals.open("./tvals_bias_large.txt", std::ios::app);
  //for (int i = 0; i < 1000; i++) {
  //  tvals << 0 << "\t" << p_sigma_d->SelectT(0) << endl;
  //  tvals << 1 << "\t" << p_sigma_d->SelectT(1) << endl;
  //  tvals << 2 << "\t" << p_sigma_d->SelectT(2) << endl;
  //  tvals << 4 << "\t" << p_sigma_el->SelectT() << endl;
  //}
  //tvals.close();
  double costheta = 1.-m_abs_t/(2.*m_pl12), sintheta = sqrt(1.-sqr(costheta));
  double pt = m_pl1*sintheta, pt2 = sqr(pt);
  double phi(2.*M_PI*ran->Get()), ptx(pt*cos(phi)), pty(pt*sin(phi));
  double pl1(m_sign1*sqrt(m_pl12-pt2)), pl2(-m_sign1*sqrt(m_pl22-pt2));
  double E1 = sqrt(sqr(ptx) + sqr(pty) + sqr(pl1));
  double E2 = sqrt(sqr(-ptx) + sqr(-pty) + sqr(pl2));
  msg_Out() << E1 << "\t" << m_p1[0] << endl;
  msg_Out() << E2 << "\t" << m_p2[0] << endl;
  p_emissions->begin()->second.SetMomentum(Vec4D(E1,ptx,pty,pl1));
  p_emissions->rbegin()->second.SetMomentum(Vec4D(E2,-ptx,-pty,pl2));
  Vec4D Pgluon1 = p_emissions->begin()->second.Momentum();
  Vec4D Pgluon2 = p_emissions->rbegin()->second.Momentum();
  Vec4D remnantP1 = (1. - m_fraction)*m_Pbeam1;
  Vec4D remnantP2 = (1. - m_fraction)*m_Pbeam2;
  Vec4D beam0stuff = remnantP1 + Pgluon1;
  Vec4D beam1stuff = remnantP2 + Pgluon2;
  //msg_Out() << "total 4momentum from beam 0: " << beam0stuff << endl;
  //msg_Out() << "inv. mass2: " << beam0stuff.Abs2() << endl;
  //msg_Out() << "gluon 4momentum from beam 0: " << Pgluon1 << endl;
  //msg_Out() << "inv. mass2: " << Pgluon1.Abs2() << endl;
  //msg_Out() << "remnant 4momentum from beam 0: " << remnantP1 << endl;
  //msg_Out() << "inv. mass2: " << remnantP1.Abs2() << endl;

  //msg_Out() << "total 4momentum from beam 1: " << beam1stuff << endl;
  //msg_Out() << "inv. mass2: " << beam1stuff.Abs2() << endl;
  //msg_Out() << "gluon 4momentum from beam 1: " << Pgluon2 << endl;
  //msg_Out() << "inv. mass2: " << Pgluon2.Abs2() << endl;
  //msg_Out() << "remnant 4momentum from beam 1: " << remnantP2 << endl;
  //msg_Out() << "inv. mass2: " << remnantP2.Abs2() << endl;

  //PRINT_VAR(sqrt(m_abs_t));
  //if (m_ana) m_histomap[std::string("Q_elastic")]->Insert(m_abs_t);
  T_Prop & prop = *p_props->begin();
  prop.SetQT2(m_qt2);
  PRINT_VAR(m_qt2);
  prop.SetQ02(m_qt2min);
  PRINT_VAR(m_qt2min);
  prop.SetQ(sqrt(m_qt2)*m_eqt);
}

void Ladder_Generator_QE::FillGluons() {
  for (size_t beam=0;beam<2;beam++) {
    m_ylimits[beam] = (beam==0? 1. : -1.) * (m_Ymax + ran->Get()*m_deltaY);
    p_ladder->AddRapidity(m_ylimits[beam],m_ladder_ends[beam]);
  }
}

void Ladder_Generator_QE::SelectPropagatorColours() {
  for (size_t i=0;i<p_emissions->size()-1;i++)
    p_props->push_back(T_Prop(colour_type::singlet,Vec4D(0.,0.,0.,0.),m_qt2min));
}

void Ladder_Generator_QE::CalculateWeight() {
  m_weight  = AlphaSWeight(m_lastk.PPerp2());
  m_weight *= m_me(p_ladder,m_qt2min); 
}
