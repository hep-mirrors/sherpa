#include "SHRiMPS/Ladders/Ladder_Generator_QE.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <list>
#include "SHRiMPS/Cross_Sections/Sigma_Elastic.H"
#include "SHRiMPS/Cross_Sections/Sigma_Partonic.H"

using namespace SHRIMPS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Ladder_Generator_QE::Ladder_Generator_QE() : 
  Ladder_Generator_Base(),
  m_beam1(ATOOLS::rpa->gen.Beam1()), m_beam2(ATOOLS::rpa->gen.Beam2()),
  m_Pbeam1((ATOOLS::rpa->gen.PBeam(0))), m_Pbeam2(ATOOLS::rpa->gen.PBeam(1)),
  m_sign1(-1+2*int(m_Pbeam1[3]>0)), m_fraction(.1), m_pdf_over_estimate(2.5)
  {}
  
//Ladder_Generator_QE::~Ladder_Generator_QE() {}

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
  //msg_Out() << *p_ladder << endl;
  return p_ladder;
}

void Ladder_Generator_QE::FixEmissionsKinematics_SD(int mode) {
  m_p1 = m_Pbeam1*m_fraction;
  m_p2 = m_Pbeam2*m_fraction;
  double etot(m_p1[0]+m_p2[0]);
}

void Ladder_Generator_QE::FixEmissionsKinematics_elastic() {
  switch (m_interactionType) {
    case 0: m_abs_t = p_sigma_d->SelectT(0);
    case 1: m_abs_t = p_sigma_d->SelectT(1);
    case 2: m_abs_t = p_sigma_d->SelectT(2);
    case 4: m_abs_t = p_sigma_el->SelectT();
  }  
  for(size_t beam = 0; beam < 2; beam++) m_x[beam] = 0.;
  for(size_t beam = 0; beam < 2; beam++) {
    m_x[beam] = ran->Get();
    double rand = ran->Get()*m_pdf_over_estimate;
    while(rand > m_partonic.PDF(beam,m_x[beam],m_abs_t,ATOOLS::Flavour(kf_gluon))) {
      m_x[beam] = ran->Get();
      rand = ran->Get()*m_pdf_over_estimate;
    }
  }
  m_p1 = m_Pbeam1*m_x[0];
  m_p2 = m_Pbeam2*m_x[1];

  //Boost to rest frame of the incoming gluons
  Poincare rest    = Poincare(m_p1 + m_p2);
  rest.Boost(m_p1);
  rest.Boost(m_p2);
  m_pl12 = Vec3D(m_p1).Sqr();
  m_pl22 = Vec3D(m_p2).Sqr();
  m_pl1 = sqrt(m_pl12);
  m_pl2 = sqrt(m_pl22);
  double costheta = 1.-m_abs_t/(2.*m_pl12), sintheta = sqrt(1.-sqr(costheta));
  double pt = m_pl1*sintheta, pt2 = sqr(pt);
  double phi(2.*M_PI*ran->Get()), ptx(pt*cos(phi)), pty(pt*sin(phi));
  double pl1(m_sign1*sqrt(m_pl12-pt2)), pl2(-m_sign1*sqrt(m_pl22-pt2));
  double E1 = sqrt(sqr(ptx) + sqr(pty) + sqr(pl1));
  double E2 = sqrt(sqr(-ptx) + sqr(-pty) + sqr(pl2));
  ATOOLS::Vec4D k1(Vec4D(E1,ptx,pty,pl1));
  ATOOLS::Vec4D k2(Vec4D(E2,-ptx,-pty,pl2));
  //Boost back!
  rest.BoostBack(m_p1);
  rest.BoostBack(m_p2);
  rest.BoostBack(k1);
  rest.BoostBack(k2);

  p_emissions->begin()->second.SetMomentum(k1);
  p_emissions->rbegin()->second.SetMomentum(k2);

  T_Prop & prop = *p_props->begin();
  prop.SetQT2(m_qt2);
  prop.SetQ02(m_qt2min);
  prop.SetQ(sqrt(m_qt2)*m_eqt);

  //std::ofstream testfile;
  //testfile.open("./q2file102.txt", std::ios::app);
  //testfile.close();

  //if (m_ana) m_histomap[std::string("Q_elastic")]->Insert(m_abs_t);
  //std::ofstream xfiles;
  //xfiles.open("./xfiles101_noxminmax.txt", std::ios::app);
  //for(size_t beam = 0; beam < 2; beam++) xfiles << m_x[beam] << std::endl;
  //xfiles.close();
  //std::ofstream xfile;
  //xfile.open("./pdf_try_xminmax.txt", std::ios::app);
  //xfile << m_x[0] << "\t" << m_abs_t << "\t" << m_partonic.PDF(0,m_x[0],m_abs_t,ATOOLS::Flavour(kf_gluon)) << std::endl;
  //xfile.close();

  //std::ofstream tvals, pdffile;
  //tvals.open("./tvals.txt", std::ios::app);
  //pdffile.open("./pdf_largerQ2.txt");
  //double delta_x(0.01), delta_q2(0.02);
  //for (double x = 0.; x < 1.; x = x + delta_x) {
  //  for (double q2 = 0.; q2 < 2.5; q2 = q2 + delta_q2) {
  //    pdffile << x << "\t" << q2 << "\t" << m_partonic.PDF(0,x,q2,ATOOLS::Flavour(kf_gluon)) << std::endl;
  //  }
  //}
  //for (int i = 0; i < 1000; i++) {
    //tvals << 0 << "\t" << p_sigma_d->SelectT(0) << endl;
    //tvals << 1 << "\t" << p_sigma_d->SelectT(1) << endl;
    //tvals << 2 << "\t" << p_sigma_d->SelectT(2) << endl;
    //tvals << 4 << "\t" << p_sigma_el->SelectT() << endl;
  //}
  //tvals.close();
  //dffile.close();
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
