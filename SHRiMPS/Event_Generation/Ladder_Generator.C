#include "SHRiMPS/Event_Generation/Ladder_Generator.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Histogram.H"

using namespace SHRIMPS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;




Ladder_Generator::Ladder_Generator():
  m_Y(MBpars.GetEikonalParameters().originalY),
  m_Ymax(MBpars.GetEikonalParameters().Ymax),
  p_ladder(0)
{  }

Ladder_Generator::~Ladder_Generator() {
  if (p_ladder) { delete p_ladder; p_ladder = 0; }
}

void Ladder_Generator::InitCollision(Omega_ik * eikonal,const double & B) {
  Reset();
  p_eikonal = eikonal;
  m_FS.SetEikonal(p_eikonal);
  m_B = B;
}

void Ladder_Generator::Reset() {
  Vec4D cms(rpa->gen.PBeam(0)+rpa->gen.PBeam(1));
  m_E[0] = m_E[1] = cms[0]/2.;
  m_mandels = cms.Abs2();
  m_FS.SetAvailableEnergy(sqrt(m_mandels));
  m_Colours.Reset();
}

bool Ladder_Generator::MakePrimaryLadder(Blob * blob) {
  blob->SetPosition(p_eikonal->SelectB1B2(m_b1,m_b2,m_B));
  m_FS.SetImpactParameters(m_b1,m_b2);
  InitLadder(blob);
  m_FS.FillPrimaryLadder();
  AddInParticles();
  m_Colours(p_ladder);
  if (CheckTotalMomentum()) {
    msg_Out()<<METHOD<<" gives\n"<<(*p_ladder)<<"\n";
    FillBlob(blob);
    return true;
  }
  return false;
}

void Ladder_Generator::InitLadder(Blob * blob) {
  if (p_ladder) { delete p_ladder; p_ladder = 0; }
  p_ladder = new Ladder(blob->Position());
  m_FS.SetLadder(p_ladder);
}

void Ladder_Generator::AddInParticles() {
  Vec4D   mom1(m_FS.Pplus()*Vec4D(1.,0.,0.,1.));
  Vec4D   mom2(m_FS.Pminus()*Vec4D(1.,0.,0.,-1.));
  Flavour flav1(kf_gluon),flav2(kf_gluon);
  Ladder_Particle * in1 = new Ladder_Particle(flav1,mom1,p_ladder->Position());
  Ladder_Particle * in2 = new Ladder_Particle(flav2,mom2,p_ladder->Position());
  p_ladder->SetInParticles(in1,in2);
}

bool Ladder_Generator::CheckTotalMomentum() {
  m_E[0] -= p_ladder->GetIn1()->m_mom[0];
  m_E[1] -= p_ladder->GetIn2()->m_mom[0];
  if (m_E[0]>0. && m_E[1]>0.) return true;
  return false;
}

void Ladder_Generator::FillBlob(Blob * blob) {
  blob->AddToInParticles(p_ladder->GetIn1()->GetParticle());
  blob->AddToInParticles(p_ladder->GetIn2()->GetParticle());
  for (LadderMap::iterator lit=p_ladder->GetEmissions()->begin();
       lit!=p_ladder->GetEmissions()->end();lit++) {
    blob->AddToOutParticles(lit->second.GetParticle());
  }
}


void Ladder_Generator::Test(const std::string & dirname) {
  InitCollision((*MBpars.GetEikonals()->begin()),0.);
  TestPositioning(dirname);
  m_FS.SetImpactParameters(0.,0.);
  m_FS.Test(dirname);
}
  
void Ladder_Generator::TestPositioning(const std::string & dirname) {
  Histogram histo_b1(0,0.,10.,100), histo_b2(0,0.,10.,100);
  for (int i=0;i<100000;i++) {
    p_eikonal->SelectB1B2(m_b1,m_b2,m_B);
    histo_b1.Insert(m_b1);
    histo_b2.Insert(m_b2);
  }
  string name1(string("PosLadders_b1.dat"));
  histo_b1.Finalize();
  histo_b1.Output(dirname+"/"+name1);
  string name2(string("PosLadders_b2.dat"));
  histo_b2.Finalize();
  histo_b2.Output(dirname+"/"+name2);
}
