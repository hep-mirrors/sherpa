#include "SHRiMPS/Event_Generation/Ladder_Generator.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHRIMPS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;




Ladder_Generator::Ladder_Generator():
  m_Y(MBpars.GetEikonalParameters().originalY),
  m_Ymax(MBpars.GetEikonalParameters().Ymax),
  p_ladder(0),
  m_Ecms(0.5*(rpa->gen.PBeam(0)+rpa->gen.PBeam(1))[0])
{
  m_histos[string("number")]      = new Histogram(0,0.,10.,10);
  m_histos[string("rapidities")]  = new Histogram(0,-m_Y,m_Y,100);
  m_histos[string("kts_all")]     = new Histogram(0,0.,10.,100);
  m_histos[string("kts_central")] = new Histogram(0,0.,10.,100);
  m_histos[string("kts_forward")] = new Histogram(0,0.,10.,100);
  m_histos[string("kts_beam")]    = new Histogram(0,0.,10.,100);
}

Ladder_Generator::~Ladder_Generator() {
  if (p_ladder) { delete p_ladder; p_ladder = 0; }
  for (map<string,Histogram *>::iterator hit=m_histos.begin();
       hit!=m_histos.end();hit++) {
    hit->second->Finalize();
    hit->second->Output("Tests/"+hit->first+".dat");
    delete hit->second;
  }
}

void Ladder_Generator::
InitCollision(Omega_ik * eikonal,const double & B,const size_t & N) {
  Reset();
  p_eikonal  = eikonal;
  m_B        = B;
  m_Nladders = N;
  m_FS.SetEikonal(p_eikonal);
  m_FS.SetNLadders(N);
}

void Ladder_Generator::Reset() {
  m_E[0] = m_E[1] = m_Ecms;
  m_FS.SetAvailableEnergies(m_E);
  m_Colours.Reset();
}

bool Ladder_Generator::MakePrimaryLadder(Blob * blob,bool isfirst) {
  if (m_E[0]<m_Ecms/10. || m_E[0]<m_Ecms/10.) {
    return false;
  }
  blob->SetPosition(p_eikonal->SelectB1B2(m_b1,m_b2,m_B));
  m_FS.SetImpactParameters(m_b1,m_b2);
  InitLadder(blob);
  bool success;
  // do loop to enforce that first ladder exists.
  do {
    m_FS.FillPrimaryLadder();
    AddInParticles();
  } while (!BreakPrimaryLadderGenerationLoop(isfirst,success));
  if (success) {
    m_Colours(p_ladder);
    UpdateInitialEnergies();
    FillBlob(blob);
    AnalyseLadder();
  }
  return success;
}

void Ladder_Generator::InitLadder(Blob * blob) {
  if (p_ladder) { delete p_ladder; p_ladder = 0; }
  p_ladder = new Ladder(blob->Position());
  m_FS.SetLadder(p_ladder);
}

void Ladder_Generator::AddInParticles() {
  Vec4D   mom1(m_FS.Pplus()/2.*Vec4D(1.,0.,0.,1.));
  Vec4D   mom2(m_FS.Pminus()/2.*Vec4D(1.,0.,0.,-1.));
  Flavour flav1(kf_gluon),flav2(kf_gluon);
  Ladder_Particle * in1 = new Ladder_Particle(flav1,mom1,p_ladder->Position());
  Ladder_Particle * in2 = new Ladder_Particle(flav2,mom2,p_ladder->Position());
  p_ladder->SetInParticles(in1,in2);
}

bool Ladder_Generator::BreakPrimaryLadderGenerationLoop(const bool & isfirst,
							bool & success) {
  success = false;
  // enough energy left in beams?
  bool oktotalmom = CheckTotalMomentum(), nottoohard(false);
  if (!oktotalmom) {
    p_ladder->Reset();
    // if first ladder: try again
    if (isfirst) {
      Reset();
      return false;
    }
    // exit loop
    return true;
  }
  // reweight with mu/t or so
  nottoohard = AcceptLadderForHardness();
  // everything is fine - can exit loop
  if (nottoohard) {
    success = true;
    return true;
  }
  // delete ladder and reset stuff
  p_ladder->Reset();
  // if first ladder: try again
  if (isfirst) Reset();
  return false;
}

bool Ladder_Generator::CheckTotalMomentum() {
  //msg_Out()<<METHOD<<"("<<m_E[0]<<" & "<<m_E[1]<<" vs.\n"
  //	   <<p_ladder->GetIn1()->m_mom<<" & "
  //	   <<p_ladder->GetIn2()->m_mom<<".\n";
  if (m_E[0]-p_ladder->GetIn1()->m_mom[0]>0. &&
      m_E[1]-p_ladder->GetIn2()->m_mom[0]>0.) return true;
  return false;
}

bool Ladder_Generator::AcceptLadderForHardness() {
  p_ladder->ExtractHardest();
  double that   = p_ladder->That(), mu2 = p_ladder->Mu2();
  double expo   = p_ladder->IsHardDiffractive()?1.:1.;
  double weight = pow(4.*mu2/that,expo);
  //msg_Out()<<METHOD<<"(that = "<<that<<", mu2 = "<<mu2<<") "
  //	   <<"--> weight = "<<weight<<"\n";
  if (that<4.*mu2) return true;
  if (weight>ran->Get()) return true;
  //msg_Out()<<(*p_ladder)<<"\n";
  return false;
}

void Ladder_Generator::UpdateInitialEnergies() {
  m_E[0] -= p_ladder->GetIn1()->m_mom[0];
  m_E[1] -= p_ladder->GetIn2()->m_mom[0];
  m_FS.SetAvailableEnergies(m_E);
}

void Ladder_Generator::FillBlob(Blob * blob) {
  blob->AddToInParticles(p_ladder->GetIn1()->GetParticle());
  blob->AddToInParticles(p_ladder->GetIn2()->GetParticle());
  for (LadderMap::iterator lit=p_ladder->GetEmissions()->begin();
       lit!=p_ladder->GetEmissions()->end();lit++) {
    blob->AddToOutParticles(lit->second.GetParticle());
  }
  if (!blob->MomentumConserved()) {
    msg_Error()<<METHOD<<" yields problem with momentum conservation: "
	       <<blob->CheckMomentumConservation()<<"\n"<<(*blob)<<"\n";
    exit(1);
  }
}


void Ladder_Generator::Test(const std::string & dirname) {
  InitCollision((*MBpars.GetEikonals()->begin()),0.,0);
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

void Ladder_Generator::AnalyseLadder() {
  m_histos[string("number")]->Insert(p_ladder->Size());
  for (LadderMap::iterator emit=p_ladder->GetEmissions()->begin();
       emit!=p_ladder->GetEmissions()->end();emit++) {
    m_histos[string("rapidities")]->Insert(emit->first);
    double pt(emit->second.m_mom.PPerp());
    m_histos[string("kts_all")]->Insert(pt);
    if (dabs(emit->first)<2.5)
      m_histos[string("kts_central")]->Insert(pt);
    else if (dabs(emit->first)<5.)
      m_histos[string("kts_forward")]->Insert(pt);
    else
      m_histos[string("kts_beam")]->Insert(pt);
  }
}
