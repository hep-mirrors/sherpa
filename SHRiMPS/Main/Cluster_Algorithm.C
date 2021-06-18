#include "SHRiMPS/Main/Cluster_Algorithm.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "ATOOLS/Phys/Cluster_Leg.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

#include <limits>

using namespace SHRIMPS;
using namespace ATOOLS;

Cluster_Algorithm::Cluster_Algorithm(const double & Ymax,const double & minkt2):
  p_ampl(NULL), p_clus(NULL),
  p_jf(new JF()), p_jets(new Soft_Jet_Criterion()),
  m_Ymax(Ymax), m_minkt2(minkt2)
{
  m_histomap[std::string("startvspt")] = new Histogram(0,0.0,100.0,100);
  m_histomap[std::string("vetovspt")] = new Histogram(0,0.0,100.0,100);
  m_histomap[std::string("nstartvspt")] = new Histogram(0,0.0,100.0,100);
  m_histomap[std::string("nvetovspt")] = new Histogram(0,0.0,100.0,100);
  p_jf->SetJetCriterion(p_jets);
}

Cluster_Algorithm::~Cluster_Algorithm()
{
  if (p_jf) delete p_jf;
  WriteOutAndDeleteHistograms();
}

int Cluster_Algorithm::ColorConnected(const ColorID &i,const ColorID &j) const
{
  return int(i.m_i==j.m_j && i.m_i!=0)+int(i.m_j==j.m_i && i.m_j!=0);
}

double Cluster_Algorithm::Mass(const Flavour &fl) const
{
  return fl.Mass(true);
}

double Cluster_Algorithm::
PTi2(const ATOOLS::Vec4D & pi,const ATOOLS::Vec4D & pbeam) const
{
  double t((pi+pbeam).Abs2());
  return t*Min(pi[0],pbeam[0])/Max(pi[0],pbeam[0]);
}

double Cluster_Algorithm::
PTij2(const ATOOLS::Vec4D & pi,const ATOOLS::Vec4D & pj) const
{
  double pti2  = pi.PPerp2(), ptj2  = pj.PPerp2();    
  if      (dabs(pi.Y())>m_Ymax && dabs(pj.Y())>m_Ymax) return Max(m_minkt2, pti2);
  else if (dabs(pi.Y())>m_Ymax || dabs(pj.Y())>m_Ymax) return Max(4.*m_minkt2, pti2);
  double ptij2 = Min(pti2,ptj2)*(cosh(pi.Eta()-pj.Eta())-
				 cos(pi.Phi()-pj.Phi()));
  return Max(4.*m_minkt2,ptij2);
}

void Cluster_Algorithm::InitLeg(Cluster_Leg * leg,const double & kt2,
				const size_t & nmaxx) {
  leg->SetNMax(nmaxx);
  leg->SetStat(0);
  leg->SetKT2(0,kt2);
  leg->SetKT2(1,kt2);
}

void Cluster_Algorithm::CreateLegs(Blob * const blob)
{
  size_t nmaxx(blob->NInP()+blob->NOutP()+1);
  for (int i(0);i<blob->NInP();++i) {
    Particle * const part(blob->GetParticle(i));
    p_ampl->CreateLeg(-part->Momentum(),part->Flav().Bar(),
		      ColorID(part->GetFlow(1),part->GetFlow(2)).Conj(),
		      1<<p_ampl->Legs().size());
    InitLeg(p_ampl->Legs().back(),0.,nmaxx);
  }
  for (int i(2);i<blob->NOutP()+2;++i) {
    Particle * const part(blob->GetParticle(i));
    p_ampl->CreateLeg(part->Momentum(),part->Flav(),
		      ColorID(part->GetFlow(1),part->GetFlow(2)),
		      1<<p_ampl->Legs().size());
    InitLeg(p_ampl->Legs().back(),0.,nmaxx);
  }
}

void Cluster_Algorithm::SetAmplitudeProperties(const double & scale) {
  p_ampl->SetNIn(2);
  p_ampl->SetOrderEW(0);
  p_ampl->SetOrderQCD(p_ampl->Legs().size()-2);
  p_ampl->SetMS(this);
  p_ampl->SetKT2(scale);
  p_ampl->SetMuQ2(scale);
  p_ampl->SetMuR2(scale);
  p_ampl->SetMuF2(scale);
  p_ampl->SetMu2(scale);
  p_ampl->SetJF(NULL);
  p_jets->SetClusterAmplitude(p_ampl);
}

double Cluster_Algorithm::SetShowerScales() {
  ClusterLeg_Vector legs(p_ampl->Legs());
  size_t nlegs(legs.size());

  double kt2, kt2test, sij, sijtest, sijmax(m_minkt2);
  int connected;
  for (size_t i=2;i<nlegs;i++) {
    kt2 = 0.;
    sij = 0.;
    Cluster_Leg * split = legs[i];
    for (size_t j=2;j<nlegs;j++) {
      if (i==j) continue;
      Cluster_Leg * spect = legs[j];
      connected = ColorConnected(split->Col(),spect->Col()); 
      if (connected==0) kt2test = m_minkt2;
      else {
	sijtest   = (split->Mom()+spect->Mom()).Abs2();
	if (sijtest>sij) sij = sijtest;
	if (connected==1) kt2test = PTij2(split->Mom(),spect->Mom());
	else kt2test = sijtest;
      }
      //msg_Out()<<"   ["<<i<<" "<<j<<"]("<<connected<<"): kt2test = "<<kt2test<<" "
      //       <<"from ptij = "<<sqrt(PTij2(split->Mom(),spect->Mom()))<<", "
      //       <<"Eij = "<<sqrt(sijtest)<<"\n";
      if (kt2test>kt2) kt2 = kt2test;
    }
    split->SetKT2(0,kt2);
    split->SetKT2(1,kt2);
    if (sij>sijmax) sijmax = sij;
  }
  return Max(4.*m_minkt2,sijmax);
}

bool Cluster_Algorithm::Cluster(Blob *const blob)
{
  p_ampl = Cluster_Amplitude::New(NULL);
  p_jets->Reset();
  CreateLegs(blob);
  double scale = SetShowerScales();
  SetAmplitudeProperties(scale);
  //msg_Out()<<METHOD<<": p_ampl = ["<<p_ampl<<"], jf = ["<<p_jf<<"], "
  //	   <<"jc = ["<<p_jf->JC()<<"].\n"
  //	   <<(*p_ampl)<<"\n";
  //p_jets->Output();
  return true;
}

void Cluster_Algorithm::WriteOutAndDeleteHistograms() {
  if (!m_histomap.empty()) {
    Histogram * histo;
    std::string name;
    for (std::map<std::string,Histogram *>::iterator 
	   hit=m_histomap.begin();hit!=m_histomap.end();hit++) {
      histo = hit->second;
      name  = std::string("Ladder_Analysis/")+hit->first+std::string(".dat");
      histo->Finalize();
      histo->Output(name);
      delete histo;
    }
    m_histomap.clear();
  }
}


