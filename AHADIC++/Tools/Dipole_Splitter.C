#include "AHADIC++/Tools/Dipole_Splitter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AHADIC;
using namespace ATOOLS;



Dipole_Splitter::
Dipole_Splitter(const leading::code & lead,const PTOrder::code & ptorder,
		const ZForm::code & zform,Strong_Coupling * as,
		const bool & analyse) :
  p_tools(new Splitting_Tools(lead,ptorder,zform,as,analyse)),
  m_analyse(analyse)
{ 
  if (m_analyse) {
    m_histograms[std::string("PT_G")]                    = new Histogram(0,0.,25.,50);
    m_histograms[std::string("PT_Q")]                    = new Histogram(0,0.,25.,50);
    m_histograms[std::string("PT_primary")]              = new Histogram(0,0.,25.,50);
    m_histograms[std::string("ZQ")]                      = new Histogram(0,0.,2.,50);
    m_histograms[std::string("ZQ")]                      = new Histogram(0,0.,2.,50);
    m_histograms[std::string("XQ")]                      = new Histogram(0,0.,1.,50);
    m_histograms[std::string("EQ")]                      = new Histogram(0,0.,10.,50);
    m_histograms[std::string("XB_bquark_before_qsplit")] = new Histogram(0,0.,1.,100);
    m_histograms[std::string("XB_bquark_after_qsplit")]  = new Histogram(0,0.,1.,100);
    m_histograms[std::string("XB_bquark_before_gsplit")] = new Histogram(0,0.,1.,100);
    m_histograms[std::string("XB_bquark_after_gsplit")]  = new Histogram(0,0.,1.,100);
    m_histograms[std::string("Masses")]                  = new Histogram(0,0.,50.,100);
  }
}


Dipole_Splitter::~Dipole_Splitter() {
  if (p_tools) { delete p_tools; p_tools=NULL; }
  if (m_analyse) {
    Histogram * histo;
    std::string name;
    for (std::map<std::string,Histogram *>::iterator hit=m_histograms.begin();
	 hit!=m_histograms.end();hit++) {
      histo = hit->second;
      name  = std::string("Fragmentation_Analysis/")+hit->first+std::string(".dat");
      histo->Output(name);
      delete histo;
    }
    m_histograms.clear();
  }
}





bool Dipole_Splitter::SplitCluster(Cluster * cluster) {
  Dipole 
    * dip1(new Dipole(new Proto_Particle((*cluster->GetTrip())),
		      new Proto_Particle((*cluster->GetAnti())))), 
    * dip2;
  bool vetodiquark(cluster->GetTrip()->m_flav.IsDiQuark() ||
		   cluster->GetAnti()->m_flav.IsDiQuark());

  bool first(false);
  if (dip1->Triplet()->m_info=='L' && dip1->AntiTriplet()->m_info=='L' &&
      dip1->Triplet()->m_kt2max==dip1->AntiTriplet()->m_kt2max) first = true;
  if (!EmitGluon(dip1,dip2,first)) return false;

  bool swapped(SelectOrder(dip1,dip2));
  
  if (!SplitDipole(dip1,false,vetodiquark) && 
      !SplitDipole(dip2,false,vetodiquark)) return false;

  Proto_Particle * out1, * out2;
  p_tools->GetNewParticles(out1,out2);

  if (!swapped) {
    dip1->SetAntiTriplet(out1);
    dip2->SetTriplet(out2);
    Produce2Clusters(cluster,dip1,dip2,false);
  }
  else {
    dip2->SetAntiTriplet(out1);
    dip1->SetTriplet(out2);
    Produce2Clusters(cluster,dip2,dip1,true);
  }

  return true;
}

bool Dipole_Splitter::
EmitGluon(Dipole * dip1,Dipole *& dip2,const bool & first) {
  p_tools->SetSpectatorAndSplitter(dip1);

  if (!p_tools->PrepareKinematics(dip1,first) || 
      !p_tools->DetermineSplitting(dip1,first)) {
    p_tools->SwapSpectatorAndSplitter(dip1);
    if (!p_tools->PrepareKinematics(dip1,first) || 
	!p_tools->DetermineSplitting(dip1,first,false)) {
      delete dip1->Triplet();
      delete dip1->AntiTriplet();
      delete dip1;
      return false;
    }
  }
  Proto_Particle * out1, * out2;
  p_tools->AftermathOfSplitting(dip1);
  p_tools->GetNewParticles(out1,out2);
  dip2 = new Dipole(out1,dip1->AntiTriplet());
  dip1->SetAntiTriplet(out1);
  dip1->Update();
  dip2->Update();

  return true;
}

bool Dipole_Splitter::
SplitDipole(Dipole * dip,const bool & first,const bool & vetodiquark) {
  p_tools->SetSpectatorAndSplitter(dip);
  if (!p_tools->PrepareKinematics(dip,first) || 
      !p_tools->DetermineSplitting(dip,first,vetodiquark)) {
    return false;
  }  
  p_tools->AftermathOfSplitting(dip);
  return true;
}

bool Dipole_Splitter::SelectOrder(Dipole *& dip1,Dipole *& dip2) {
  if (ran->Get()<
      (dip1->Mass2()-dip1->Triplet()->m_mom.Abs2())/
      (dip1->Mass2()-dip1->Triplet()->m_mom.Abs2()+
       dip2->Mass2()-dip2->AntiTriplet()->m_mom.Abs2())) {
    Dipole * swap(dip1);
    dip1 = dip2;
    dip2 = swap;
    return true;
  }
  return false;
}

bool Dipole_Splitter::Produce2Clusters(Cluster * cluster,Dipole *& dip1,Dipole *& dip2,
				       const bool & swap) {
  dip1->Update();
  dip2->Update();

  Cluster * clu1 = new Cluster(dip1->Triplet(),dip1->AntiTriplet());
  Cluster * clu2 = new Cluster(dip2->Triplet(),dip2->AntiTriplet());
  if (!swap) {
    cluster->SetLeft(clu1);
    cluster->SetRight(clu2);
  }
  else {
    cluster->SetLeft(clu2);
    cluster->SetRight(clu1);
  }
  clu1->SetPrev(cluster);
  clu2->SetPrev(cluster);

  delete dip1;
  delete dip2;
  if (m_analyse) AnalyseClusterSplitting(cluster);

  return true;
}









void Dipole_Splitter::AnalyseClusterSplitting(Cluster * cluster) {
  Cluster * left(cluster->GetLeft()), * right(cluster->GetRight());
  Histogram * histo;
  histo = (m_histograms.find(std::string("Masses")))->second;
  histo->Insert(left->Mass());
  histo->Insert(right->Mass());
  
  double Eq    = left->GetAnti()->m_mom[0];
  double Eqbar = right->GetTrip()->m_mom[0];
  histo = (m_histograms.find(std::string("EQ")))->second;
  histo->Insert(Eq);
  histo->Insert(Eqbar);
  double Q     = cluster->Mass();
  histo = (m_histograms.find(std::string("XQ")))->second;
  histo->Insert(Eq/Q);
  histo->Insert(Eqbar/Q);
  double zq = 
    (left->GetTrip()->m_mom*left->GetAnti()->m_mom)/
    (left->GetTrip()->m_mom*right->GetAnti()->m_mom);
  double zqbar = 
    (right->GetTrip()->m_mom*right->GetAnti()->m_mom)/
    (left->GetTrip()->m_mom*right->GetAnti()->m_mom);
  histo = (m_histograms.find(std::string("ZQ")))->second;
  histo->Insert(zq);
  histo->Insert(zqbar);
}

void Dipole_Splitter::AnalyseKinematics(const ATOOLS::Vec4D q1,const ATOOLS::Vec4D q2,
					const ATOOLS::Vec4D q3,const bool glusplit) {
  //double Ebeam = rpa->gen.Ecms()/2.;
  Histogram * histo;
  double kt2(p_tools->PT2()), kt(sqrt(kt2)), z(p_tools->Z());
  if (glusplit) {
    histo = (m_histograms.find(std::string("PT_Gluon_Splitting")))->second;
    histo->Insert(kt);
    histo = (m_histograms.find(std::string("PT2_Gluon_Splitting")))->second;
    histo->Insert(kt2);
    histo = (m_histograms.find(std::string("Z_Gluon_Splitting")))->second;
    histo->Insert(z);
  }
  else {
    histo = (m_histograms.find(std::string("PT_Gluon_Emission")))->second;
    histo->Insert(kt);
    histo = (m_histograms.find(std::string("PT2_Gluon_Emission")))->second;
    histo->Insert(kt2);
    histo = (m_histograms.find(std::string("Z_Gluon_Emission")))->second;
    histo->Insert(z);
  }
}



