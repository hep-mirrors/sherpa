#include"Amegic_Apacic_Interface.H"

#include "XS_Selector.H" 
#include "Data_Read.H"
#include "Message.H"
#include "Random.H"
#include "Running_AlphaS.H"

using namespace SHERPA;
using namespace EXTRAXS;
using namespace AMEGIC;
using namespace APACIC;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;


using namespace EXTRAXS;

// static -- can be used for consitency plots
namespace  SHERPA {
  double amegic_apacic_interface_last_hard_scale=0.;
}

Amegic_Apacic_Interface::Amegic_Apacic_Interface(Matrix_Element_Handler * me,
						 Shower_Handler * shower) :
  Perturbative_Interface(me,shower), p_blob_psme_IS(0), p_blob_psme_FS(0)
{
  p_jf      = 0;
  p_cluster = 0;
  p_xs      = 0;

  p_two2two = new XS_Group(2,2,"Core processes");
  p_fl      = new Flavour[4];
  p_moms    = new Vec4D[4];

  m_ycut    = rpa.gen.Ycut();

  if (rpa.gen.Beam1().IsLepton() && rpa.gen.Beam2().IsLepton()) {
    msg.Debugging()<<" Jet_Finder in Amegic_Apacic_Interface set up  to deal with lepton-lepton collisions "<<endl;
    m_type = 1;
  }
  else if ((!rpa.gen.Beam1().IsLepton() && !rpa.gen.Beam2().IsLepton())) {
    msg.Debugging()<<" Jet_Finder in Amegic_Apacic_Interface set up  to deal with hadron-hadron collisions "<<endl;
    m_type = 4;
  }
  else {
    m_type = 4;
  }

  p_jf       = new ATOOLS::Jet_Finder(m_ycut,m_type);
  p_cluster  = new Cluster_Partons(p_me,p_jf,m_maxjetnumber,p_shower->GetISRHandler()->On(),p_shower->ISROn(),p_shower->FSROn());
  m_jetscale = m_ycut * sqr(rpa.gen.Ecms());
}  

Amegic_Apacic_Interface::~Amegic_Apacic_Interface() 
{
  if (p_two2two) { delete p_two2two; p_two2two = NULL; }
  if (p_jf)      { delete p_jf; p_jf = NULL; }
  if (p_cluster) { delete p_cluster; p_cluster = NULL; }
  // note :
  //  p_shower and p_me are deleted in Jet_Evolution
  //  p_fl an p_moms are deleted in Perturbative_Interface

}


bool Amegic_Apacic_Interface::ClusterConfiguration(Blob * blob)
{
  if (!(p_cluster->ClusterConfiguration(blob,p_me->GetISR_Handler()->X1(),p_me->GetISR_Handler()->X2()))) {
    return 0; // Failure!
  }
  for (int i=0;i<4;i++) {
    p_fl[i]   = p_cluster->Flav(i); 
    p_moms[i] = p_cluster->Momentum(i);
  }

  // prepare Blob , will be inserted later
  if (p_blob_psme_IS) {
    delete p_blob_psme_IS; p_blob_psme_IS = 0;
  }
  if (p_blob_psme_FS) {
    delete p_blob_psme_FS; p_blob_psme_FS = 0;
  }
  if (p_shower->ISROn()) {
    p_blob_psme_IS = new Blob();
    p_blob_psme_IS->SetType(string("ME PS Interface (Sherpa, IS)"));
    p_blob_psme_IS->SetStatus(1);
    for (int i=0;i<blob->NInP();++i) {
      p_blob_psme_IS->AddToOutPartons(blob->InParton(i));
      blob->InParton(i)->SetProductionBlob(p_blob_psme_IS);
      p_blob_psme_IS->SetId(-1);
    }
  }
  if (p_shower->FSROn()) {
    p_blob_psme_FS = new Blob();
    p_blob_psme_FS->SetType(string("ME PS Interface (Sherpa, FS)"));
    p_blob_psme_FS->SetStatus(1);
    for (int i=0;i<blob->NOutP();++i) {
      p_blob_psme_FS->AddToInPartons(blob->OutParton(i));
      blob->OutParton(i)->SetDecayBlob(p_blob_psme_FS);
      p_blob_psme_FS->SetId(-2);
    }
  }

  return 1;  // OK!
}

bool Amegic_Apacic_Interface::DefineInitialConditions(ATOOLS::Blob * blob)
{
  ClusterConfiguration(blob);

  p_xs = 0;

  if (!(XS_Selector::FindInGroup(p_two2two,p_xs,2,2,p_fl))) {
    p_xs = XS_Selector::GetXS(2,2,p_fl);
    if (p_xs) p_two2two->Add(p_xs);
  }

  if (!p_xs) {
    p_cluster->SetColours(p_moms,p_fl);
  }
  else {
    if (!(p_xs->SetColours(p_moms))) return 0;
  }

  if (m_type==1) { // e+ e-
    double sprime = (p_me->Momenta()[0]+p_me->Momenta()[1]).Abs2();
    m_jetscale    = m_ycut * sprime;
  }

  double scale;
  if (p_xs) scale=p_xs->Scale();
  else scale=p_cluster->Scale();
  // save hard scale to be used in plots!
  amegic_apacic_interface_last_hard_scale=scale;

  p_cluster->CalculateWeight(scale,m_jetscale);

  m_weight=p_cluster->Weight();
  if (p_me->Weight()==1. && p_me->UseSudakovWeight()) {
    if (m_weight>ran.Get()) {
      p_cluster->FillTrees(p_shower->GetIniTrees(),p_shower->GetFinTree(),p_xs);

      m_weight=1.;
      return 1;
    }
    m_weight=1.;
  }
  else {
    p_cluster->FillTrees(p_shower->GetIniTrees(),p_shower->GetFinTree(),p_xs);
    return 1;
  }
  return 0;
}

bool   Amegic_Apacic_Interface::FillBlobs(ATOOLS::Blob_List * bl)
{
  if (p_blob_psme_IS) {
    p_blob_psme_IS->SetId(bl->size());
    bl->push_back(p_blob_psme_IS);  
    p_blob_psme_IS=0;
  }
  if (p_blob_psme_FS) {
    p_blob_psme_FS->SetId(bl->size());
    bl->push_back(p_blob_psme_FS);  
    p_blob_psme_FS=0;
  }
  return 1;
}
