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
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace MODEL;
using namespace std;


using namespace EXTRAXS;

// static -- can be used for consitency plots
namespace  SHERPA {
  double amegic_apacic_interface_last_hard_scale=0.;
}

Amegic_Apacic_Interface::Amegic_Apacic_Interface(Matrix_Element_Handler * me,
						 Shower_Handler * shower) :
  Perturbative_Interface(me,shower)
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

  p_jf       = new APHYTOOLS::Jet_Finder(m_ycut,m_type);
  p_cluster  = new Cluster_Partons(p_me,p_jf,m_maxjetnumber,p_shower->ISROn(),p_shower->FSROn());
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
  if (!(p_cluster->ClusterConfiguration(blob))) {
    return 0; // Failure!
  }
  for (int i=0;i<4;i++) {
    p_fl[i]   = p_cluster->Flav(i); 
    p_moms[i] = p_cluster->Momentum(i);
  }
  return 1;  // OK!
}

bool Amegic_Apacic_Interface::DefineInitialConditions(APHYTOOLS::Blob * blob)
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
  if (p_me->Weight()==1.) {
    if (m_weight>ran.Get()) {
      if (p_shower->GetIniTrees() || p_shower->GetFinTree())
	p_cluster->FillTrees(p_shower->GetIniTrees(),p_shower->GetFinTree(),p_xs);

      m_weight=1.;
      return 1;
    }
    m_weight=1.;
  }
  else {
    if (p_shower->GetIniTrees() || p_shower->GetFinTree())
      p_cluster->FillTrees(p_shower->GetIniTrees(),p_shower->GetFinTree(),p_xs);
    return 1;
  }
  return 0;
}
