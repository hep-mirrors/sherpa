#include"Amegic_Apacic_Interface.H"

#include "Data_Read.H"
#include "Message.H"
#include "Random.H"
#include "Running_AlphaS.H"

using namespace SHERPA;
using namespace AMEGIC;
using namespace APACIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace MODEL;
using namespace std;


using namespace EXTRAXS;

Amegic_Apacic_Interface::Amegic_Apacic_Interface(Matrix_Element_Handler * me,
						 Shower_Handler * shower) :
  Perturbative_Interface(me,shower)  //, p_flavs(NULL)
{
  p_jf = 0;
  p_cluster = 0;
  p_xs = 0;
  p_two2two = 0;
  p_xsselector = 0;
  p_blob  = 0;

  p_two2two = new XS_Group(2,2,"Core processes");
  p_fl      = new Flavour[4];
  p_moms    = new Vec4D[4];

  m_ycut    = rpa.integ.Ycut();

  int mode = 4;
  /*
  if (rpa.gen.Beam1().IsLepton() && rpa.gen.Beam2().IsLepton()) {
    msg.Out()<<" Jet_Finder in ME_Interface set up  to deal with lepton-lepton collisions "<<endl;
    mode = 1;
  }
  else if ((!rpa.gen.Beam1().IsLepton() && !rpa.gen.Beam2().IsLepton())) {
    msg.Out()<<" Jet_Finder in ME_Interface set up  to deal with hadron-hadron collisions "<<endl;
    mode = 4;
  }
  else {
    cout<<"ERROR: ME_PS_Interface - DIS is not yet implemented in the Jetfinder "<<endl;
  }
  */

  p_jf       = new APHYTOOLS::Jet_Finder(m_ycut,mode);
  p_cluster  = new Cluster_Partons(p_me,p_jf,m_maxjetnumber);
  m_jetscale = m_ycut * sqr(rpa.gen.Ecms());
}  

Amegic_Apacic_Interface::~Amegic_Apacic_Interface() 
{
  cout<<"in  Amegic_Apacic_Interface::~Amegic_Apacic_Interface()"<<endl;
  if (p_two2two) { delete p_two2two; p_two2two = NULL; }
  if (p_jf)      { delete p_jf; p_jf = NULL; }
  if (p_cluster) { delete p_cluster; p_cluster = NULL; }
  // note
  //  p_shower and p_me are deleted in Jet_Evolution
  cout<<"out  Amegic_Apacic_Interface::~Amegic_Apacic_Interface()"<<endl;

}


bool Amegic_Apacic_Interface::ClusterConfiguration(Blob * blob)
{
  if (!(p_cluster->ClusterConfiguration(blob))) {
    msg.Debugging()<<"Clustering failed !"<<std::endl;
    return 0; // Failure!
  }
  for (int i=0;i<4;i++) {
    p_fl[i]   = p_cluster->Flav(i); 
    p_moms[i] = p_cluster->Momentum(i);
    msg.Tracking()<<" Hard Process : "<<endl;
    msg.Tracking()<<"  "<<p_fl[i]<<" "<<p_moms[i]<<std::endl;
  }
  return 1;  // OK!
}

bool Amegic_Apacic_Interface::DefineInitialConditions(APHYTOOLS::Blob * blob)
{
  ClusterConfiguration(blob);

  p_xs = 0;
  if (!(p_xsselector->FindInGroup(p_two2two,p_xs,2,2,p_fl))) {
    p_xs = p_xsselector->GetXS(2,2,p_fl);
    if (p_xs) p_two2two->Add(p_xs);
  }

  if (!p_xs) {
    msg.Tracking()<<" no xs found "<<endl;
    p_cluster->SetColours(p_moms,p_fl);
  }
  else {
    if (!(p_xs->SetColours(p_moms))) return 0;
  }

  // *AS* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  int type=0; 
  // *AS* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  if (type==1) { // e+ e-
    double sprime = (p_me->Momenta()[0]+p_me->Momenta()[1]).Abs2();
    m_jetscale      = m_ycut * sprime;
  }

  double scale;
  if (p_xs) scale=p_xs->Scale();
  else scale=p_cluster->Scale();
  p_cluster->CalculateWeight(scale,m_jetscale);

  m_weight=p_cluster->Weight();
  if (m_weight>ran.Get()) {
    p_cluster->FillTrees(p_shower->GetIniTrees(),p_shower->GetFinTree(),p_xs);
    return 1;
  }
  msg.Tracking()<<" Reject event due to sudakov weight "<<std::endl;
  return 0;
}


