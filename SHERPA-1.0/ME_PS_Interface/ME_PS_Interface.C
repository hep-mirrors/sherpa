#include "ME_PS_Interface.H"
#include "Run_Parameter.H"
#include "Random.H"

using namespace SHERPA;
using namespace APACIC;
using namespace AMEGIC;
using namespace EXTRAXS;
using namespace ISR;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;



ME_PS_Interface::ME_PS_Interface(ISR_Handler * _isr,int _number) :
  isr(_isr), maxjetnumber(_number)
{
  ps      = new Hard_Interface(isr,maxjetnumber,1);

  two2two = new XS_Group(2,2,"Core processes");
  fl      = new Flavour[4];
  p       = new Vec4D[4];

  ycut      = rpa.integ.Ycut();
  jf        = new APHYTOOLS::Jet_Finder(ycut,1);
  cluster   = new Cluster_Partons(jf,maxjetnumber);
  jetscale  = ycut * sqr(rpa.gen.Ecms());

  ini       = 0;
  fin       = 1;
  NLLweight = 1;
}
  
ME_PS_Interface::~ME_PS_Interface() {}

bool ME_PS_Interface::Treat(Process_Base * proc,Blob * blob,int type) 
{
  if (proc->Nout()>=2) {
    if (!(cluster->ClusterConfiguration(proc,blob))) {
      msg.Debugging()<<"Clustering failed !"<<std::endl;
      return 0;
    }
    for (int i=0;i<4;i++) {
      fl[i]   = cluster->Flav(i); 
      p[i]    = cluster->Momentum(i);
      msg.Debugging()<<"  "<<fl[i]<<" "<<p[i]<<std::endl;
    } 
  }

  xs = 0;
  if (!(xsselector->FindInGroup(two2two,xs,2,2,fl))) {
    xs = xsselector->GetXS(2,2,fl);
    if (xs) two2two->Add(xs,false);
  }

  if (!xs) {
    cout<<" no xs found "<<endl;
    //    cluster->SetColours(p);
    abort();

  }
  else {

    if (!(xs->SetColours(p))) return 0;
  }

  // *AS* 
  //  NLLweight =0;

  if (NLLweight) {
    if (type==1) {
      double sprime = (proc->Momenta()[0]+proc->Momenta()[1]).Abs2();
      jetscale      = ycut * sprime;
    }
    cluster->CalculateWeight(xs->Scale(),jetscale);
    //    cluster->CalculateWeight(sprime,jetscale);
    if (cluster->Weight()>ran.Get()) return 1;
    msg.Debugging()<<" Reject event due to sudakov weight "<<std::endl;
    return 0;
//     msg.Debugging()<<"In principle : Reject event, but for demo : shower !"<<std::endl;
//     return 1;
  }
  return 1;
}

int ME_PS_Interface::PerformShower(Process_Base * proc,int type)
{
  ps->PrepareTrees();
  cluster->FillTrees(ps->IniTrees(),ps->FinTree(),xs);
  int stat=ps->PerformShowers(ini,fin);
  if (!stat) {
    msg.Error()<<"ERROR in ME_PS_Interface::PerformShower."<<std::endl
	       <<"   Parton shower did not work out !!!"<<std::endl;
    return 0;
  }
  return stat;
} 


bool ME_PS_Interface::ExtractPartons(Blob_List * bl,Parton_List * pl) {
  return ps->ExtractPartons(ini,fin,bl,pl);
}

bool ME_PS_Interface::ExtractFinalPartons(Parton_List * pl) {
  return ps->ExtractPartons(0,1,0,pl);
}


