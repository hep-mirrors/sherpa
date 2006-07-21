#include "Ahadic.H"
#include "Cluster.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Ahadic::Ahadic(string path,string file,bool ana) 
{
  hadpars.Init(path,file);
  ana=false;

  p_cformhandler = new Cluster_Formation_Handler(ana);
  p_cdechandler  = new Cluster_Decay_Handler(p_cformhandler->GetClusterTransformer(),ana);
  msg.Tracking()<<"Initialisation of Ahadic complete."<<endl;
}

Ahadic::~Ahadic() 
{
  if (p_cdechandler)  { delete p_cdechandler;  p_cdechandler=NULL;  }
  if (p_cformhandler) { delete p_cformhandler; p_cformhandler=NULL; }
}

bool Ahadic::Hadronize(ATOOLS::Blob_List * blobs,ATOOLS::Particle_List *)
{
  //cout<<"#############################################################"<<endl
  //   <<"#############################################################"<<endl
  //   <<"#############################################################"<<endl
  //    <<METHOD<<"."<<endl;
  Blob * blob = p_cformhandler->FormClusters(blobs);
  p_cdechandler->DecayClusters(p_cformhandler->GetClusters(),blob);
  blob->SetStatus(2);
  blob->SetType(btp::Fragmentation);
  //std::cout<<METHOD<<" : "<<std::endl<<(*blob)<<endl<<endl;
  return true;
}
