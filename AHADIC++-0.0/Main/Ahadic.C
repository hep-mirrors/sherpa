#include "Ahadic.H"
#include "Cluster.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Ahadic::Ahadic(string path,string file,bool ana) 
{
  hadpars.Init(path,file);

  p_cformhandler = new Cluster_Formation_Handler(ana);
  p_cdechandler  = new Cluster_Decay_Handler(p_cformhandler->GetClusterTransformer(),ana);
  msg.Tracking()<<"Initialisation of Ahadic complete."<<endl;
}

Ahadic::~Ahadic() 
{
  if (p_cdechandler)  { delete p_cdechandler;  p_cdechandler=NULL;  }
  if (p_cformhandler) { delete p_cformhandler; p_cformhandler=NULL; }
}

void Ahadic::Hadronize(ATOOLS::Blob_List * blobs)
{
  msg.Tracking()<<"In Ahadic::Hadronize."<<endl;
  Blob * blob = p_cformhandler->FormClusters(blobs);
  p_cdechandler->DecayClusters(p_cformhandler->GetClusters(),blob);
  msg.Events()<<(*blob)<<endl<<endl;
}
