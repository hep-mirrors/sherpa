#include "Ahadic.H"
#include "Cluster.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Ahadic::Ahadic(string path,string file,bool ana)  :
  m_maxtrials(3)
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

Return_Value::code Ahadic::Hadronize(ATOOLS::Blob_List * blobs,ATOOLS::Particle_List *)
{
  //cout<<"#############################################################"<<endl
  //   <<"#############################################################"<<endl
  //   <<"#############################################################"<<endl
  //    <<METHOD<<"."<<endl;
  Blob * blob = new Blob();
  blob->SetType(btp::Cluster_Formation);
  blob->SetId();
  
  for (int i=0;i<m_maxtrials;i++) {
    switch (int(p_cformhandler->FormClusters(blob,blobs))) {
    case int(Return_Value::Retry_Method) :
      rvalue.IncRetryMethod(METHOD);
      if (blob) { blob->RemoveInParticles(); blob->RemoveOutParticles(); }
      continue;
    case int(Return_Value::Success) :
    case int(Return_Value::Warning) :
      //std::cout<<METHOD<<" 1 : Success."<<endl;
      break;
    default:
      msg.Error()<<"Error in "<<METHOD<<": "<<endl
		 <<"   Unknown return value."<<endl;
      abort();
      break;
    }

    switch (int(p_cdechandler->DecayClusters(p_cformhandler->GetClusters(),blob))) {
    case int(Return_Value::Retry_Method) :
      rvalue.IncRetryMethod(METHOD);
      if (blob) { blob->RemoveInParticles(); blob->RemoveOutParticles(); }
      continue;
    case int(Return_Value::Success) :
    case int(Return_Value::Warning) :
      blob->SetStatus(2);
      blob->SetType(btp::Fragmentation);
      std::cout<<METHOD<<" 2 : Success."<<endl;
      return Return_Value::Success;
    default:
      msg.Error()<<"Error in "<<METHOD<<": "<<endl
		 <<"   Unknown return value."<<endl;
      abort();
      break;
    }
  }
  return Return_Value::Retry_Event;
}  
