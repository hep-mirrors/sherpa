#include "Ahadic.H"
#include "Cluster.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Ahadic::Ahadic(string path,string file,bool ana)  :
  m_fullinfo(true), m_maxtrials(3)
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

Return_Value::code Ahadic::Hadronize(ATOOLS::Blob_List * blobs)
{
  cout<<"##########################################################################"<<endl
      <<"##########################################################################"<<endl
      <<"##########################################################################"<<endl
      <<"##########################################################################"<<endl
      <<"##########################################################################"<<endl;

  Blob * blob(NULL);
  Cluster clus;
  clus.ResetClusterNumber();
  for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) {
    if ((*blit)->Has(blob_status::needs_hadronization) &&
	(*blit)->Type()==btp::Fragmentation) {
      blob = (*blit);
      blob->SetType(btp::Cluster_Formation);
      blob->SetTypeSpec("AHADIC-1.0");
      for (int i=0;i<m_maxtrials;i++) {
	switch (int(p_cformhandler->FormClusters(blob))) {
	case int(Return_Value::Retry_Method) :
	  rvalue.IncRetryMethod(METHOD);
	  if (blob) { blob->RemoveInParticles(); blob->RemoveOutParticles(); }
	  continue;
	case int(Return_Value::Success) :
	case int(Return_Value::Warning) :
	  break;
	default:
	  msg.Error()<<"Error in "<<METHOD<<": "<<endl
		     <<"   Unknown return value."<<endl;
	  abort();
	  break;
	}
	//cout<<(*blob)<<endl;
	//return Return_Value::Success;
	
	switch (int(p_cdechandler->DecayClusters(p_cformhandler->GetClusters(),blobs))) {
	case int(Return_Value::Retry_Method) :
	  rvalue.IncRetryMethod(METHOD);
	  if (blob) { blob->RemoveInParticles(); blob->RemoveOutParticles(); }
	  continue;
	case int(Return_Value::Success) :
	case int(Return_Value::Warning) :
	  blob->AddStatus(blob_status::needs_hadrondecays);
	  if (!m_fullinfo) {
	    blobs->MergeSubsequentTypeRecursively(btp::Cluster_Formation,btp::Cluster_Decay);
	    blob->SetStatus(blob_status::needs_hadrondecays);
	    blob->SetType(btp::Fragmentation);
	  }
	  cout<<"##############################################################"<<endl
	      <<"##############################################################"<<endl
	      <<METHOD<<" : "<<clus.RemainingClusters()<<" remaining clusters."<<endl<<(*blobs)
	      <<endl<<"##############################################################"<<endl;

	  return Return_Value::Success;
	default:
	  msg.Error()<<"Error in "<<METHOD<<": "<<endl
		     <<"   Unknown return value."<<endl;
	  abort();
	  break;
	}
      }
    }
  }

  return Return_Value::Retry_Event;
}  
