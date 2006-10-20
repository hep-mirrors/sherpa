#include "Ahadic.H"
#include "Soft_Cluster_Handler.H"
#include "Cluster.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Ahadic::Ahadic(string path,string file,bool ana)  :
  m_fullinfo(false), m_maxtrials(3)
{
  hadpars.Init(path,file);
  ana=false;

  p_cformhandler = new Cluster_Formation_Handler(ana);
  p_cdechandler  = new Cluster_Decay_Handler(p_cformhandler->GetSoftClusterHandler(),ana);
  msg.Tracking()<<"Initialisation of Ahadic complete."<<endl;
}

Ahadic::~Ahadic() 
{
  if (p_cdechandler)  { delete p_cdechandler;  p_cdechandler=NULL;  }
  if (p_cformhandler) { delete p_cformhandler; p_cformhandler=NULL; }
}

Return_Value::code Ahadic::Hadronize(ATOOLS::Blob_List * blobs)
{
  //cout<<"##########################################################################"<<endl
  //   <<"###################################### IN ################################"<<endl;

  Blob * blob(NULL);
  Cluster clus;
  clus.ResetClusterNumber();

  control::s_AHAparticles=0;

  for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) {
    if ((*blit)->Has(blob_status::needs_hadronization) &&
	(*blit)->Type()==btp::Fragmentation) {
      blob = (*blit);
      blob->SetType(btp::Cluster_Formation);
      blob->SetTypeSpec("AHADIC-1.0");
      for (int i=0;i<m_maxtrials;i++) {
	switch (int(p_cformhandler->FormClusters(blob,blobs))) {
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
	/*
	  cout<<METHOD<<" (1) : "<<clus.RemainingClusters()<<" remaining clusters."<<endl
	  <<control::s_AHAblobs<<"/"<<control::s_AHAparticles<<" vs. "<<blob->NOutP()
	  <<"   : "<<blob->CheckMomentumConservation()<<endl
	  <<(*blob)<<endl
	  <<(*p_cformhandler->GetClusters())<<endl;
	*/
	switch (int(p_cdechandler->DecayClusters(p_cformhandler->GetClusters(),blobs))) {
	case int(Return_Value::Retry_Method) :
	  rvalue.IncRetryMethod(METHOD);
	  if (blob) { blob->RemoveInParticles(); blob->RemoveOutParticles(); }
	  continue;
	case int(Return_Value::Success) :
	case int(Return_Value::Warning) :
	  blob->AddStatus(blob_status::needs_hadrondecays);
	  if (!m_fullinfo) {
	    blobs->MergeSubsequentTypeRecursively(btp::Cluster_Formation,btp::Cluster_Decay,
						  control::s_AHAblobs,
						  control::s_AHAparticles);
	    blob->SetStatus(blob_status::needs_hadrondecays);
	    blob->SetType(btp::Fragmentation);
	  }
	  blob->UnsetStatus(blob_status::needs_hadronization);
	  /*
	    cout<<METHOD<<" : "<<clus.RemainingClusters()<<" remaining clusters."<<endl<<(*blob)
	    <<"##############################################################"<<endl
	    <<"##############################################################"<<endl
	    <<"##############################################################"<<endl;
	  */
	  if (dabs(blob->CheckMomentumConservation().Abs2())>1.e-6 ||
	      clus.RemainingClusters()!=1 ||
	      control::s_AHAblobs!=0 || 
	      control::s_AHAparticles!=blob->NOutP() ||
	      control::s_AHAprotoparticles!=0) {
	    cout<<METHOD<<" (2) : "<<clus.RemainingClusters()<<" remaining clusters."<<endl
		<<control::s_AHAblobs<<"/"<<control::s_AHAprotoparticles<<"/"
		<<control::s_AHAparticles<<" vs. "<<blob->NOutP()
		<<"   : "<<blob->CheckMomentumConservation()<<endl
		<<(*blob)<<endl
		<<"################################## OUT ###################################"<<endl
		<<"##########################################################################"<<endl;
	  }
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
