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
  ana=true;

  p_cformhandler = new Cluster_Formation_Handler(ana);
  p_cdechandler  = new Cluster_Decay_Handler(ana);
  msg_Tracking()<<"Initialisation of Ahadic complete."<<endl;
}

Ahadic::~Ahadic() 
{
  if (p_cdechandler)  { delete p_cdechandler;  p_cdechandler=NULL;  }
  if (p_cformhandler) { delete p_cformhandler; p_cformhandler=NULL; }
}

Return_Value::code Ahadic::Hadronize(ATOOLS::Blob_List * blobs)
{
  if (msg->LevelIsDebugging()) {
    msg_Out()<<"##########################################################################"<<endl
	     <<"###################################### IN ################################"<<endl;
    for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) {
      if ((*blit)->Has(blob_status::needs_hadronization) &&
	  (*blit)->Type()==btp::Fragmentation) 
	msg_Out()<<"##########################################################################"<<endl
		 <<(**blit)<<std::endl
		 <<"##########################################################################"<<endl;
    }
  }
  
  Blob * blob(NULL);
  Cluster clus;
  clus.ResetClusterNumber();

  control::s_AHAparticles=0;

  bool moveon(false);
  for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) {
    if ((*blit)->Has(blob_status::needs_hadronization) &&
	(*blit)->Type()==btp::Fragmentation) {
      blob = (*blit);
      blob->SetType(btp::Cluster_Formation);
      blob->SetTypeSpec("AHADIC-1.0");
      double norm2 = blob->CheckMomentumConservation().Abs2();
#ifdef AHAmomcheck
      msg_Out()<<"=============================================================="<<endl
      	       <<METHOD<<" (0) momentum of incoming in blob : "
      	       <<blob->CheckMomentumConservation()<<"   ("<<norm2<<") "<<endl;
#endif
      for (int i=0;i<m_maxtrials;i++) {
	switch (p_cformhandler->FormClusters(blob,blobs)) {
	case -1 : return Return_Value::Retry_Event;
	case  0 :
	  msg_Error()<<"Warning in "<<METHOD<<" : "<<std::endl
		     <<"   Cluster formation did not work out properly in the "<<i<<"th attempt,"<<std::endl
		     <<"   retry it "<<m_maxtrials<<" times."<<std::endl;
	  rvalue.IncRetryMethod(METHOD);
	  if (blob) { blob->RemoveOutParticles(); }
	  continue;
	case 1 :
	  moveon = true;
	  break;
	}
	if (moveon) break;
      }
#ifdef AHAmomcheck
      if (dabs(blob->CheckMomentumConservation().Abs2())/norm2>1.e-12) {
	msg_Out()<<"=============================================================="<<endl
		 <<METHOD<<" (1) momentum violation for : "<<endl
		 <<"   "<<clus.RemainingClusters()<<" remaining clusters and blobs/particles = "
		 <<control::s_AHAblobs<<"/"<<control::s_AHAparticles<<" vs. "<<blob->NOutP()
		 <<", "<<blob->CheckMomentumConservation()<<"   ("
		 <<blob->CheckMomentumConservation().Abs2()<<") "<<endl
		 <<(*blob)<<endl
		 <<(*p_cformhandler->GetClusters())<<endl;
      }
      else msg_Out()<<"=============================================================="<<endl
      		    <<METHOD<<" (1) momentum conservation for : "<<endl
      		    <<"   "<<clus.RemainingClusters()<<" remaining clusters and blobs/particles = "
      		    <<control::s_AHAblobs<<"/"<<control::s_AHAparticles<<" vs. "<<blob->NOutP()
      		    <<", "<<blob->CheckMomentumConservation().Abs2()<<endl
      		    <<(*blob)<<endl;
#endif

      for (int i=0;i<m_maxtrials;i++) {
	switch (p_cdechandler->DecayClusters(p_cformhandler->GetClusters(),blobs)) {
	case -1 : return Return_Value::Retry_Event;
	case  0 :
	  msg_Error()<<"Warning in "<<METHOD<<" : "<<std::endl
		     <<"   Cluster decays did not work out properly in the "<<i<<"th attempt,"<<std::endl
		     <<"   retry it "<<m_maxtrials<<" times."<<std::endl;
	  rvalue.IncRetryMethod(METHOD);
	  if (blob) { blob->RemoveOutParticles(); }
	  continue;
	case  1 :
	  blob->AddStatus(blob_status::needs_hadrondecays);
	  if (!m_fullinfo) {
	    blobs->MergeSubsequentTypeRecursively(btp::Cluster_Formation,btp::Cluster_Decay,
						  control::s_AHAblobs,
						  control::s_AHAparticles);
	    blob->SetStatus(blob_status::needs_hadrondecays);
	    blob->SetType(btp::Fragmentation);
	  }
	  blob->UnsetStatus(blob_status::needs_hadronization);
	  break;
	}
	if (dabs(blob->CheckMomentumConservation().Abs2())/norm2>1.e-12 ||
	    clus.RemainingClusters()!=1 ||
	    control::s_AHAblobs!=0 || 
	    control::s_AHAparticles!=blob->NOutP() ||
	    control::s_AHAprotoparticles!=0) {
	  msg_Out()<<"ERROR in "<<METHOD<<" : "<<endl
		   <<"   Momentum violation for "<<(clus.RemainingClusters()-1)<<" remaining clusters."<<endl
		   <<" Blobs =  = "<<control::s_AHAblobs<<"/ protos = "<<control::s_AHAprotoparticles
		   <<"/ parts = "<<control::s_AHAparticles<<" vs. "<<blob->NOutP()
		   <<"   : "<<blob->CheckMomentumConservation()<<endl
		   <<(*blob)<<endl;
	  return Return_Value::Retry_Event;
	}
	else {
#ifdef AHAmomcheck
	  msg_Out()<<"Momentum conservation at the end : "<<blob->CheckMomentumConservation()<<endl
		   <<(*blob)<<endl
		   <<"##########################  OUT : No Error ###############################"<<endl
		   <<"##########################################################################"<<endl;
#else
	  if (msg->LevelIsDebugging()) {
	    msg_Out()<<METHOD<<" : "<<(clus.RemainingClusters()-1)<<" remaining clusters."<<endl
		     <<(*blob)<<(*blobs)
		     <<"##############################################################"<<endl
		     <<"##############################################################"<<endl
		     <<"##############################################################"<<endl;
	  }
#endif
	}
      }
    }	  
  }    
  return Return_Value::Success;
}  
