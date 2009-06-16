#include <cassert>
#include "AHADIC++/Main/Ahadic.H"
#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "AHADIC++/Tools/Cluster.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "AHADIC++/Tools/Hadron_Init.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Ahadic::Ahadic(string path,string file,bool ana)  :
  m_fullinfo(false), m_maxtrials(3), m_clulist(), m_prilist()
{
  
  Data_Reader dr(" ",";","!","=");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(path);
  dr.SetInputFile(file);

  hadpars.Init(path,file);
  ana=false;

  p_cformhandler = new Cluster_Formation_Handler(&m_clulist,&m_prilist,ana);
  p_cdechandler  = new Cluster_Decay_Handler(&m_clulist,ana);
  msg_Tracking()<<"Initialisation of Ahadic complete."<<endl;
}

Ahadic::~Ahadic() 
{
  if (p_cdechandler)  { delete p_cdechandler;  p_cdechandler=NULL;  }
  if (p_cformhandler) { delete p_cformhandler; p_cformhandler=NULL; }
}

Return_Value::code Ahadic::Hadronize(Blob_List * blobs)
{
  static std::string mname(METHOD);
  rvalue.IncCall(mname);
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
  bool moveon(false);
  double norm2(sqr(rpa.gen.Ecms()));
  Return_Value::code result;
  for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();) {
    if ((*blit)->Has(blob_status::needs_hadronization) &&
	(*blit)->Type()==btp::Fragmentation) {
      blob   = (*blit);
      moveon = false;
      Reset();
      for (short int i=0;i<m_maxtrials;i++) {
	try {
	  result = Hadronize(blob,i);
	} catch (Return_Value::code ret) {
	  msg_Error()<<"ERROR in "<<METHOD<<" : "<<std::endl
		     <<"   Hadronization for blob "
		     <<"("<<blob<<"; "<<blob->NInP()<<" -> "<<blob->NOutP()<<") "
		     <<"did not work out,"<<std::endl
		     <<"   will trigger retrying the event."<<std::endl;
	  CleanUp(blob);
	  return Return_Value::Retry_Event;
	}

	switch (result) {
	case Return_Value::Success : 
	  ++blit;
	  moveon = true;
	  break;
	case Return_Value::Retry_Event : 
	  msg_Tracking()<<"ERROR in "<<METHOD<<" : "<<std::endl
			<<"   Hadronization for blob "
			<<"("<<blob<<"; "<<blob->NInP()<<" -> "<<blob->NOutP()<<") "
			<<"did not work out,"<<std::endl
			<<"   will trigger retrying the event."<<std::endl;
	  CleanUp(blob);
	  return result;
	case Return_Value::Retry_Method :
	  msg_Tracking()<<"Warning in "<<METHOD<<" : "<<std::endl
			<<"   Hadronization for blob "
			<<"("<<blob<<"; "<<blob->NInP()<<" -> "<<blob->NOutP()<<") "
			<<"did not work out properly in the "<<(i+1)<<"th attempt,"<<std::endl
			<<"   retry it "<<m_maxtrials<<" times."<<std::endl;
	  rvalue.IncRetryMethod(mname);
	  CleanUp(blob);
	  break;
	case Return_Value::Nothing :
	default:
	  msg_Tracking()<<"Warning in "<<METHOD<<":"<<std::endl
			<<"   Calling Hadronization for Blob("<<blob<<") yields "
			<<int(result)<<"."<<std::endl
			<<"   Continue and hope for the best."<<std::endl;
	  moveon = true;
	  ++blit;
	  break;
	}
	if (moveon) break;
      }

      CleanUp();
      moveon = moveon && SanityCheck(blob,norm2);	    
      if (moveon) {
	blob->SetStatus(blob_status::needs_hadrondecays);
	blob->SetType(btp::Fragmentation);
      }
      else {
	CleanUp(blob);
	return Return_Value::Retry_Event;
      }
    }
    else blit++;
  }
  return Return_Value::Success;
}  



Return_Value::code Ahadic::Hadronize(Blob * blob,int retry) {
  assert(m_clulist.empty() && m_prilist.empty());
  blob->SetType(btp::Cluster_Formation);
  blob->SetTypeSpec("AHADIC-1.0");

  msg_Debugging()<<"In "<<METHOD<<" with "<<std::endl<<(*blob)<<std::endl;
  switch (p_cformhandler->FormClusters(blob)) {
  case -1 : 
    msg_Tracking()<<"ERROR in "<<METHOD<<" :"<<std::endl
		  <<"   Will retry event."<<std::endl;
    p_cformhandler->Reset();
    return Return_Value::Retry_Event;
  case  0 :
    msg_Tracking()<<"ERROR in "<<METHOD<<" :"<<std::endl
		  <<"   Will retry method."<<std::endl;
    p_cformhandler->Reset();
    return Return_Value::Retry_Method;
  case 1 :
    if (retry>0) msg_Out()<<"   Passed cluster formation now ("<<retry<<"th trial)."<<std::endl;    
    break;
  }
  
  msg_Debugging()<<METHOD<<": finally the cluster list :"<<std::endl<<(m_clulist)<<std::endl;
  
  switch (p_cdechandler->DecayClusters(blob)) {
  case -1 : 
    msg_Tracking()<<"ERROR in "<<METHOD<<" :"<<std::endl
		  <<"   Will retry event."<<std::endl;
    return Return_Value::Retry_Event;
  case  0 :
    msg_Tracking()<<"ERROR in "<<METHOD<<" :"<<std::endl
		  <<"   Will retry method."<<std::endl;
    return Return_Value::Retry_Method;
  case  1 :
    if (retry) msg_Out()<<"   Passed cluster decays now ("<<retry<<"th trial)."<<std::endl;
    break;
  }


  if (msg->LevelIsDebugging()) {
    msg_Out()<<"Momentum conservation at the end : "
#ifdef AHAmomcheck
	     <<blob->CheckMomentumConservation()
#endif
	     <<endl<<(*blob)<<endl
	     <<"##########################  OUT : No Error ###############################"<<endl
	     <<"##########################################################################"<<endl;
  }

  assert(m_clulist.empty());
  return Return_Value::Success;

}


void Ahadic::Reset() {
  assert(!Cluster::RemainingClusters());
  Cluster::ResetClusterNumber();
  control::s_AHAparticles=0;
}

bool Ahadic::SanityCheck(Blob * blob,double norm2) {
  if (dabs(blob->CheckMomentumConservation().Abs2())/norm2>1.e-12 ||
      (norm2<0. && norm2>0.) ||
      Cluster::RemainingClusters()!=0 ||
      control::s_AHAparticles!=blob->NOutP()/* ||
					       control::s_AHAprotoparticles!=0*/) {
    msg_Out()<<"ERROR in "<<METHOD<<" : "<<endl
	     <<"   Momentum/particle-blob number violation for "<<(Cluster::RemainingClusters())
	     <<" remaining clusters (norm2 = "<<norm2<<")."<<endl
	     <<"   Protoparticles = "<<control::s_AHAprotoparticles
	     <<"/ parts = "<<control::s_AHAparticles<<" vs. "<<blob->NOutP()
	     <<"   : "<<blob->CheckMomentumConservation()<<endl
	     <<(*blob)<<endl;
    abort();
    return false;
  }
  msg_Debugging()<<"Passed "<<METHOD<<" with "
		 <<"protoparticles = "<<control::s_AHAprotoparticles
		 <<"/ parts = "<<control::s_AHAparticles<<" vs. "<<blob->NOutP()
		 <<"   : "<<blob->CheckMomentumConservation()<<endl
		 <<(*blob)<<endl;
  return true;
}

void Ahadic::CleanUp(Blob * blob) {
  m_clulist.clear();
  while(!m_prilist.empty()) {
    if(msg_LevelIsDebugging()) m_prilist.front()->Print();
    m_prilist.front()->Delete();
    m_prilist.pop_front();
  }
  if(blob) blob->DeleteOutParticles(0);
}
