#include "SHERPA/SoftPhysics/Soft_Collision_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "SHERPA/SoftPhysics/Cluster_Algorithm.H"
#include "SINIC++/Main/Sinic.H"

#ifdef PROFILE__all
#define PROFILE__Soft_Collision_Handler
#endif
#ifdef PROFILE__Soft_Collision_Handler
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace SHERPA;
using namespace SINIC;
using namespace ATOOLS;
using namespace std;

Soft_Collision_Handler::Soft_Collision_Handler(string _dir,string _file,
					       BEAM::Beam_Spectra_Handler *beam,
					       PDF::ISR_Handler *isr):
  m_dir(_dir), m_file(_file), m_mode(0), p_sinic(NULL)
{
  Data_Reader dr(" ",";","!","=");
  dr.AddWordSeparator("\t");
  dr.AddIgnore("[");
  dr.AddIgnore("]");
  dr.SetInputPath(m_dir);
  dr.SetInputFile(m_file);
  m_softcollisionmodel=dr.GetValue<string>("SOFT_COLLISIONS",string("Off"));
  if (m_softcollisionmodel==string("Sinic")) {
    m_sfile=dr.GetValue<string>("SINIC_FILE",string("Sinic.dat"));
    p_sinic = new Sinic(&dr,beam,isr);
    m_cluster.SetShowerParams(p_sinic->ShowerMode(),p_sinic->ShowerMinKT2());
    m_mode=1;
    exh->AddTerminatorObject(this);
    return;
  }
  else if (m_softcollisionmodel==string("Off")) return;
  THROW(critical_error,"Soft_Collision model not implemented.");
}
   
Soft_Collision_Handler::~Soft_Collision_Handler() 
{
  if (p_sinic) delete p_sinic;
  exh->RemoveTerminatorObject(this);
}

void Soft_Collision_Handler::CleanUp() {
  if (p_sinic) p_sinic->CleanUp();
} 

void Soft_Collision_Handler::PrepareTerminate() 
{
  std::string path(rpa->gen.Variable("SHERPA_STATUS_PATH"));
  if (path=="") return;
  Copy(m_dir+"/"+m_sfile,path+"/"+m_sfile);
}

ATOOLS::Return_Value::code 
Soft_Collision_Handler::GenerateMinimumBiasEvent(ATOOLS::Blob_List * blobs,
						 double & weight)
{
  PROFILE_HERE;
  int outcome;
  switch (m_mode) {
  case 1: 
    weight = 1.;
    //msg_Out()<<"#################################"<<std::endl
    //	     <<METHOD<<"("<<blobs->size()<<")"<<std::endl;
    outcome = p_sinic->GenerateEvent(blobs);
    //msg_Out()<<(*blobs)<<"\n";
    //msg_Out()<<"####################### yields "<<outcome<<"."<<std::endl
    //	     <<"#################################"<<std::endl<<std::endl;
    switch (outcome) {
    case 1: return Return_Value::Success;
    case 0: return Return_Value::Nothing;
    default:
      msg_Tracking()<<"Error in "<<METHOD<<":"<<std::endl
		    <<"   Did not manage to produce a Sinic event."<<std::endl;
      return Return_Value::New_Event;
    }
  default:
    break;
  }
  return Return_Value::Nothing;
}

Cluster_Amplitude *Soft_Collision_Handler::ClusterConfiguration(Blob *const bl)
{
  m_cluster.SetMinKT2(p_sinic->Smin());
  if (!m_cluster.Cluster(bl)) {
    msg_Error()<<"Error in "<<METHOD<<": could not cluster blob.\n"
	       <<(*bl)<<"\n";
    return NULL;
  }
  Cluster_Amplitude *ampl(m_cluster.Amplitude());
  return ampl;
}


