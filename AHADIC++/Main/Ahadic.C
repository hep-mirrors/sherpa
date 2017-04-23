#include <cassert>
#include "AHADIC++/Main/Ahadic.H"
#include "AHADIC++/Tools/Hadron_Init.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "AHADIC++/Tools/Cluster.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Ahadic::Ahadic(string path, string file) :
  m_beamparticles(Beam_Particles_Shifter(&m_singlet_list)),
  m_sformer(Singlet_Former(&m_singlet_list)),
  m_softclusters(Soft_Cluster_Handler(&m_hadron_list)),
  m_singletchecker(Singlet_Checker(&m_singlet_list, &m_softclusters)),
  m_gluondecayer(Gluon_Decayer(&m_cluster_list, &m_softclusters)),
  m_clusterdecayer(Cluster_Decayer(&m_cluster_list, &m_softclusters))
{  
  hadpars = new Hadronisation_Parameters();
  hadpars->Init(path,file);

  m_beamparticles.Init();
  m_softclusters.Init();
  m_singletchecker.Init();
  m_gluondecayer.Init();
  m_clusterdecayer.Init();
  msg_Tracking()<<"Initialisation of Ahadic complete."<<endl;
}


Ahadic::~Ahadic() 
{
  Reset();
}

Return_Value::code Ahadic::Hadronize(Blob_List * blobs)
{
  static std::string mname(METHOD);
  Return_Value::IncCall(mname);
  for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();) {
    if ((*blit)->Has(blob_status::needs_hadronization) &&
	(*blit)->Type()==btp::Fragmentation) {
      Blob * blob = (*blit);
      switch (Hadronize(blob)) {
      case Return_Value::Success : 
	//msg_Out()<<"****************************************************\n"
	//	 <<"****************************************************\n"
	//	 <<"****************************************************\n"
	//	 <<(*blob)<<"\n";
	break;
      case Return_Value::Retry_Event : 
	blobs->ColorConservation();
	msg_Tracking()<<"ERROR in "<<METHOD<<" :\n"
		      <<"   Hadronization for blob "
		      <<"("<<blob<<"; "
		      <<blob->NInP()<<" -> "<<blob->NOutP()<<") "
		      <<"did not work out,\n"
		      <<"   will trigger retrying the event.\n";
	CleanUp(blob);
	if (rpa->gen.Beam1().IsLepton() ||
	    rpa->gen.Beam2().IsLepton()) {
	  msg_Tracking()<<METHOD<<": Non-hh collision.\n"
			<<"   Request new event instead.\n";
	  return Return_Value::New_Event;
	}
	return Return_Value::Retry_Event;
      case Return_Value::Nothing :
      default:
	msg_Tracking()<<"Warning in "<<METHOD<<":\n"
		      <<"   Calling Hadronization for Blob("<<blob<<").\n"
		      <<"   Continue and hope for the best.\n";
      }
    }
    blit++;      
  }
  return Return_Value::Success;
}  

Return_Value::code Ahadic::Hadronize(Blob * blob, int retry) {
  //msg_Out()<<"********************************************************\n"
  //	   <<"********************************************************\n"
  //	   <<"********************************************************\n"
  //	   <<"In "<<METHOD<<" for \n"<<(*blob)<<"\n";
  Reset();
  if (!ExtractSinglets(blob) || !ShiftBeamParticles() || !CheckSinglets() ||
      !DecayGluons() ||!DecayClusters()) {
    Reset();
    msg_Out()<<"ERROR in "<<METHOD<<" (formation):"<<std::endl
	     <<"   Will retry event."<<std::endl;
    return Return_Value::Retry_Event;
  }
  blob->UnsetStatus(blob_status::needs_hadronization);
  blob->SetStatus(blob_status::needs_hadrondecays);
  blob->SetType(btp::Fragmentation);
  blob->SetTypeSpec("AHADIC-1.0");
  FillOutgoingParticles(blob);
  if (blob->CheckMomentumConservation()[0]>1.e-3) {
    msg_Error()<<METHOD<<" yields "<<blob->CheckMomentumConservation()
	       <<" ("<<blob->CheckMomentumConservation().Abs2()<<")\n"
	       <<(*blob)<<".\n";
    exit(1);
    return Return_Value::Retry_Event;
  }
  //msg_Out()<<"Success.\n"
  //	   <<"********************************************************\n"
  //	   <<"********************************************************\n"
  //	   <<"********************************************************\n";
  return Return_Value::Success;
}
  
bool Ahadic::ExtractSinglets(Blob * blob)
{
  return m_sformer.Extract(blob);
}

bool Ahadic::ShiftBeamParticles()
{
  return m_beamparticles();
}

bool Ahadic::CheckSinglets()
{
  return m_singletchecker();
}

bool Ahadic::DecayGluons() {
  while (!m_singlet_list.empty()) {
    if (m_gluondecayer(m_singlet_list.front()))
      m_singlet_list.pop_front();
    else return false;
  }
  return true;
}

bool Ahadic::DecayClusters() {
  return m_clusterdecayer();
}

void Ahadic::FillOutgoingParticles(Blob * blob) {
  while (!m_hadron_list.empty()) {
    Particle * part = (*m_hadron_list.front())();
    part->SetNumber();
    blob->AddToOutParticles(part);
    delete m_hadron_list.front();
    m_hadron_list.pop_front();
  }
}

void Ahadic::Reset() {
  while (!m_singlet_list.empty()) {
    while (!m_singlet_list.front()->empty()) {
      delete m_singlet_list.front()->front();
      m_singlet_list.front()->pop_front();
    }
    delete m_singlet_list.front();
    m_singlet_list.pop_front();
  }
  m_singlet_list.clear();

  m_cluster_list.clear();
}

bool Ahadic::SanityCheck(Blob * blob,double norm2) {
  Vec4D checkmom(blob->CheckMomentumConservation());
  if (dabs(checkmom.Abs2())/norm2>1.e-12 ||
      (norm2<0. && norm2>0.)) {
    msg_Error()<<"ERROR in "<<METHOD<<" :\n"
	       <<"   Momentum violation in blob: "
	       <<checkmom<<" ("<<sqrt(Max(0.,checkmom.Abs2()))<<")\n"
	       <<(*blob)<<"\n";
    return false;
  }
  return true;
}

void Ahadic::CleanUp(Blob * blob) {
  if(blob) blob->DeleteOutParticles(0);
}
