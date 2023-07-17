#include <cassert>
#include "AHADIC++/Main/Ahadic.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "AHADIC++/Tools/Cluster.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Shell_Tools.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Ahadic::Ahadic(string shower) :
  m_softclusters(Soft_Cluster_Handler(&m_hadron_list)),
  m_beamparticles(Beam_Particles_Shifter(&m_singlet_list, &m_softclusters)),
  m_sformer(Singlet_Former(&m_singlet_list)),
  m_singletchecker(Singlet_Checker(&m_singlet_list, &m_softclusters)),
  m_gluondecayer(Gluon_Decayer(&m_cluster_list, &m_softclusters)),
  m_clusterdecayer(Cluster_Decayer(&m_cluster_list, &m_softclusters))
{
  ReadMassParameters();
  hadpars = new Hadronisation_Parameters();
  hadpars->Init(shower);
  m_sformer.Init();
  m_beamparticles.Init();
  m_softclusters.Init();
  m_singletchecker.Init();
  m_gluondecayer.Init();
  m_clusterdecayer.Init();
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
      const auto result = Hadronize(blob);
      switch (result) {
      case Return_Value::Success :
	break;
      case Return_Value::Retry_Event :
      case Return_Value::New_Event:
	blobs->ColorConservation();
	msg_Tracking()<<"ERROR in "<<METHOD<<" :\n"
		      <<"   Hadronization for blob "
	  //	      <<"("<<blobs<<"; "
	  //	      <<blob->NInP()<<" -> "<<blob->NOutP()<<") "*/
		      <<"did not work out,";
        if (result==Return_Value::New_Event)
          msg_Tracking()<<" due to momentum problems,";
        msg_Tracking()<<"\n   will trigger "<<result<<":\n"
                      <<(*blobs);
	CleanUp(blob);
	if (result == Return_Value::Retry_Event &&
	    (rpa->gen.Beam1().IsLepton() || rpa->gen.Beam2().IsLepton())) {
	  msg_Tracking()<<METHOD<<": Non-hh collision.\n"
			<<"   Request new event instead.\n";
	  return Return_Value::New_Event;
	}
	return result;
      case Return_Value::Nothing :
      default:
	msg_Tracking()<<"Warning in "<<METHOD<<":\n"
		      <<"   Calling Hadronization for Blob("<<blob<<").\n"
		      <<"   Continue and hope for the best.\n";
      }
    }
    blit++;
  }
  if (m_shrink) Shrink(blobs);
  return Return_Value::Success;
}

Return_Value::code Ahadic::Hadronize(Blob * blob, int retry) {
  Reset();
  m_totmom = blob->CheckMomentumConservation();
  if (!ExtractSinglets(blob) || !ShiftBeamParticles() || !CheckSinglets() ||
      !DecayGluons() ||!DecayClusters()) {
    //msg_Error()<<"ERROR in "<<METHOD<<": Will retry event!\n"
    //	       <<(*blob);
    Reset(blob);
    Reset();
    return Return_Value::New_Event;
  }
  blob->UnsetStatus(blob_status::needs_hadronization);
  blob->SetStatus(blob_status::needs_hadrondecays);
  blob->SetType(btp::Fragmentation);
  blob->SetTypeSpec("AHADIC-1.0");
  FillOutgoingParticles(blob);
  if (dabs(blob->CheckMomentumConservation()[0])>1.e-3) {
    msg_Error()<<"\n"<<METHOD<<" violates four-momentum conservation by "
	       <<blob->CheckMomentumConservation()
	       <<" ("<<blob->CheckMomentumConservation().Abs2()<<")\n";
    Reset(blob);
    return Return_Value::Retry_Event;
  }

  //Ask for weight vector and add to blob
  const auto wgts = m_clusterdecayer.get_variationweights();
  if(wgts.size() != 0) {
    // add weight vector to blob
    msg_Out()<<"wgts.size = " << wgts.size() << std::endl;

    // somehow reset the weight vectors
  }
  return Return_Value::Success;
}

bool Ahadic::ExtractSinglets(Blob * blob)
{
  if (!m_sformer.Extract(blob)) {
    //msg_Error()<<METHOD<<" could not extract singlet.\n";
    return false;
  }
  return true;
}

bool Ahadic::ShiftBeamParticles()
{
  if (!m_beamparticles()) {
    //msg_Error()<<METHOD<<" could not shift beam particles on mass shells.\n";
    return false;
  }
  return true;
}

bool Ahadic::CheckSinglets()
{
  if (!m_singletchecker()) {
    //msg_Error()<<METHOD<<" singlets did not check out.\n";
    return false;
  }
  return true;
}

bool Ahadic::DecayGluons() {
  m_gluondecayer.ResetN();
  while (!m_singlet_list.empty()) {
    if (m_gluondecayer(m_singlet_list.front())) {
      m_singlet_list.pop_front();
    }
    else {
      //msg_Error()<<METHOD<<" could not decay all gluons.\n";
      return false;
    }
  }
  m_gluondecayer.FillNs(m_hadron_list.size());
  return true;
}

bool Ahadic::DecayClusters() {
  bool success = m_clusterdecayer();
  //if (!success) msg_Error()<<METHOD<<" could not decay all clusters.\n";
  return success;
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

void Ahadic::Reset(Blob * blob) {
  m_beamparticles.Reset();
  m_softclusters.Reset();
  m_singletchecker.Reset();
  m_gluondecayer.Reset();
  m_clusterdecayer.Reset();
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
  if (blob) blob->DeleteOutParticles();
  Cluster::Reset();
  Proto_Particle::Reset();
}

bool Ahadic::SanityCheck(Blob * blob,double norm2) {
  Vec4D checkmom(blob->CheckMomentumConservation());
  if (dabs(checkmom.Abs2())/norm2>1.e-12 ||
      (norm2<0. && norm2>0.)) {
    //msg_Error()<<"ERROR in "<<METHOD<<" :\n"
    //	       <<"   Momentum violation in blob: "
    //	       <<checkmom<<" ("<<sqrt(Max(0.,checkmom.Abs2()))<<")\n"
    //	       <<(*blob)<<"\n";
    return false;
  }
  return true;
}

void Ahadic::CleanUp(Blob * blob) {
  Reset();
  if(blob) blob->DeleteOutParticles(0);
}

DEFINE_FRAGMENTATION_GETTER(Ahadic, "Ahadic")
